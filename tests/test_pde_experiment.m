function tests = test_pde_experiment
% Native MATLAB integration tests, using explicit small-run overrides.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
testCase.TestData.root = root;
testCase.TestData.problem = pde_problem(0);
testCase.TestData.smoke = struct('rt',0,'N',1,'B',2,'save_results',false, ...
    'options',struct('max_steps',2));
end

function testFormalDefaultsAndSmokeLabel(testCase)
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
steps = [2,1,0.5,0.725,0.25,1,0.5];
taus = [3.25,3.25,3.2,4.1,3.4,3.3,3.35];
for k = 1:numel(methods)
    config = pde_experiment_config(methods{k});
    verifyEqual(testCase,[config.rt,config.N,config.B],[2,1000,500]);
    verifyEqual(testCase,config.options.step,steps(k));
    verifyEqual(testCase,config.options.tau,taus(k));
    verifyTrue(testCase,config.is_manuscript_configuration);
    smoke = pde_experiment_config(methods{k},testCase.TestData.smoke);
    verifyFalse(testCase,smoke.is_manuscript_configuration);
end
config = pde_experiment_config('ARk');
verifyEqual(testCase,config.options.kappa,1.5);
config = pde_experiment_config('Nesterov');
verifyEqual(testCase,config.options.omega,5);
verifyError(testCase,@() pde_experiment_config('Landweber',struct('typo',1)), ...
    'PDEExperiment:UnknownOverride');
end

function testStreamsAreUniqueStableAndDoNotChangeGlobalRng(testCase)
before = rng;
indices = [];
draws = [];
for outer = 1:3
    for inner = 0:3
        for branch = 0:2
            [draw,index] = pde_noise(123,outer,inner,branch,8);
            indices(end+1) = index; %#ok<AGROW>
            draws(:,end+1) = draw; %#ok<AGROW>
        end
    end
end
verifyEqual(testCase,numel(unique(indices)),numel(indices));
verifyEqual(testCase,size(unique(draws.','rows'),1),size(draws,2));
verifyEqual(testCase,pde_noise(123,2,3,1,8),pde_noise(123,2,3,1,8));
verifyEqual(testCase,rng,before);
end

function testObserveUsesFittedSourceAndEffectiveNoise(testCase)
problem = testCase.TestData.problem;
center = problem.pe*0.4;
stream = RandStream('mrg32k3a','Seed',123);
[gD,gN,u,noise_norm,details] = Observe(0,0.002,center,stream);
verifyEqual(testCase,gD(:,3),problem.F*center,'AbsTol',1e-12);
verifyEqual(testCase,u,problem.solve(problem.M0*center),'AbsTol',1e-12);
verifyEqual(testCase,gD(:,4)-gD(:,3),gN,'AbsTol',1e-16);
verifyEqual(testCase,details.effective_noise,problem.Q*gN,'AbsTol',1e-15);
verifyEqual(testCase,details.z,gD(:,4)-problem.B*gN,'AbsTol',1e-12);
verifyEqual(testCase,noise_norm,sqrt((problem.Q*gN)'*problem.W*(problem.Q*gN)), ...
    'AbsTol',1e-14);
end

function testSevenMethodSmokeAndFittedBootstrapCenters(testCase)
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
problem = testCase.TestData.problem;
for method = methods
    result = run_pde_experiment(method{1},testCase.TestData.smoke);
    verifyEqual(testCase,result.status,'complete');
    verifyEqual(testCase,result.completed_replications,1);
    verifyEqual(testCase,numel(fieldnames(result.summary)),8);
    verifyEqual(testCase,size(result.bootstrap.steps),[1,2,2]);
    verifyTrue(testCase,all(isfinite(cell2mat(struct2cell(result.summary)))));
    config = result.config;
    centers = {result.estimates.original(:,1),result.estimates.debiased(:,1)};
    variances = [result.replications.sigma2_original(1), ...
        result.replications.sigma2_debiased(1)];
    for branch = 1:2
        noise = problem.Q*(sqrt(variances(branch))* ...
            pde_noise(config.seed,1,1,branch,numel(problem.boundary)));
        z = problem.F*centers{branch}+noise;
        options = config.options;
        options.compute_debiased = branch == 2;
        fit = pde_fit(problem,z,method{1},options,sqrt(noise'*problem.W*noise));
        if branch == 1
            expected = max(abs(problem.source_mass_root*(fit.original-centers{branch})));
            observed = result.bootstrap.error_original(1,1);
        else
            expected = max(abs(problem.source_mass_root*(fit.debiased-centers{branch})));
            observed = result.bootstrap.error_debiased(1,1);
        end
        verifyEqual(testCase,observed,expected,'AbsTol',1e-12);
        verifyEqual(testCase,result.bootstrap.steps(1,1,branch),fit.steps);
    end
    error = centers{1}-problem.pe;
    verifyEqual(testCase,result.replications.E2_original,sqrt(error'*problem.M0*error), ...
        'AbsTol',1e-12);
    verifyEqual(testCase,result.summary.mean_E2_original,mean(result.replications.E2_original));
end
end

function testReproducibleRunAndOuterStreams(testCase)
overrides = testCase.TestData.smoke;
overrides.N = 2;
first = run_pde_experiment('Landweber',overrides);
second = run_pde_experiment('Landweber',overrides);
verifyEqual(testCase,first.replications,second.replications);
verifyEqual(testCase,first.bootstrap,second.bootstrap);
verifyNotEqual(testCase,first.bootstrap.substream(1,:,:),first.bootstrap.substream(2,:,:));
verifyTrue(testCase,all(first.replications.cap_hit == 1));
verifyEqual(testCase,first.diagnostics.outer_cap_count,2);
verifyEqual(testCase,first.diagnostics.outer_cap_fraction,1);
verifyEqual(testCase,first.source_error_summary.run_status,'complete');
verifyEqual(testCase,first.diagnostics.bootstrap_cap_fraction, ...
    reshape(mean(mean(first.bootstrap.cap_hit,1),2),1,2));
end

function testBootstrapSelectsItsOwnStoppingSteps(testCase)
overrides = testCase.TestData.smoke;
overrides.B = 4;
overrides.options.max_steps = 2000;
result = run_pde_experiment('Landweber',overrides);
verifyEqual(testCase,result.replications.cap_hit,0);
verifyTrue(testCase,all(result.bootstrap.cap_hit(:) == 0));
verifyTrue(testCase,any(result.bootstrap.steps(:) ~= result.replications.steps));
verifyGreaterThan(testCase,numel(unique(result.bootstrap.steps(:))),1);
end

function testLegacyDriverEntrypointsAreSelfContained(testCase)
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
pde_overrides = testCase.TestData.smoke;
for method = methods
    M0 = -1;
    evalc('run(fullfile(testCase.TestData.root,[''Bootstrap_'' method{1} ''.m'']));');
    verifyEqual(testCase,pde_result.status,'complete');
    verifyEqual(testCase,M0,testCase.TestData.problem.M0);
    verifyEqual(testCase,source_error_summary,pde_result.source_error_summary);
end
end

function testSaveContainsRebuildAndPlotProvenance(testCase)
overrides = testCase.TestData.smoke;
overrides.save_results = true;
overrides.output_root = fullfile(testCase.TestData.root,'outputs','native-tests');
overrides.checkpoint_every = 1;
first = run_pde_experiment('Landweber',overrides);
second = run_pde_experiment('Landweber',overrides);
verifyNotEqual(testCase,first.result_file,second.result_file);
saved = load(first.result_file);
verifyTrue(testCase,all(isfield(saved, ...
    {'result','config','M0','mesh','pe','node0','n','source_error_summary'})));
verifyEqual(testCase,saved.source_error_summary.mass_matrix,saved.M0);
verifyEqual(testCase,saved.source_error_summary.original_estimate,first.estimates.original(:,end));
verifyEqual(testCase,saved.source_error_summary.original_radius,first.replications.quantile_original(end));
verifyEqual(testCase,saved.result.problem_data.d0,testCase.TestData.problem.d0);
verifyTrue(testCase,ismember('run_pde_experiment.m',saved.result.provenance.source_files));
verifyFalse(testCase,isfield(saved.result,'problem'));
end

function testFailureIsRecordedAndNotDropped(testCase)
overrides = testCase.TestData.smoke;
overrides.save_results = true;
test_output_root = fullfile(testCase.TestData.root,'outputs','native-tests');
if ~exist(test_output_root,'dir')
    mkdir(test_output_root);
end
overrides.output_root = tempname(test_output_root);
overrides.options.step = realmax;
verifyError(testCase,@() run_pde_experiment('Landweber',overrides), ...
    'PDEExperiment:ReplicationFailure');
files = dir(fullfile(overrides.output_root,'*','result.mat'));
verifyEqual(testCase,numel(files),1);
saved = load(fullfile(files(1).folder,files(1).name),'result');
verifyEqual(testCase,saved.result.status,'failed');
verifyEqual(testCase,saved.result.completed_replications,0);
verifyTrue(testCase,all(isnan(saved.result.replications.E2_original)));
verifyTrue(testCase,isfield(saved.result,'failure'));
end
