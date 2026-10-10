function tests = test_paper_table6
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'benchmarks','pde_gpu_preflight'), ...
    fullfile(root,'benchmarks','pde_formal_cpu'), ...
    fullfile(root,'benchmarks','pde_far_formal'));
testCase.TestData.root = root;
end

function testPublishedParametersAndSeeds(testCase)
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
taus = [3.25,3.25,3.2,4.1,3.4,3.3,3.35];
steps = [2,1,.5,.725,.25,1,.5];
seeds = [20261011,20261012,20261013,20261015,20261009,20261010,20261014];
for index = 1:numel(methods)
    config = paper_table6_config(methods{index});
    verifyEqual(testCase,[config.rt config.N config.B],[2 1000 500]);
    verifyEqual(testCase,[config.sigma config.confidence],[.002 .95]);
    verifyEqual(testCase,config.options.tau,taus(index));
    verifyEqual(testCase,config.options.step,steps(index));
    verifyEqual(testCase,config.options.max_steps,100000);
    verifyEqual(testCase,config.seed,seeds(index));
    verifyTrue(testCase,config.is_manuscript_configuration);
end
end

function testSmallRunUsesFormalBatchedReplications(testCase)
overrides = struct('rt',0,'N',1,'B',4,'save_results',false, ...
    'options',struct('max_steps',5));
state = warning('off','PDEPaper:IterationCap');
restore = onCleanup(@()warning(state)); %#ok<NASGU>
result = run_paper_table6('all',overrides);
problem = pde_problem(0);
methods = fieldnames(result);
for index = 1:numel(methods)
    run = result.(methods{index});
    verifyEqual(testCase,run.status,'complete');
    verifyEqual(testCase,run.completed_replications,1);
    verifyFalse(testCase,run.config.is_manuscript_configuration);
    if strcmp(methods{index},'FAR')
        row = pde_far_replication(problem,run.config,1,struct('backend','far-spectral','payload',[]));
    else
        row = pde_formal_replication(problem,run.config,1);
    end
    names = fieldnames(row.metrics);
    for key = 1:numel(names)
        verifyEqual(testCase,run.replications.(names{key}),double(row.metrics.(names{key})));
    end
    verifyEqual(testCase,run.estimates.original,row.original);
    verifyEqual(testCase,run.estimates.debiased,row.debiased);
    verifyEqual(testCase,run.bootstrap.error_original,row.bootstrap.error_original.');
    verifyEqual(testCase,run.bootstrap.error_debiased,row.bootstrap.error_debiased.');
    verifyEqual(testCase,reshape(run.bootstrap.steps,4,2),row.bootstrap.steps);
end
end

function testRunSavesReplicationsAndRebuildInputs(testCase)
directory = tempname;
overrides = struct('rt',0,'N',1,'B',2,'output_root',directory, ...
    'options',struct('max_steps',5));
state = warning('off','PDEPaper:IterationCap');
restore = onCleanup(@()warning(state)); %#ok<NASGU>
result = run_paper_table6('Landweber',overrides);
saved = load(result.result_file,'result','config','M0','mesh','source_error_summary');
verifyEqual(testCase,saved.result.replications,result.replications);
verifyEqual(testCase,saved.source_error_summary.run_status,'complete');
verifyTrue(testCase,exist(fullfile(result.run_directory,'problem.mat'),'file') == 2);
row = load(fullfile(result.run_directory,'replications','outer_0001.mat'),'row');
verifyEqual(testCase,row.row.original,result.estimates.original);
verifyEqual(testCase,row.row.config.seed,20261011);
end
