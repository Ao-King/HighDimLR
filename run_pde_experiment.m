function result = run_pde_experiment(method,overrides)
%RUN_PDE_EXPERIMENT Paired Monte Carlo with fitted-center PDE bootstrap.
if nargin < 2
    overrides = struct();
end
config = pde_experiment_config(method,overrides);
problem = pde_problem(config.rt);
n = size(problem.M0,1);
n_boundary = numel(problem.boundary);
result = struct('version',config.version,'config',config,'status','running', ...
    'completed_replications',0,'run_directory','','result_file','');
result.provenance = source_provenance();
result.provenance.matlab_version = version;
result.provenance.started_utc = char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd HH:mm:ss Z'));
result.provenance.boundary_nodes = n_boundary;
result.provenance.d0 = problem.d0;
result.problem_data = struct('CM',problem.CM,'C',problem.C,'M0',problem.M0, ...
    'boundary',problem.boundary,'W',problem.W,'Q',problem.Q, ...
    'F',problem.F,'d0',problem.d0);
names = {'E2_original','E2_debiased','Einf_original','Einf_debiased', ...
    'sigma2_original','sigma2_debiased','variance_error_original', ...
    'variance_error_debiased','quantile_original','quantile_debiased', ...
    'coverage_original','coverage_debiased','steps','cap_hit','noise_norm', ...
    'substream'};
for k = 1:numel(names)
    result.replications.(names{k}) = nan(config.N,1);
end
result.estimates.original = nan(n,config.N);
result.estimates.debiased = nan(n,config.N);
result.bootstrap.error_original = nan(config.N,config.B);
result.bootstrap.error_debiased = nan(config.N,config.B);
result.bootstrap.steps = nan(config.N,config.B,2);
result.bootstrap.cap_hit = nan(config.N,config.B,2);
result.bootstrap.noise_norm = nan(config.N,config.B,2);
result.bootstrap.substream = nan(config.N,config.B,2);
result.source_error_summary = struct();
result.figure_data = struct('M0',problem.M0,'mesh',problem.mesh, ...
    'pe',problem.pe,'node0',problem.node0,'n',n);
if config.save_results
    if ~exist(config.output_root,'dir')
        mkdir(config.output_root);
    end
    result.run_directory = tempname(config.output_root);
    mkdir(result.run_directory);
    result.result_file = fullfile(result.run_directory,'result.mat');
end
clock_start = tic;
outer = 0;
inner = 0;
branch = 0;
try
    for outer = 1:config.N
        inner = 0;
        branch = 0;
        [standard_noise,substream] = pde_noise(config.seed,outer,0,0,n_boundary);
        effective_noise = problem.Q*(config.sigma*standard_noise);
        noise_norm = weighted_norm(problem.W,effective_noise);
        z = problem.F*problem.pe+effective_noise;
        fit = pde_fit(problem,z,config.method,config.options,noise_norm);
        check_fit(fit,true);
        result.estimates.original(:,outer) = fit.original;
        result.estimates.debiased(:,outer) = fit.debiased;
        centers = {fit.original,fit.debiased};
        variances = [fit.sigma2_original,fit.sigma2_debiased];
        for inner = 1:config.B
            for branch = 1:2
                [standard_noise,substream_star] = pde_noise( ...
                    config.seed,outer,inner,branch,n_boundary);
                effective_noise_star = problem.Q*(sqrt(variances(branch))*standard_noise);
                noise_norm_star = weighted_norm(problem.W,effective_noise_star);
                z_star = problem.F*centers{branch}+effective_noise_star;
                options_star = config.options;
                options_star.compute_debiased = branch == 2;
                % Every bootstrap sample selects k*. Its second fit fixes k*.
                fit_star = pde_fit(problem,z_star,config.method,options_star,noise_norm_star);
                check_fit(fit_star,branch == 2);
                if branch == 1
                    difference = fit_star.original-centers{branch};
                    result.bootstrap.error_original(outer,inner) = ...
                        max(abs(problem.source_mass_root*difference));
                else
                    difference = fit_star.debiased-centers{branch};
                    result.bootstrap.error_debiased(outer,inner) = ...
                        max(abs(problem.source_mass_root*difference));
                end
                result.bootstrap.steps(outer,inner,branch) = fit_star.steps;
                result.bootstrap.cap_hit(outer,inner,branch) = fit_star.cap_hit;
                result.bootstrap.noise_norm(outer,inner,branch) = noise_norm_star;
                result.bootstrap.substream(outer,inner,branch) = substream_star;
            end
        end
        original_error = fit.original-problem.pe;
        debiased_error = fit.debiased-problem.pe;
        q_original = quantile(result.bootstrap.error_original(outer,:),config.confidence);
        q_debiased = quantile(result.bootstrap.error_debiased(outer,:),config.confidence);
        row = struct('E2_original',sqrt(max(0,original_error'*problem.M0*original_error)), ...
            'E2_debiased',sqrt(max(0,debiased_error'*problem.M0*debiased_error)), ...
            'Einf_original',max(abs(problem.source_mass_root*original_error)), ...
            'Einf_debiased',max(abs(problem.source_mass_root*debiased_error)), ...
            'sigma2_original',variances(1),'sigma2_debiased',variances(2), ...
            'variance_error_original',abs(variances(1)-config.sigma^2), ...
            'variance_error_debiased',abs(variances(2)-config.sigma^2), ...
            'quantile_original',q_original,'quantile_debiased',q_debiased, ...
            'coverage_original',0,'coverage_debiased',0,'steps',fit.steps, ...
            'cap_hit',fit.cap_hit,'noise_norm',noise_norm,'substream',substream);
        row.coverage_original = row.Einf_original <= q_original;
        row.coverage_debiased = row.Einf_debiased <= q_debiased;
        for k = 1:numel(names)
            value = row.(names{k});
            if ~isfinite(value)
                error('PDEExperiment:NonfiniteMetric','Nonfinite metric: %s',names{k});
            end
            result.replications.(names{k})(outer) = value;
        end
        result.source_error_summary = struct('transform','full-source-mass-sqrt-v1', ...
            'mass_matrix',problem.M0,'original_estimate',fit.original, ...
            'original_radius',q_original,'debiased_estimate',fit.debiased, ...
            'debiased_radius',q_debiased,'experiment_version',config.version, ...
            'method',config.method,'replication',outer,'run_status','running');
        result.completed_replications = outer;
        result.elapsed_seconds = toc(clock_start);
        if config.save_results && mod(outer,config.checkpoint_every) == 0
            persist_result(result);
        end
    end
    result.summary = summarize(result.replications);
    result.diagnostics = struct('outer_cap_count',sum(result.replications.cap_hit), ...
        'outer_cap_fraction',mean(result.replications.cap_hit), ...
        'bootstrap_cap_count',reshape(sum(sum(result.bootstrap.cap_hit,1),2),1,2), ...
        'bootstrap_cap_fraction',reshape(mean(mean(result.bootstrap.cap_hit,1),2),1,2));
    result.status = 'complete';
    result.source_error_summary.run_status = 'complete';
    result.elapsed_seconds = toc(clock_start);
    if config.save_results
        persist_result(result);
    end
catch exception
    result.status = 'failed';
    if isfield(result.source_error_summary,'run_status')
        result.source_error_summary.run_status = 'failed';
    end
    result.elapsed_seconds = toc(clock_start);
    result.failure = struct('outer',outer,'inner',inner,'branch',branch, ...
        'identifier',exception.identifier,'message',exception.message);
    if config.save_results
        persist_result(result);
    end
    contextual = MException('PDEExperiment:ReplicationFailure', ...
        '%s failed at outer=%d, inner=%d, branch=%d. Partial result: %s', ...
        config.method,outer,inner,branch,result.result_file);
    throw(addCause(contextual,exception));
end
end

function value = weighted_norm(W,vector)
value = sqrt(max(0,vector'*W*vector));
if ~isfinite(value)
    error('PDEExperiment:NonfiniteNoise','Effective noise norm is nonfinite.');
end
end

function check_fit(fit,needs_debiased)
values = [fit.original(:);fit.sigma2_original;fit.steps;double(fit.cap_hit)];
if needs_debiased
    values = [values;fit.debiased(:);fit.sigma2_debiased];
end
if any(~isfinite(values)) || fit.sigma2_original < 0 || ...
        (needs_debiased && fit.sigma2_debiased < 0)
    error('PDEExperiment:NonfiniteFit','A reconstruction or variance is invalid.');
end
end

function summary = summarize(rows)
% These are means of norms, not root-mean-square errors.
summary = struct('coverage_original',mean(rows.coverage_original), ...
    'coverage_debiased',mean(rows.coverage_debiased), ...
    'mean_E2_original',mean(rows.E2_original),'mean_E2_debiased',mean(rows.E2_debiased), ...
    'mean_Einf_original',mean(rows.Einf_original),'mean_Einf_debiased',mean(rows.Einf_debiased), ...
    'mean_absolute_variance_error_original',mean(rows.variance_error_original), ...
    'mean_absolute_variance_error_debiased',mean(rows.variance_error_debiased));
end

function provenance = source_provenance()
root = fileparts(mfilename('fullpath'));
files = dir(fullfile(root,'*.m'));
provenance = struct('source_files',{{files.name}},'source_text',{{}});
for k = 1:numel(files)
    provenance.source_text{k} = fileread(fullfile(root,files(k).name));
end
end

function persist_result(result)
% Save plotting inputs together, so an old MAT file cannot supply a new tag.
config = result.config;
source_error_summary = result.source_error_summary;
M0 = result.figure_data.M0;
mesh = result.figure_data.mesh;
pe = result.figure_data.pe;
node0 = result.figure_data.node0;
n = result.figure_data.n;
temporary = fullfile(result.run_directory,'checkpoint.partial.mat');
save(temporary,'result','config','source_error_summary','M0','mesh','pe','node0','n','-v7.3');
movefile(temporary,result.result_file,'f');
end
