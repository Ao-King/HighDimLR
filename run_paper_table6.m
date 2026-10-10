function result = run_paper_table6(method,overrides)
%RUN_PAPER_TABLE6 Reproduce one method, or all seven methods, in Table 6.
%   run_paper_table6('SOAR') runs N=1000 and B=500 with the published settings.
%   run_paper_table6('all',struct('N',1,'B',8)) runs a small paired experiment.
%   Each bootstrap fit chooses its own discrepancy stopping step. Debiasing
%   replays that same step count. FAR uses the validated full-mode filters.
if nargin < 1 || isempty(method), method = 'all'; end
if nargin < 2, overrides = struct(); end
root = fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'benchmarks','pde_gpu_preflight'), ...
    fullfile(root,'benchmarks','pde_formal_cpu'), ...
    fullfile(root,'benchmarks','pde_far_formal'));
if strcmpi(method,'all')
    methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
    result = struct();
    for index = 1:numel(methods)
        result.(methods{index}) = run_paper_table6(methods{index},overrides);
    end
    return
end
config = paper_table6_config(method,overrides);
problem = pde_problem(config.rt);
n = size(problem.M0,1);
result = struct('version',config.version,'config',config,'status','running', ...
    'completed_replications',0,'run_directory','','result_file','');
result.provenance = struct('matlab_version',version,'execution_backend', ...
    config.execution_backend,'started_utc',utc_now(),'boundary_nodes', ...
    numel(problem.boundary),'d0',problem.d0);
result.problem_data = struct('CM',problem.CM,'C',problem.C,'M0',problem.M0, ...
    'boundary',problem.boundary,'W',problem.W,'Q',problem.Q,'F',problem.F,'d0',problem.d0);
result.figure_data = struct('M0',problem.M0,'mesh',problem.mesh, ...
    'pe',problem.pe,'node0',problem.node0,'n',n);
metric_fields = {'E2_original','E2_debiased','Einf_original','Einf_debiased', ...
    'sigma2_original','sigma2_debiased','variance_error_original', ...
    'variance_error_debiased','quantile_original','quantile_debiased', ...
    'coverage_original','coverage_debiased','steps','cap_hit','noise_norm','substream'};
for index = 1:numel(metric_fields)
    result.replications.(metric_fields{index}) = nan(config.N,1);
end
result.estimates.original = nan(n,config.N);
result.estimates.debiased = nan(n,config.N);
result.bootstrap.error_original = nan(config.N,config.B);
result.bootstrap.error_debiased = nan(config.N,config.B);
for key = {'steps','cap_hit','noise_norm','substream'}
    result.bootstrap.(key{1}) = nan(config.N,config.B,2);
end
result.source_error_summary = struct();
if config.save_results
    if ~exist(config.output_root,'dir'), mkdir(config.output_root); end
    result.run_directory = tempname(config.output_root);
    mkdir(result.run_directory);
    mkdir(fullfile(result.run_directory,'replications'));
    result.result_file = fullfile(result.run_directory,'result.mat');
    problem_data = result.problem_data; %#ok<NASGU>
    save(fullfile(result.run_directory,'problem.mat'),'problem_data','-v7.3');
end
clock_start = tic;
cache = struct('backend',config.far_backend,'payload',[]);
outer = 0;
try
    for outer = 1:config.N
        if strcmp(config.method,'FAR')
            [row,cache] = pde_far_replication(problem,config,outer,cache);
        else
            row = pde_formal_replication(problem,config,outer);
        end
        for index = 1:numel(metric_fields)
            name = metric_fields{index};
            result.replications.(name)(outer) = row.metrics.(name);
        end
        result.estimates.original(:,outer) = row.original;
        result.estimates.debiased(:,outer) = row.debiased;
        result.bootstrap.error_original(outer,:) = row.bootstrap.error_original.';
        result.bootstrap.error_debiased(outer,:) = row.bootstrap.error_debiased.';
        for key = {'steps','cap_hit','noise_norm','substream'}
            result.bootstrap.(key{1})(outer,:,:) = reshape(row.bootstrap.(key{1}),1,config.B,2);
        end
        result.source_error_summary = struct('transform','full-source-mass-sqrt-v1', ...
            'mass_matrix',problem.M0,'original_estimate',row.original, ...
            'original_radius',row.metrics.quantile_original,'debiased_estimate',row.debiased, ...
            'debiased_radius',row.metrics.quantile_debiased,'experiment_version',config.version, ...
            'method',config.method,'replication',outer,'run_status','running');
        result.completed_replications = outer;
        result.elapsed_seconds = toc(clock_start);
        if config.save_results
            row.config = config;
            save(fullfile(result.run_directory,'replications',sprintf('outer_%04d.mat',outer)), ...
                'row','-v7.3');
            if mod(outer,config.checkpoint_every) == 0, persist_result(result); end
        end
        if row.metrics.cap_hit || any(row.bootstrap.cap_hit,'all')
            warning('PDEPaper:IterationCap', ...
                '%s replication %d contains fits at the iteration cap; all samples were retained.', ...
                config.method,outer);
        end
        fprintf('%s: %d/%d replications, %.1f seconds\n', ...
            config.method,outer,config.N,result.elapsed_seconds);
    end
    result.summary = summarize(result.replications);
    result.diagnostics = struct('outer_cap_count',sum(result.replications.cap_hit), ...
        'outer_cap_fraction',mean(result.replications.cap_hit), ...
        'bootstrap_cap_count',reshape(sum(sum(result.bootstrap.cap_hit,1),2),1,2), ...
        'bootstrap_cap_fraction',reshape(mean(mean(result.bootstrap.cap_hit,1),2),1,2));
    result.status = 'complete';
    result.source_error_summary.run_status = 'complete';
    result.provenance.finished_utc = utc_now();
    result.elapsed_seconds = toc(clock_start);
    if config.save_results, persist_result(result); end
catch exception
    result.status = 'failed';
    if isfield(result.source_error_summary,'run_status')
        result.source_error_summary.run_status = 'failed';
    end
    result.elapsed_seconds = toc(clock_start);
    result.failure = struct('outer',outer,'identifier',exception.identifier,'message',exception.message);
    if config.save_results, persist_result(result); end
    contextual = MException('PDEPaper:ReplicationFailure', ...
        '%s failed at outer=%d. Partial result: %s',config.method,outer,result.result_file);
    throw(addCause(contextual,exception));
end
end

function summary = summarize(rows)
summary = struct('coverage_original',mean(rows.coverage_original), ...
    'coverage_debiased',mean(rows.coverage_debiased), ...
    'mean_E2_original',mean(rows.E2_original),'mean_E2_debiased',mean(rows.E2_debiased), ...
    'mean_Einf_original',mean(rows.Einf_original),'mean_Einf_debiased',mean(rows.Einf_debiased), ...
    'mean_absolute_variance_error_original',mean(rows.variance_error_original), ...
    'mean_absolute_variance_error_debiased',mean(rows.variance_error_debiased));
end

function value = utc_now()
value = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss Z'));
end

function persist_result(result)
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
