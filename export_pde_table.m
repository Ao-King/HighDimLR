function exported = export_pde_table(result_files, output_root, allow_test)
%EXPORT_PDE_TABLE Validate seven saved runs and export the eight Table 6 rows.
% allow_test=true permits smaller runs; parameters remain in the manifest.
if nargin < 2 || isempty(output_root)
    output_root = fullfile(fileparts(mfilename('fullpath')), 'outputs', 'tables');
end
if nargin < 3
    allow_test = false;
end
validateattributes(allow_test, {'logical'}, {'scalar'});
if isstring(result_files)
    result_files = cellstr(result_files);
end
if ~iscell(result_files) || numel(result_files) ~= 7
    error('PDETable:MethodSet', 'Supply exactly seven result.mat paths, one per method.');
end
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
row_fields = {'coverage_original','coverage_debiased','E2_original','E2_debiased', ...
    'Einf_original','Einf_debiased','variance_error_original','variance_error_debiased'};
summary_fields = {'coverage_original','coverage_debiased','mean_E2_original', ...
    'mean_E2_debiased','mean_Einf_original','mean_Einf_debiased', ...
    'mean_absolute_variance_error_original','mean_absolute_variance_error_debiased'};
values = zeros(8,7);
caps = zeros(7,6);
configs = cell(1,7);
sources = cell(1,7);
reference = [];
for index = 1:7
    filename = char(result_files{index});
    [exists, info] = fileattrib(filename);
    if ~exists || info.directory
        error('PDETable:InputFile', 'Result file does not exist: %s', filename);
    end
    saved = load(info.Name, 'result');
    require_fields(saved, {'result'});
    r = saved.result;
    require_fields(r, {'version','config','status','completed_replications', ...
        'replications','bootstrap','summary','problem_data','figure_data', ...
        'source_error_summary'});
    require_fields(r.config, {'method'});
    column = find(strcmp(r.config.method, methods), 1);
    if isempty(column) || ~isempty(configs{column})
        error('PDETable:MethodSet', 'Each recognized method must occur exactly once.');
    end
    defaults = paper_table6_config(methods{column});
    validate_run(r, defaults, allow_test);
    if isempty(reference)
        reference = r;
    else
        compare_runs(reference, r);
    end
    for row = 1:8
        values(row,column) = mean(r.replications.(row_fields{row}));
        require_fields(r.summary, summary_fields(row));
        assert_close(r.summary.(summary_fields{row}), values(row,column), 'summary');
    end
    counts = [sum(r.replications.cap_hit), ...
        reshape(sum(sum(r.bootstrap.cap_hit,1),2),1,2)];
    fractions = counts ./ [r.config.N, r.config.N*r.config.B, r.config.N*r.config.B];
    caps(column,:) = [counts, fractions];
    configs{column} = r.config;
    sources{column} = info.Name;
end

if allow_test
    prefix = 'Table6';
    mode = 'sample';
else
    prefix = 'Table6';
    mode = 'paper';
end
output_root = char(output_root);
if ~exist(output_root, 'dir')
    mkdir(output_root);
end
directory = tempname(output_root);
mkdir(directory);
[~, info] = fileattrib(directory);
directory = info.Name;
table_data = array2table(values, 'VariableNames', methods);
table_data = addvars(table_data, summary_fields.', 'Before', 1, 'NewVariableNames', 'Metric');
diagnostics = array2table(caps, 'VariableNames', {'OuterCapCount', ...
    'OriginalBootstrapCapCount','DebiasedBootstrapCapCount','OuterCapFraction', ...
    'OriginalBootstrapCapFraction','DebiasedBootstrapCapFraction'});
diagnostics = addvars(diagnostics, methods.', 'Before', 1, 'NewVariableNames', 'Method');
exported = struct('version', defaults.version, 'mode', mode, 'directory', directory, ...
    'csv_file', fullfile(directory,[prefix '.csv']), ...
    'latex_file', fullfile(directory,[prefix '.tex']), ...
    'diagnostics_file', fullfile(directory,[prefix '_diagnostics.csv']), ...
    'manifest_file', fullfile(directory,[prefix '_manifest.mat']), ...
    'source_files', {sources}, 'configs', {configs}, 'table', table_data, ...
    'diagnostics', diagnostics, 'coverage_units', 'probability in CSV; percent in LaTeX');
writetable(table_data, exported.csv_file);
writetable(diagnostics, exported.diagnostics_file);
write_latex(exported.latex_file, values, reference.config, configs);
save(exported.manifest_file, 'exported');
if any(caps(:,1:3) > 0, 'all')
    warning('PDETable:CapHits', 'Cap-hit estimates were retained. Inspect %s.', exported.diagnostics_file);
end
end

function validate_run(r, defaults, allow_test)
c = r.config;
require_fields(c, {'version','rt','N','B','sigma','confidence','seed','options', ...
    'is_manuscript_configuration','noise_model','stopping_rule','variance_rule','random_stream'});
if ~strcmp(r.version, defaults.version) || ~strcmp(c.version, defaults.version)
    error('PDETable:Version', 'Only current calibrated PDE v2 results can be exported.');
end
validateattributes(c.N, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(c.B, {'numeric'}, {'scalar','integer','>=',2,'finite'});
validateattributes(c.rt, {'numeric'}, {'scalar','integer','nonnegative','finite'});
validateattributes(c.sigma, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(c.confidence, {'numeric'}, {'scalar','real','>',0,'<',1});
validateattributes(c.seed, {'numeric'}, {'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(c.is_manuscript_configuration, {'logical'}, {'scalar'});
require_fields(c.options, {'max_steps','compute_debiased'});
validateattributes(c.options.max_steps, {'numeric'}, {'scalar','integer','positive','finite'});
if ~isequal(c.options.compute_debiased, true)
    error('PDETable:Configuration', 'Both estimators must have been computed.');
end
actual_formal = isequal(c.rt, defaults.rt) && isequal(c.N, defaults.N) && ...
    isequal(c.B, defaults.B) && isequal(c.sigma, defaults.sigma) && ...
    isequal(c.confidence, defaults.confidence) && isequal(c.options, defaults.options) && ...
    isequal(c.seed, defaults.seed);
if ~isequal(c.is_manuscript_configuration, actual_formal)
    error('PDETable:Configuration', 'The manuscript-configuration flag disagrees with saved parameters.');
end
if ~allow_test && ~actual_formal
    error('PDETable:NonManuscript', 'Nondefault results require explicit allow_test=true.');
end
for field = {'noise_model','stopping_rule','variance_rule','random_stream'}
    if ~isequal(c.(field{1}), defaults.(field{1}))
        error('PDETable:Configuration', 'Unexpected experiment contract: %s.', field{1});
    end
end
if ~strcmp(r.status, 'complete') || ~isequal(r.completed_replications, c.N)
    error('PDETable:Incomplete', 'Every input run must be complete with all N replications.');
end
require_fields(r.source_error_summary, {'transform','experiment_version','run_status','mass_matrix'});
s = r.source_error_summary;
if ~strcmp(s.transform,'full-source-mass-sqrt-v1') || ...
        ~strcmp(s.experiment_version, defaults.version) || ~strcmp(s.run_status,'complete')
    error('PDETable:Version', 'The saved transformation or plotting result is not current and complete.');
end
require_fields(r.problem_data, {'CM','C','M0','boundary','W','Q','F','d0'});
require_fields(r.figure_data, {'M0','mesh','pe','node0','n'});
if ~isequal(s.mass_matrix,r.problem_data.M0) || ~isequal(s.mass_matrix,r.figure_data.M0)
    error('PDETable:Configuration', 'Saved mass matrices disagree.');
end
fields = {'E2_original','E2_debiased','Einf_original','Einf_debiased', ...
    'sigma2_original','sigma2_debiased','variance_error_original','variance_error_debiased', ...
    'quantile_original','quantile_debiased','coverage_original','coverage_debiased', ...
    'steps','cap_hit','noise_norm','substream'};
require_fields(r.replications, fields);
for field = fields
    valid_array(r.replications.(field{1}), [c.N,1], field{1});
end
require_fields(r.bootstrap, {'error_original','error_debiased','steps','cap_hit','noise_norm','substream'});
for field = {'error_original','error_debiased'}
    valid_array(r.bootstrap.(field{1}), [c.N,c.B], field{1});
end
for field = {'steps','cap_hit','noise_norm','substream'}
    valid_array(r.bootstrap.(field{1}), [c.N,c.B,2], field{1});
end
for suffix = {'original','debiased'}
    name = suffix{1};
    rows = r.replications;
    assert_close(rows.(['variance_error_' name]), abs(rows.(['sigma2_' name])-c.sigma^2), 'variance error');
    radius = quantile(r.bootstrap.(['error_' name]), c.confidence, 2);
    assert_close(rows.(['quantile_' name]), radius, 'bootstrap quantile');
    if ~isequal(rows.(['coverage_' name]), double(rows.(['Einf_' name]) <= radius))
        error('PDETable:Metrics', 'Coverage does not match the closed confidence region.');
    end
end
counts = [r.replications.steps(:);r.bootstrap.steps(:)];
flags = [r.replications.cap_hit(:);r.bootstrap.cap_hit(:)];
if any(counts ~= floor(counts) | counts > c.options.max_steps) || any(flags ~= 0 & flags ~= 1)
    error('PDETable:Metrics', 'Invalid iteration counts or cap-hit flags.');
end
end

function compare_runs(a, b)
for field = {'rt','N','B','sigma','confidence','noise_model', ...
        'stopping_rule','variance_rule','random_stream'}
    if ~isequal(a.config.(field{1}), b.config.(field{1}))
        error('PDETable:ConfigurationMismatch', 'The seven runs differ in %s.', field{1});
    end
end
if ~isequal(a.config.options.max_steps,b.config.options.max_steps) || ...
        ~isequal(a.problem_data,b.problem_data) || ~isequal(a.figure_data,b.figure_data)
    error('PDETable:ConfigurationMismatch', 'Iteration caps, operators, meshes, and true sources must match.');
end
end

function require_fields(value, fields)
if ~isstruct(value) || ~isscalar(value) || ~all(isfield(value,fields))
    error('PDETable:Schema', 'A saved result is missing required fields.');
end
end

function valid_array(value, dimensions, name)
if ~isnumeric(value) || ~isreal(value) || ~isequal(size(value),dimensions) || ...
        any(~isfinite(value(:))) || any(value(:) < 0)
    error('PDETable:Metrics', 'Invalid or incomplete numeric array: %s.', name);
end
end

function assert_close(actual, expected, name)
if ~isnumeric(actual) || ~isreal(actual) || ~isequal(size(actual),size(expected)) || ...
        any(~isfinite(actual(:))) || any(abs(actual(:)-expected(:)) > 64*eps(max(1,abs(expected(:)))))
    error('PDETable:Metrics', 'Saved %s disagrees with the underlying replications.', name);
end
end

function write_latex(filename, values, config, configs)
labels = {'Coverage I','Coverage II','$\overline{E_2(\hat f)}$', ...
    '$\overline{E_2(\tilde f)}$','$\overline{E_\infty(\hat f)}$', ...
    '$\overline{E_\infty(\tilde f)}$','$\overline{|\widehat\sigma^2-\sigma^2|}$', ...
    '$\overline{|\widetilde\sigma^2-\sigma^2|}$'};
fid = fopen(filename, 'w');
if fid == -1
    error('PDETable:OutputFile', 'Cannot create LaTeX output: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%% Generated from seven validated pde-calibrated-bootstrap-v2 runs.\n');
caption = 'Reconstruction errors, variance-estimation errors, and empirical coverage in the inverse source experiment.';
fprintf(fid, '\\begin{table*}[t]\n\\centering\n\\footnotesize\n');
fprintf(fid, '\\caption{%s $N=%d$, $B=%d$, nominal coverage $%.1f\\%%$.}\n', ...
    caption, config.N, config.B, 100*config.confidence);
fprintf(fid, '\\begin{tabular}{lrrrrrrr}\n\\hline\n');
fprintf(fid, ' & Landweber & Showalter & HBF & $\\mathrm{AR}^{\\kappa}$ & SOAR & Nesterov & FAR \\\\\n\\hline\n');
fprintf(fid,'$\\varsigma$');
for column = 1:7
    fprintf(fid,' & %.4g',configs{column}.options.tau);
end
fprintf(fid,' \\\\\n');
for row = 1:8
    fprintf(fid, '%s', labels{row});
    for column = 1:7
        value = values(row,column);
        if row <= 2
            formatted = sprintf('%.1f\\%%',100*value);
        elseif row <= 6
            formatted = sprintf('%.4f',value);
        elseif value == 0
            formatted = '$0$';
        else
            exponent = floor(log10(value));
            formatted = sprintf('$%.3f\\times 10^{%d}$', value/10^exponent, exponent);
        end
        fprintf(fid, ' & %s', formatted);
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\hline\n\\end{tabular}\n\\end{table*}\n');
end
