function tests = test_pde_table_export
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
output_root = fullfile(root,'outputs','native-tests');
if ~exist(output_root,'dir')
    mkdir(output_root);
end
directory = tempname(output_root);
mkdir(directory);
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
files = cell(1,7);
for k = 1:7
    result = run_pde_experiment(methods{k},struct('rt',0,'N',2,'B',3, ...
        'options',struct('max_steps',2),'save_results',false));
    files{k} = fullfile(directory,[methods{k} '.mat']);
    save(files{k},'result');
end
testCase.TestData.files = files;
testCase.TestData.directory = directory;
testCase.TestData.methods = methods;
end

function testExplicitTestExportMatchesEightMetricsAndMethodOrder(testCase)
state = warning('off','PDETable:CapHits');
restore = onCleanup(@() warning(state));
e = export_pde_table(testCase.TestData.files(end:-1:1),testCase.TestData.directory,true);
verifyEqual(testCase,e.mode,'sample');
verifyEqual(testCase,e.table.Properties.VariableNames(2:end),testCase.TestData.methods);
verifyEqual(testCase,size(e.table),[8,8]);
verifyTrue(testCase,contains(e.csv_file,'Table6.csv'));
verifyTrue(testCase,contains(e.latex_file,'Table6.tex'));
verifyTrue(testCase,contains(e.diagnostics_file,'Table6_diagnostics'));
verifyTrue(testCase,contains(e.manifest_file,'Table6_manifest'));
for k = 1:7
    saved = load(testCase.TestData.files{k},'result');
    expected = cell2mat(struct2cell(saved.result.summary));
    verifyEqual(testCase,e.table{:,k+1},expected);
    verifyEqual(testCase,e.diagnostics.OuterCapCount(k),sum(saved.result.replications.cap_hit));
end
csv = readtable(e.csv_file);
verifyEqual(testCase,csv{:,2:end},e.table{:,2:end},'AbsTol',1e-14);
latex = fileread(e.latex_file);
verifyTrue(testCase,contains(latex,'Reconstruction errors'));
verifyTrue(testCase,contains(latex,'N=2'));
verifyTrue(testCase,contains(latex,'B=3'));
verifyTrue(testCase,contains(latex,'\mathrm{AR}^{\kappa}'));
verifyTrue(testCase,contains(latex,'95.0\%'));
verifyTrue(testCase,contains(latex,'\times 10^{'));
verifyEqual(testCase,numel(strfind(latex,' \\')),10);
manifest = load(e.manifest_file,'exported');
verifyEqual(testCase,manifest.exported.table,e.table);
verifyEqual(testCase,manifest.exported.source_files,testCase.TestData.files);
end

function testTestRunsRequireExplicitOptInBeforeCreatingOutput(testCase)
output = fullfile(testCase.TestData.directory,'must-not-exist');
verifyError(testCase,@() export_pde_table(testCase.TestData.files,output), ...
    'PDETable:NonManuscript');
verifyFalse(testCase,exist(output,'dir') == 7);
end

function testUniqueOutputDirectoriesAndCapDiagnostics(testCase)
verifyWarning(testCase,@() export_pde_table(testCase.TestData.files, ...
    testCase.TestData.directory,true),'PDETable:CapHits');
state = warning('off','PDETable:CapHits');
restore = onCleanup(@() warning(state));
a = export_pde_table(testCase.TestData.files,testCase.TestData.directory,true);
b = export_pde_table(testCase.TestData.files,testCase.TestData.directory,true);
verifyNotEqual(testCase,a.directory,b.directory);
verifyEqual(testCase,a.table,b.table);
verifyTrue(testCase,all(a.diagnostics.OuterCapFraction >= 0 & a.diagnostics.OuterCapFraction <= 1));
end

function testRequiresAllSevenDistinctMethods(testCase)
files = testCase.TestData.files;
verifyError(testCase,@() export_pde_table(files(1:6),testCase.TestData.directory,true), ...
    'PDETable:MethodSet');
files{2} = files{1};
verifyError(testCase,@() export_pde_table(files,testCase.TestData.directory,true), ...
    'PDETable:MethodSet');
end

function testRejectsOldAndIncompleteResults(testCase)
r = read_first(testCase);
r.version = 'pde-calibrated-bootstrap-v1';
reject(testCase,r,'PDETable:Version');
r = read_first(testCase);
r.source_error_summary.experiment_version = 'pde-calibrated-bootstrap-v1';
reject(testCase,r,'PDETable:Version');
r = read_first(testCase);
r.status = 'failed';
reject(testCase,r,'PDETable:Incomplete');
r = read_first(testCase);
r.completed_replications = r.config.N-1;
reject(testCase,r,'PDETable:Incomplete');
r = read_first(testCase);
r.bootstrap.error_original(1,1) = NaN;
reject(testCase,r,'PDETable:Metrics');
r = read_first(testCase);
r.bootstrap.error_debiased = r.bootstrap.error_debiased(:,1:2);
reject(testCase,r,'PDETable:Metrics');
end

function testRejectsConfigurationMismatchesAndForgedFormalFlag(testCase)
r = read_first(testCase);
r.config.rt = r.config.rt+1;
reject(testCase,r,'PDETable:ConfigurationMismatch');
r = read_first(testCase);
r.figure_data.pe(1) = r.figure_data.pe(1)+1;
reject(testCase,r,'PDETable:ConfigurationMismatch');
r = read_first(testCase);
r.config.is_manuscript_configuration = true;
reject(testCase,r,'PDETable:Configuration');
r = read_first(testCase);
r.config = pde_experiment_config('Landweber');
r.completed_replications = r.config.N;
reject(testCase,r,'PDETable:Metrics');
r = read_first(testCase);
r.config.variance_rule = 'wrong-denominator';
reject(testCase,r,'PDETable:Configuration');
end

function testRejectsInconsistentMetrics(testCase)
r = read_first(testCase);
r.summary.mean_E2_original = r.summary.mean_E2_original+0.1;
reject(testCase,r,'PDETable:Metrics');
r = read_first(testCase);
r.replications.variance_error_debiased(1) = r.replications.variance_error_debiased(1)+0.1;
reject(testCase,r,'PDETable:Metrics');
r = read_first(testCase);
r.replications.quantile_original(1) = r.replications.quantile_original(1)+0.1;
reject(testCase,r,'PDETable:Metrics');
r = read_first(testCase);
r.replications.coverage_debiased(1) = 1-r.replications.coverage_debiased(1);
reject(testCase,r,'PDETable:Metrics');
r = read_first(testCase);
r.bootstrap.steps(1) = r.config.options.max_steps+1;
reject(testCase,r,'PDETable:Metrics');
end

function result = read_first(testCase)
saved = load(testCase.TestData.files{1},'result');
result = saved.result;
end

function reject(testCase,result,identifier)
filename = [tempname(testCase.TestData.directory) '.mat'];
save(filename,'result');
files = testCase.TestData.files;
files{1} = filename;
verifyError(testCase,@() export_pde_table(files,testCase.TestData.directory,true),identifier);
end
