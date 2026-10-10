function config = pde_experiment_config(method,overrides)
%PDE_EXPERIMENT_CONFIG Final Table 6 parameters and explicit run overrides.
if nargin < 2
    overrides = struct();
end
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
index = find(strcmpi(method,methods),1);
if isempty(index)
    error('PDEExperiment:UnknownMethod','Unknown PDE method: %s',method);
end
method = methods{index};
config = struct('version','pde-calibrated-bootstrap-v2','method',method, ...
    'rt',2,'N',1000,'B',500,'sigma',0.002,'confidence',0.95, ...
    'seed',20261011,'save_results',true,'checkpoint_every',10, ...
    'output_root',fullfile(fileparts(mfilename('fullpath')),'outputs', ...
        'paper-table6'));
options = struct('step',2,'tau',3.25,'max_steps',100000, ...
    'compute_debiased',true);
switch method
    case 'Landweber'
        config.seed = 20261011;
    case 'Showalter'
        config.seed = 20261012;
        options.step = 1;
    case 'HBF'
        config.seed = 20261013;
        options.step = 0.5;
        options.tau = 3.2;
        options.eta = 5;
    case 'ARk'
        config.seed = 20261015;
        options.step = 0.725;
        options.tau = 4.1;
        options.kappa = 1.5;
    case 'SOAR'
        config.seed = 20261009;
        options.step = 0.25;
        options.tau = 3.4;
        options.s = 1.5;
        options.initial_time = 1;
    case 'Nesterov'
        config.seed = 20261010;
        options.step = 1;
        options.tau = 3.3;
        options.omega = 5;
    case 'FAR'
        config.seed = 20261014;
        options.step = 0.5;
        options.tau = 3.35;
        options.theta = 1.25;
end
default_seed = config.seed;
config.options = options;
validateattributes(overrides,{'struct'},{'scalar'});
names = fieldnames(overrides);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(config,name) || ismember(name,{'method','version'})
        error('PDEExperiment:UnknownOverride','Unsupported override: %s',name);
    end
    if strcmp(name,'options')
        option_names = fieldnames(overrides.options);
        for j = 1:numel(option_names)
            field = option_names{j};
            if ~isfield(options,field) || strcmp(field,'compute_debiased')
                error('PDEExperiment:UnknownOverride','Unsupported option: %s',field);
            end
            config.options.(field) = overrides.options.(field);
        end
    else
        config.(name) = overrides.(name);
    end
end
validateattributes(config.rt,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(config.N,{'numeric'},{'scalar','integer','positive'});
validateattributes(config.B,{'numeric'},{'scalar','integer','>=',2});
validateattributes(config.sigma,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(config.confidence,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(config.seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(config.checkpoint_every,{'numeric'},{'scalar','integer','positive'});
validateattributes(config.save_results,{'logical'},{'scalar'});
validateattributes(config.options.max_steps,{'numeric'},{'scalar','integer','positive'});
validateattributes(config.options.step,{'numeric'},{'scalar','finite','positive'});
validateattributes(config.options.tau,{'numeric'},{'scalar','finite','>',1});
config.is_manuscript_configuration = isequal(config.rt,2) && ...
    isequal(config.N,1000) && isequal(config.B,500) && ...
    isequal(config.sigma,0.002) && isequal(config.confidence,0.95) && ...
    isequal(config.options,options) && isequal(config.seed,default_seed);
config.noise_model = 'same-node-noise-on-both-boundaries';
config.stopping_rule = 'realized-effective-noise-W-norm-oracle-TDP';
config.variance_rule = 'weighted-residual-energy-divided-by-trace-QWQ';
config.random_stream = 'mrg32k3a: Cantor(outer-1,inner)*3+branch+1';
config.overrides = overrides;
end
