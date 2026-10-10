function row = pde_formal_replication(problem,config,outer)
%PDE_FORMAL_REPLICATION One CPU replication using unchanged batched kernels.
% Uses every configured bootstrap sample in each of the two fitted-center
% branches. Random substreams, per-sample stops and fixed-k debias replay are
% inherited from pde_noise and the validated pde_fit_batched implementation.
% The caller supplies the assembled problem and experiment configuration.
% No result files, GPU arrays, worker pools or algorithm overrides are created.
clock_start = tic;
validateattributes(config,{'struct'},{'scalar'});
validateattributes(outer,{'numeric'},{'scalar','integer','finite','positive','<=',config.N});
validateattributes(config.B,{'numeric'},{'scalar','integer','finite','>=',2});
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov'};
if ~any(strcmpi(config.method,methods))
    error('PDEFormal:Method','This formal CPU path supports the six non-FAR methods only.');
end
if ~config.options.compute_debiased
    error('PDEFormal:Debias','Both original and debiased estimates are required.');
end
cpu_fields = {'F','G','W','Q','pe','M0','source_mass_root','d0'};
if isfield(problem,'normal_boundary'), cpu_fields{end+1} = 'normal_boundary'; end
for j = 1:numel(cpu_fields)
    value = problem.(cpu_fields{j});
    if isa(value,'gpuArray') || ~isa(value,'double') || ~isreal(value)
        error('PDEFormal:Precision','problem.%s must be a real CPU double array.',cpu_fields{j});
    end
end
b = size(problem.F,1);
B = config.B;
[standard_noise,substream] = pde_noise(config.seed,outer,0,0,b);
effective_noise = problem.Q*(config.sigma*standard_noise);
noise_norm = weighted_norm(problem.W,effective_noise);
z = problem.F*problem.pe+effective_noise;
fit = pde_fit_batched(problem,z,config.method,config.options,noise_norm);
check_fit(fit,true);
centers = {fit.original,fit.debiased};
variances = [fit.sigma2_original,fit.sigma2_debiased];
bootstrap = struct('error_original',nan(B,1),'error_debiased',nan(B,1), ...
    'steps',nan(B,2),'cap_hit',false(B,2),'noise_norm',nan(B,2), ...
    'substream',nan(B,2));
for branch = 1:2
    observations = zeros(b,B);
    deltas = zeros(1,B);
    fitted_center = problem.F*centers{branch};
    for inner = 1:B
        [standard_noise,substream_star] = pde_noise(config.seed,outer,inner,branch,b);
        effective_noise_star = problem.Q*(sqrt(variances(branch))*standard_noise);
        deltas(inner) = weighted_norm(problem.W,effective_noise_star);
        observations(:,inner) = fitted_center+effective_noise_star;
        bootstrap.substream(inner,branch) = substream_star;
    end
    options_star = config.options;
    options_star.compute_debiased = branch == 2;
    % The batched solver selects each k* and replays that column's own k*.
    fit_star = pde_fit_batched(problem,observations,config.method,options_star,deltas);
    check_fit(fit_star,branch == 2);
    if branch == 1
        difference = fit_star.original-centers{branch};
        bootstrap.error_original = max(abs(problem.source_mass_root*difference),[],1).';
    else
        difference = fit_star.debiased-centers{branch};
        bootstrap.error_debiased = max(abs(problem.source_mass_root*difference),[],1).';
    end
    bootstrap.steps(:,branch) = fit_star.steps(:);
    bootstrap.cap_hit(:,branch) = fit_star.cap_hit(:);
    bootstrap.noise_norm(:,branch) = deltas(:);
end
original_error = fit.original-problem.pe;
debiased_error = fit.debiased-problem.pe;
q_original = quantile(bootstrap.error_original,config.confidence);
q_debiased = quantile(bootstrap.error_debiased,config.confidence);
% Preserve the production quadratic-mass E2 and complete mass-root Einf.
metrics = struct('E2_original',sqrt(max(0,original_error'*problem.M0*original_error)), ...
    'E2_debiased',sqrt(max(0,debiased_error'*problem.M0*debiased_error)), ...
    'Einf_original',max(abs(problem.source_mass_root*original_error)), ...
    'Einf_debiased',max(abs(problem.source_mass_root*debiased_error)), ...
    'sigma2_original',variances(1),'sigma2_debiased',variances(2), ...
    'variance_error_original',abs(variances(1)-config.sigma^2), ...
    'variance_error_debiased',abs(variances(2)-config.sigma^2), ...
    'quantile_original',q_original,'quantile_debiased',q_debiased, ...
    'coverage_original',0,'coverage_debiased',0,'steps',fit.steps, ...
    'cap_hit',fit.cap_hit,'noise_norm',noise_norm,'substream',substream);
metrics.coverage_original = metrics.Einf_original <= q_original;
metrics.coverage_debiased = metrics.Einf_debiased <= q_debiased;
names = fieldnames(metrics);
for j = 1:numel(names)
    value = metrics.(names{j});
    if ~isscalar(value) || ~isfinite(value)
        error('PDEFormal:NonfiniteMetric','Invalid scalar metric: %s',names{j});
    end
end
names = fieldnames(bootstrap);
for j = 1:numel(names)
    if any(~isfinite(bootstrap.(names{j})),'all')
        error('PDEFormal:NonfiniteBootstrap','Invalid bootstrap field: %s',names{j});
    end
end
row = struct('outer',outer,'metrics',metrics,'original',fit.original, ...
    'debiased',fit.debiased,'bootstrap',bootstrap,'seconds',toc(clock_start));
end

function value = weighted_norm(W,vector)
value = sqrt(max(0,vector'*W*vector));
if ~isfinite(value)
    error('PDEFormal:NonfiniteNoise','Effective noise norm is nonfinite.');
end
end

function check_fit(fit,needs_debiased)
values = [fit.original(:);fit.sigma2_original(:);fit.steps(:);double(fit.cap_hit(:))];
if needs_debiased
    values = [values;fit.debiased(:);fit.sigma2_debiased(:)];
end
if any(~isfinite(values)) || any(fit.sigma2_original(:) < 0) || ...
        (needs_debiased && any(fit.sigma2_debiased(:) < 0))
    error('PDEFormal:NonfiniteFit','A reconstruction or variance is invalid.');
end
end
