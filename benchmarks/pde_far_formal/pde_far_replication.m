function [row,cache] = pde_far_replication(problem,config,outer,cache)
%PDE_FAR_REPLICATION One FAR step-variant replication with an explicit backend.
% Uses every configured bootstrap sample in each of the two fitted-center
% branches. Random substreams, per-sample stops and fixed-k debias replay are
% inherited from pde_noise and the chosen validated fitting implementation.
% cache is an envelope with backend and payload; identity metadata, if present,
% is preserved. The worker checks run/manifest identity, not this row helper.
% The mathematical config does not contain cache data or backend switches.
clock_start = tic;
if nargin < 4 || isempty(cache)
    cache = struct('backend','cpu-boundary-batch','payload',[]);
end
if ~isstruct(cache) || ~isscalar(cache) || ~all(isfield(cache,{'backend','payload'})) || ...
        ~any(strcmp(cache.backend,{'far-spectral','cpu-boundary-batch'}))
    error('PDEFARFormal:Backend','Provide cache.backend and cache.payload for an explicitly selected backend.');
end
validateattributes(config,{'struct'},{'scalar'});
validateattributes(outer,{'numeric'},{'scalar','integer','finite','positive','<=',config.N});
validateattributes(config.B,{'numeric'},{'scalar','integer','finite','>=',2});
if ~strcmpi(config.method,'FAR')
    error('PDEFARFormal:Method','This replication helper supports FAR only.');
end
if ~config.options.compute_debiased
    error('PDEFARFormal:Debias','Both original and debiased estimates are required.');
end
cpu_fields = {'F','G','W','Q','pe','M0','source_mass_root','d0'};
if isfield(problem,'normal_boundary'), cpu_fields{end+1} = 'normal_boundary'; end
for j = 1:numel(cpu_fields)
    value = problem.(cpu_fields{j});
    if isa(value,'gpuArray') || ~isa(value,'double') || ~isreal(value)
        error('PDEFARFormal:Precision','problem.%s must be a real CPU double array.',cpu_fields{j});
    end
end
b = size(problem.F,1);
B = config.B;
[standard_noise,substream] = pde_noise(config.seed,outer,0,0,b);
effective_noise = problem.Q*(config.sigma*standard_noise);
noise_norm = weighted_norm(problem.W,effective_noise);
z = problem.F*problem.pe+effective_noise;
[fit,cache] = fit_samples(problem,z,config.options,noise_norm,cache);
check_fit(fit,true);
fit_diagnostics = struct('backend',cache.backend,'outer',fallback_info(fit), ...
    'bootstrap',{{}});
centers = {fit.original,fit.debiased};
variances = [fit.sigma2_original,fit.sigma2_debiased];
bootstrap = struct('error_original',nan(B,1),'error_debiased',nan(B,1), ...
    'steps',nan(B,2),'cap_hit',false(B,2),'noise_norm',nan(B,2), ...
    'substream',nan(B,2),'residual_norm',nan(B,2));
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
    % Both backends select each column's first k* and replay its own fixed k*.
    [fit_star,cache] = fit_samples(problem,observations,options_star,deltas,cache);
    check_fit(fit_star,branch == 2);
    fit_diagnostics.bootstrap{branch} = fallback_info(fit_star);
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
    bootstrap.residual_norm(:,branch) = fit_star.residual_norm_original(:);
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
        error('PDEFARFormal:NonfiniteMetric','Invalid scalar metric: %s',names{j});
    end
end
names = fieldnames(bootstrap);
for j = 1:numel(names)
    if any(~isfinite(bootstrap.(names{j})),'all')
        error('PDEFARFormal:NonfiniteBootstrap','Invalid bootstrap field: %s',names{j});
    end
end
row = struct('outer',outer,'metrics',metrics,'original',fit.original, ...
    'debiased',fit.debiased,'bootstrap',bootstrap,'seconds',toc(clock_start), ...
    'fit_diagnostics',fit_diagnostics);
end

function [fit,cache] = fit_samples(problem,Z,options,deltas,cache)
switch cache.backend
    case 'far-spectral'
        [fit,cache.payload] = pde_fit_far_spectral(problem,Z,options,deltas,cache.payload);
    case 'cpu-boundary-batch'
        fit = pde_fit_batched(problem,Z,'FAR',options,deltas);
    otherwise
        error('PDEFARFormal:Backend','Unsupported FAR execution backend.');
end
end

function value = fallback_info(fit)
value = [];
if isfield(fit,'spectral_fallback'), value = fit.spectral_fallback; end
end

function value = weighted_norm(W,vector)
value = sqrt(max(0,vector'*W*vector));
if ~isfinite(value)
    error('PDEFARFormal:NonfiniteNoise','Effective noise norm is nonfinite.');
end
end

function check_fit(fit,needs_debiased)
values = [fit.original(:);fit.sigma2_original(:);fit.steps(:);double(fit.cap_hit(:))];
if needs_debiased
    values = [values;fit.debiased(:);fit.sigma2_debiased(:)];
end
if any(~isfinite(values)) || any(fit.sigma2_original(:) < 0) || ...
        (needs_debiased && any(fit.sigma2_debiased(:) < 0))
    error('PDEFARFormal:NonfiniteFit','A reconstruction or variance is invalid.');
end
end
