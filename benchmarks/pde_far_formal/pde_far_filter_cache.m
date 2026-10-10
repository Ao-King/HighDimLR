function cache = pde_far_filter_cache(problem,options,required_steps,cache)
%PDE_FAR_FILTER_CACHE Full-mode FAR filters with an extendable, saveable table.
% cache = pde_far_filter_cache(problem,opt,K,cache) extends through step K.
% Omit K for min(4000,opt.max_steps); omit cache to construct a new cache.
% Explicit K is never clipped to opt.max_steps (fixed-step replay may exceed it).
% Tables contain k=0 in column 1. No eigenvalues or history terms are dropped.
% Use save(...,'cache','-v7.3') and load normally; no hidden persistent state.
if nargin < 2 || isempty(options), options = struct(); end
if ~isfield(options,'step'), options.step = 1; end
if ~isfield(options,'theta'), options.theta = 1.25; end
if ~isfield(options,'max_steps'), options.max_steps = 100000; end
validateattributes(options.step,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(options.theta,{'numeric'},{'scalar','real','finite','>',0,'<',2});
validateattributes(options.max_steps,{'numeric'},{'scalar','integer','finite','nonnegative'});
if nargin < 3 || isempty(required_steps), required_steps = min(4000,options.max_steps); end
validateattributes(required_steps,{'numeric'},{'scalar','integer','finite','nonnegative'});
if nargin < 4, cache = []; end
required_steps = double(required_steps);
F = problem.F; G = problem.G; W = problem.W;
if isfield(problem,'normal_boundary'), A = problem.normal_boundary; else, A = F*G; end
if isempty(cache)
    clock_start = tic;
    for value = {F,G,W,A}
        assert(isa(value{1},'double') && ~isa(value{1},'gpuArray') && ...
            isreal(value{1}) && all(isfinite(value{1}),'all'), ...
            'FARFilter:CPU','Operators must be finite real CPU double arrays.');
    end
    b = size(F,1);
    assert(isequal(size(A),[b,b]) && isequal(size(W),[b,b]) && ...
        isequal(size(G),[size(F,2),b]),'FARFilter:Dimensions','Inconsistent operator dimensions.');
    W_skew = norm(W-W.','fro')/max(norm(W,'fro'),realmin);
    symmetry_tolerance = 1000*eps*max(1,b);
    assert(W_skew <= symmetry_tolerance,'FARFilter:WeightSymmetry', ...
        'W is not symmetric to rounding accuracy.');
    [R,flag] = chol((W+W.')/2);
    assert(flag == 0 && rcond(R) > eps,'FARFilter:WeightSPD', ...
        'W must be numerically positive definite.');
    T = R*A/R;
    T_scale = max(norm(T,'fro'),realmin);
    T_skew = norm(T-T.','fro')/T_scale;
    assert(T_skew <= symmetry_tolerance,'FARFilter:OperatorSymmetry', ...
        'R*A/R is not symmetric to rounding accuracy; spectral backend refused.');
    symmetric_T = (T+T.')/2;
    [U,lambda] = eig(symmetric_T,'vector');
    [lambda,order] = sort(lambda);
    U = U(:,order);
    psd_tolerance = 1000*eps*max(1,b)*max(max(abs(lambda)),realmin);
    assert(all(isfinite(lambda)) && min(lambda) >= -psd_tolerance, ...
        'FARFilter:OperatorPSD','The transformed operator has a significant negative eigenvalue.');
    % Retain signed roundoff-level negative eigenvalues and ALL tiny modes.
    % No lambda division is used anywhere in the recurrence or reconstruction.
    V = R\U;
    L = U.'*R;
    reconstructed_A = (V.*lambda.')*L;
    factorization = struct('W_skew_relative',W_skew,'T_skew_relative',T_skew, ...
        'symmetry_tolerance',symmetry_tolerance,'condition_R',cond(R), ...
        'operator_relative_error',norm(reconstructed_A-A,'fro')/max(norm(A,'fro'),realmin), ...
        'eigen_relative_residual',norm(symmetric_T*U-U.*lambda.','fro')/T_scale, ...
        'basis_inverse_residual',norm(V*L-eye(b),'fro')/sqrt(b), ...
        'minimum_eigenvalue',min(lambda),'maximum_eigenvalue',max(lambda), ...
        'psd_tolerance',psd_tolerance,'retained_modes',b,'truncated_modes',0);
    assert(factorization.operator_relative_error <= 10*symmetry_tolerance && ...
        factorization.basis_inverse_residual <= 10*symmetry_tolerance, ...
        'FARFilter:Factorization','Spectral factorization is insufficiently accurate.');
    cache = struct('version','pde-far-spectral-cache-v1','step',double(options.step), ...
        'theta',double(options.theta),'computed_steps',0,'s',zeros(b,1), ...
        'residual_factors',ones(b,1),'lambda',lambda,'L',L,'V',V, ...
        'source_basis',G*V,'R',R,'U',U,'operator_A',A, ...
        'reference_F',F,'reference_G',G,'reference_W',W,'factorization',factorization, ...
        'creation_seconds',toc(clock_start),'extension_seconds',0);
else
    fields = {'version','step','theta','computed_steps','s','residual_factors', ...
        'lambda','L','V','source_basis','R','U','operator_A','reference_F','reference_G', ...
        'reference_W','factorization','creation_seconds','extension_seconds'};
    assert(isstruct(cache) && isscalar(cache) && all(isfield(cache,fields)) && ...
        strcmp(cache.version,'pde-far-spectral-cache-v1'), ...
        'FARFilter:CacheVersion','Unsupported or incomplete FAR cache.');
    assert(isequal(cache.step,double(options.step)) && isequal(cache.theta,double(options.theta)) && ...
        isequal(cache.reference_F,F) && isequal(cache.reference_G,G) && ...
        isequal(cache.reference_W,W) && isequal(cache.operator_A,A), ...
        'FARFilter:CacheIdentity','FAR cache does not match the operators, step or theta.');
    validateattributes(cache.computed_steps,{'numeric'},{'scalar','integer','finite','nonnegative'});
    b = size(F,1);
    assert(isequal(size(cache.s),[b,cache.computed_steps+1]) && ...
        isequal(size(cache.residual_factors),size(cache.s)) && ...
        isequal(size(cache.lambda),[b,1]) && cache.factorization.retained_modes == b && ...
        cache.factorization.truncated_modes == 0,'FARFilter:CacheShape','Invalid FAR table dimensions.');
end
if required_steps <= cache.computed_steps, return; end
extension_clock = tic;
first = cache.computed_steps+1;
cache.s(:,required_steps+1) = 0;
cache.residual_factors(:,required_steps+1) = 0;
h = cache.step;
theta = cache.theta;
gamma_theta = gamma(theta);
scale = h^theta/(theta*(theta+1));
for i = first:required_steps
    % Exactly the production FAR predictor/corrector weights, applied to all
    % eigenmodes simultaneously. Column j stores r_(j-1), not r_j.
    history = cache.residual_factors(:,1:i);
    lag = i-1:-1:0;
    weights_b = h^theta/theta*((lag+1).^theta-lag.^theta);
    predicted_s = (history*weights_b.')/gamma_theta;
    weights_a = zeros(1,i);
    weights_a(1) = (i-1)^(theta+1)-(i-1-theta)*i^theta;
    if i > 1
        lag = i-2:-1:0;
        weights_a(2:end) = (lag+2).^(theta+1)+lag.^(theta+1)-2*(lag+1).^(theta+1);
    end
    next_s = scale*((1-cache.lambda.*predicted_s)+history*weights_a.')/gamma_theta;
    next_r = 1-cache.lambda.*next_s;
    assert(all(isfinite(next_s)) && all(isfinite(next_r)), ...
        'FARFilter:Nonfinite','The FAR filter became nonfinite at step %d.',i);
    cache.s(:,i+1) = next_s;
    cache.residual_factors(:,i+1) = next_r;
end
cache.computed_steps = required_steps;
cache.extension_seconds = cache.extension_seconds+toc(extension_clock);
end
