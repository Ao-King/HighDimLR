function [result,cache] = pde_fit_far_spectral(problem,Z,options,noise_norm,cache,fixed_steps)
%PDE_FIT_FAR_SPECTRAL Full-mode CPU FAR with exact first-hit scan and guards.
% [fit,cache] = pde_fit_far_spectral(problem,Z,opt,delta,cache,fixed_steps)
% Columns are independent samples. The cache can be saved/loaded and grows
% when needed; max_steps remains the actual supplied ceiling, never 4000.
% Every integer k is scanned, including k=0, without a monotonicity assumption.
% A wide first pass records ALL near-boundary candidates before the first hit.
% A second pass reconstructs their source-space residuals and estimates the
% floating-point guard margin. Unresolved cases or inconsistent decisions
% fall back to the unchanged source-space pde_fit reference implementation.
if nargin < 3 || isempty(options), options = struct(); end
if nargin < 4 || isempty(noise_norm), noise_norm = 0; end
if nargin < 5, cache = []; end
if nargin < 6, fixed_steps = []; end
options = complete_options(options);
assert(isa(Z,'double') && ~isa(Z,'gpuArray') && isreal(Z) && ismatrix(Z) && ...
    all(isfinite(Z),'all'),'FARSpectral:Observation','Z must be a finite real CPU double matrix.');
[b,m] = size(Z);
assert(b == size(problem.F,1) && m >= 1,'FARSpectral:ObservationSize','Invalid observation dimensions.');
noise_norm = sample_row(noise_norm,m,false);
is_fixed = ~isempty(fixed_steps);
if is_fixed
    fixed_steps = sample_row(fixed_steps,m,true);
    limit = max(fixed_steps);
else
    limit = options.max_steps;
end
if isempty(cache)
    cache = pde_far_filter_cache(problem,options,min(4000,limit));
else
    % A zero request validates a loaded cache without rebuilding its table.
    cache = pde_far_filter_cache(problem,options,0,cache);
end
coordinates = cache.L*Z;
threshold = options.tau*noise_norm;
threshold_squared = threshold.^2;
near_steps = zeros(1,0);
near_columns = zeros(1,0);
near_hit = false(1,0);
cap_hit = false(1,m);
ambiguity_relative_tolerance = max(1e-10,1000*cache.factorization.operator_relative_error);
if is_fixed
    cache = pde_far_filter_cache(problem,options,limit,cache);
    steps = fixed_steps;
else
    steps = nan(1,m);
    coordinate_energy = coordinates.^2;
    initial_energy = sum(coordinate_energy,1);
    ambiguity_band = ambiguity_relative_tolerance*max([initial_energy;threshold_squared;repmat(realmin,1,m)],[],1);
    next_scan = 0;
    while any(isnan(steps))
        available = min(cache.computed_steps,limit);
        while next_scan <= available && any(isnan(steps))
            last_scan = min(available,next_scan+511);
            indices = find(isnan(steps));
            factors = cache.residual_factors(:,next_scan+1:last_scan+1);
            energies = (factors.^2).'*coordinate_energy(:,indices);
            assert(all(isfinite(energies),'all'),'FARSpectral:ResidualEnergy','Nonfinite spectral residual energy.');
            hits = energies <= threshold_squared(indices);
            [found,first_hit] = max(hits,[],1);
            found = logical(found);
            prefix_last = first_hit;
            prefix_last(~found) = size(energies,1);
            prefix = (1:size(energies,1)).' <= prefix_last;
            near = abs(energies-threshold_squared(indices)) <= ambiguity_band(indices);
            [near_row,near_column] = find(near & prefix);
            if ~isempty(near_row)
                near_linear = sub2ind(size(hits),near_row,near_column);
                near_steps = [near_steps,reshape(next_scan+near_row-1,1,[])]; %#ok<AGROW>
                near_columns = [near_columns,reshape(indices(near_column),1,[])]; %#ok<AGROW>
                near_hit = [near_hit,reshape(hits(near_linear),1,[])]; %#ok<AGROW>
            end
            steps(indices(found)) = next_scan+first_hit(found)-1;
            next_scan = last_scan+1;
        end
        unresolved = isnan(steps);
        if ~any(unresolved), break; end
        if available == limit
            steps(unresolved) = limit;
            cap_hit(unresolved) = true;
            break
        end
        next_limit = min(limit,max(cache.computed_steps+1,max(1,2*cache.computed_steps)));
        cache = pde_far_filter_cache(problem,options,next_limit,cache);
    end
end
[fallback,guard] = refine_candidates(problem,Z,coordinates,threshold,cache, ...
    near_steps,near_columns,near_hit);
selected_filters = cache.s(:,steps+1);
original = reconstruct(problem,cache,selected_filters.*coordinates);
residual_original = Z-problem.F*original;
norm_original = boundary_norm(problem,residual_original);
reconstruction_guard = false(1,m);
if ~is_fixed
    % Check source-space residuals at the selected k and its predecessor.
    % These checks also detect a capped result that actually meets the rule.
    reconstruction_guard = (norm_original <= threshold) ~= ~cap_hit;
    positive = steps > 0;
    if any(positive)
        previous_filters = cache.s(:,steps(positive));
        previous = reconstruct(problem,cache,previous_filters.*coordinates(:,positive));
        previous_norm = boundary_norm(problem,Z(:,positive)-problem.F*previous);
        reconstruction_guard(positive) = reconstruction_guard(positive) | previous_norm <= threshold(positive);
    end
end
fallback = fallback | reconstruction_guard;
debiased = [];
residual_debiased = [];
norm_debiased = nan(1,m);
if options.compute_debiased
    % In exact arithmetic this is diag(2*s-lambda*s.^2)*coordinates.
    % Replay the actual F*original response to match production association.
    fitted_response = problem.F*original;
    refit = reconstruct(problem,cache,selected_filters.*(cache.L*fitted_response));
    debiased = 2*original-refit;
end
for column = reshape(find(fallback),1,[])
    if is_fixed, replay = fixed_steps(column); else, replay = []; end
    reference = pde_fit(problem,Z(:,column),'FAR',options,noise_norm(column),replay);
    original(:,column) = reference.original;
    steps(column) = reference.steps;
    cap_hit(column) = reference.cap_hit;
    if options.compute_debiased, debiased(:,column) = reference.debiased; end
end
assert(all(isfinite(original),'all'),'FARSpectral:Source','Nonfinite reconstructed source.');
residual_original = Z-problem.F*original;
norm_original = boundary_norm(problem,residual_original);
if options.compute_debiased
    assert(all(isfinite(debiased),'all'),'FARSpectral:Source','Nonfinite debiased source.');
    residual_debiased = Z-problem.F*debiased;
    norm_debiased = boundary_norm(problem,residual_debiased);
end
result = struct('original',original,'debiased',debiased,'steps',steps, ...
    'residual_original',residual_original,'residual_debiased',residual_debiased, ...
    'residual_norm',norm_original,'residual_norm_original',norm_original, ...
    'residual_norm_debiased',norm_debiased,'sigma2_original',norm_original.^2/problem.d0, ...
    'sigma2_debiased',norm_debiased.^2/problem.d0,'cap_hit',cap_hit, ...
    'method','far','options',options,'noise_norm',noise_norm,'spectral_fallback',fallback, ...
    'spectral',struct('cache_version',cache.version,'cache_steps',cache.computed_steps, ...
        'fallback_count',nnz(fallback),'ambiguity_relative_tolerance',ambiguity_relative_tolerance, ...
        'guard_version','source-residual-two-stage-v2', ...
        'source_reconstruction','cached-GV-full-mode', ...
        'near_candidate_count',numel(near_steps), ...
        'near_candidate_column_count',numel(unique(near_columns)), ...
        'refined_guard_count',guard.refined_guard_count, ...
        'refined_guard_column_count',guard.refined_guard_column_count, ...
        'refined_sign_mismatch_count',guard.sign_mismatch_count, ...
        'reconstruction_guard_column_count',nnz(reconstruction_guard), ...
        'maximum_refined_bound',guard.maximum_bound, ...
        'minimum_refined_margin',guard.minimum_margin, ...
        'roundoff_safety_factor',10,'factorization_safety_factor',1000, ...
        'retained_modes',b,'factorization',cache.factorization));
end

function [fallback,diagnostics] = refine_candidates(problem,Z,coordinates,threshold,cache,steps,columns,spectral_hit)
% The componentwise estimates below cover the explicit L, P=G*V, F, R
% products and subtraction in source reconstruction, including the rounding
% when the cache formed P. Combining G and V before multiplying mode
% coefficients preserves cancellation within each source-basis column.
% The extra factorization
% term is a conservative numerical guard, not a rigorous global forward-error
% theorem for the full fractional-history recurrence. Production comparison
% remains a prerequisite for using this backend.
[b,m] = size(Z);
n = size(problem.F,2);
fallback = false(1,m);
diagnostics = struct('refined_guard_count',0,'refined_guard_column_count',0, ...
    'sign_mismatch_count',0,'maximum_bound',0,'minimum_margin',0);
if isempty(steps), return; end
gamma_b = rounding_gamma(b);
gamma_n = rounding_gamma(n);
gamma_energy = rounding_gamma(2*b+1);
abs_L = abs(cache.L); abs_V = abs(cache.V); abs_R = abs(cache.R);
abs_G = abs(problem.G); abs_F = abs(problem.F); abs_W = abs(problem.W);
abs_P = abs(cache.source_basis);
source_basis_error = gamma_b*(abs_G*abs_V);
coordinate_error = gamma_b*(abs_L*abs(Z));
initial_norm = sqrt(sum((cache.R*Z).^2,1));
factorization_sum = cache.factorization.operator_relative_error+ ...
    cache.factorization.basis_inverse_residual+cache.factorization.eigen_relative_residual;
factorization_bound = 1000*factorization_sum*initial_norm;
diagnostics.minimum_margin = Inf;
for first = 1:256:numel(steps)
    selected = first:min(first+255,numel(steps));
    sample = columns(selected);
    filters = cache.s(:,steps(selected)+1);
    coefficients = filters.*coordinates(:,sample);
    source = cache.source_basis*coefficients;
    fitted_response = problem.F*source;
    residual = Z(:,sample)-fitted_response;
    residual_norm = boundary_norm(problem,residual);
    actual_hit = residual_norm <= threshold(sample);
    coefficient_error = abs(filters).*coordinate_error(:,sample)+eps*abs(coefficients);
    source_error = abs_P*coefficient_error+gamma_b*(abs_P*abs(coefficients))+ ...
        source_basis_error*abs(coefficients);
    residual_error = abs_F*source_error+gamma_n*(abs_F*abs(source))+ ...
        eps*(abs(Z(:,sample))+abs(fitted_response));
    norm_error = sqrt(sum((abs_R*residual_error).^2,1))+ ...
        gamma_b*sqrt(sum((abs_R*abs(residual)).^2,1));
    % Also cover the W quadratic form actually used by boundary_norm rather
    % than assuming it was evaluated by norm(R*residual).
    energy_error = gamma_energy*sum(abs(residual).*(abs_W*abs(residual)),1);
    energy_norm_error = energy_error./max(sqrt(residual_norm.^2+energy_error)+residual_norm,realmin);
    bound = 10*(norm_error+energy_norm_error+eps*abs(threshold(sample)))+factorization_bound(sample);
    assert(all(isfinite(bound)) && all(bound >= 0),'FARSpectral:Guard','Invalid refined roundoff guard.');
    margin = abs(residual_norm-threshold(sample));
    sign_mismatch = actual_hit ~= spectral_hit(selected);
    needs_reference = sign_mismatch | margin <= bound;
    fallback(sample(needs_reference)) = true;
    diagnostics.refined_guard_count = diagnostics.refined_guard_count+nnz(needs_reference);
    diagnostics.sign_mismatch_count = diagnostics.sign_mismatch_count+nnz(sign_mismatch);
    diagnostics.maximum_bound = max(diagnostics.maximum_bound,max(bound));
    diagnostics.minimum_margin = min(diagnostics.minimum_margin,min(margin));
end
diagnostics.refined_guard_column_count = nnz(fallback);
end

function value = rounding_gamma(count)
assert(count*eps < 1,'FARSpectral:Guard','Dimensions exceed the floating-point guard range.');
value = count*eps/(1-count*eps);
end

function source = reconstruct(problem,cache,coefficients)
% Every eigenmode is retained. This association avoids forming large
% boundary-coordinate components that subsequently cancel through G.
source = cache.source_basis*coefficients;
end

function value = boundary_norm(problem,residual)
energy = sum(residual.*(problem.W*residual),1);
scale = sum(residual.^2,1)*norm(problem.W,'fro');
assert(all(isfinite(energy)) && all(energy >= -100*eps*max(1,scale)), ...
    'FARSpectral:ResidualEnergy','Invalid source-space residual energy.');
value = sqrt(max(0,energy));
end

function row = sample_row(value,m,integer_required)
assert(~isa(value,'gpuArray'),'FARSpectral:CPU','Sample metadata must be on the CPU.');
validateattributes(value,{'numeric'},{'vector','real','finite','nonnegative'});
row = double(reshape(value,1,[]));
if isscalar(row), row = repmat(row,1,m); end
assert(numel(row) == m && (~integer_required || all(row == fix(row))), ...
    'FARSpectral:SampleParameter','Expected a scalar or one valid value per observation column.');
end

function options = complete_options(options)
validateattributes(options,{'struct'},{'scalar'});
if isfield(options,'alpha') && ~isfield(options,'omega'), options.omega = options.alpha; end
defaults = struct('step',1,'tau',1.1,'max_steps',100000,'eta',1, ...
    's',1.5,'initial_time',1,'omega',3,'kappa',1,'theta',1.25,'compute_debiased',true);
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(options,names{j}), options.(names{j}) = defaults.(names{j}); end
end
for name = {'step','tau','initial_time'}
    validateattributes(options.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(options.max_steps,{'numeric'},{'scalar','integer','finite','nonnegative'});
validateattributes(options.eta,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(options.s,{'numeric'},{'scalar','real','finite','>',-0.5});
validateattributes(options.omega,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(options.kappa,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(options.theta,{'numeric'},{'scalar','real','finite','>',0,'<',2});
validateattributes(options.compute_debiased,{'numeric','logical'},{'scalar','real','finite'});
options.compute_debiased = logical(options.compute_debiased);
end
