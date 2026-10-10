function result = pde_fit(problem, z, method, options, noise_norm, fixed_steps)
%PDE_FIT Fit a PDE source and replay the same fixed-step map for debiasing.
if nargin < 4 || isempty(options)
    options = struct;
end
if nargin < 5
    noise_norm = 0;
end
if nargin < 6
    fixed_steps = [];
end
method = lower(char(method));
if ~ismember(method, {'landweber', 'showalter', 'hbf', 'soar', 'nesterov', 'far', 'ark'})
    error('pde_fit:Method', 'Unknown PDE method.');
end
options = complete_options(options);
validateattributes(noise_norm, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(z, {'numeric'}, {'vector', 'real', 'finite'});
z = z(:);
if numel(z) ~= size(problem.F, 1)
    error('pde_fit:ObservationSize', 'z must have one value per boundary node.');
end
if ~isempty(fixed_steps)
    validateattributes(fixed_steps, {'numeric'}, {'scalar', 'integer', 'finite', 'nonnegative'});
end
[original, steps, cap_hit] = run_kernel(problem, z, method, options, noise_norm, fixed_steps);
residual_original = z - problem.F * original;
norm_original = boundary_norm(problem, residual_original);
debiased = [];
residual_debiased = [];
norm_debiased = NaN;
if options.compute_debiased
    fitted_response = problem.F * original;
    refit = run_kernel(problem, fitted_response, method, options, 0, steps);
    debiased = 2 * original - refit;
    residual_debiased = z - problem.F * debiased;
    norm_debiased = boundary_norm(problem, residual_debiased);
end
result = struct('original', original, 'debiased', debiased, 'steps', steps, ...
    'residual_original', residual_original, 'residual_debiased', residual_debiased, ...
    'residual_norm', norm_original, 'residual_norm_original', norm_original, ...
    'residual_norm_debiased', norm_debiased, 'sigma2_original', norm_original^2 / problem.d0, ...
    'sigma2_debiased', norm_debiased^2 / problem.d0, 'cap_hit', cap_hit, ...
    'method', method, 'options', options, 'noise_norm', noise_norm);
end

function options = complete_options(options)
if ~isstruct(options) || ~isscalar(options)
    error('pde_fit:Options', 'options must be a scalar structure.');
end
if isfield(options, 'alpha') && ~isfield(options, 'omega')
    options.omega = options.alpha;
end
defaults = struct('step', 1, 'tau', 1.1, 'max_steps', 100000, 'eta', 1, ...
    's', 1.5, 'initial_time', 1, 'omega', 3, 'kappa', 1, 'theta', 1.25, ...
    'compute_debiased', true);
fields = fieldnames(defaults);
for j = 1:numel(fields)
    key = fields{j};
    if ~isfield(options, key)
        options.(key) = defaults.(key);
    end
end
positive = {'step', 'tau', 'initial_time'};
for j = 1:numel(positive)
    validateattributes(options.(positive{j}), {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
end
validateattributes(options.max_steps, {'numeric'}, {'scalar', 'integer', 'finite', 'nonnegative'});
validateattributes(options.eta, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(options.s, {'numeric'}, {'scalar', 'real', 'finite', '>', -0.5});
validateattributes(options.omega, {'numeric'}, {'scalar', 'real', 'finite', '>', -1});
validateattributes(options.kappa, {'numeric'}, {'scalar', 'real', 'finite', '>', -1});
validateattributes(options.theta, {'numeric'}, {'scalar', 'real', 'finite', '>', 0, '<', 2});
validateattributes(options.compute_debiased, {'logical', 'numeric'}, {'scalar', 'real', 'finite'});
options.compute_debiased = logical(options.compute_debiased);
end

function [p, k, cap_hit] = run_kernel(problem, z, method, opt, noise_norm, fixed_steps)
n = size(problem.F, 2);
p = zeros(n, 1);
previous = p;
velocity = p;
history = zeros(size(problem.F, 1), 0);
k = 0;
h = opt.step;
threshold = opt.tau * noise_norm;
cap_hit = false;
if isempty(fixed_steps)
    limit = opt.max_steps;
else
    limit = fixed_steps;
end
if strcmp(method, 'ark')
    if isfield(problem, 'normal_boundary')
        normal_boundary = problem.normal_boundary;
    else
        normal_boundary = problem.F * problem.G;
    end
end
for i = 1:limit
    if isempty(fixed_steps) && boundary_norm(problem, z - problem.F * p) <= threshold
        break
    end
    switch method
        case 'landweber'
            p = p + h * force(problem, z, p);
        case 'showalter'
            k1 = force(problem, z, p);
            k2 = force(problem, z, p + h * k1 / 2);
            k3 = force(problem, z, p + h * k2 / 2);
            k4 = force(problem, z, p + h * k3);
            p = p + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
        case 'hbf'
            p1 = velocity;
            v1 = force(problem, z, p) - opt.eta * velocity;
            p2 = velocity + h * v1 / 2;
            v2 = force(problem, z, p + h * p1 / 2) - opt.eta * p2;
            p3 = velocity + h * v2 / 2;
            v3 = force(problem, z, p + h * p2 / 2) - opt.eta * p3;
            p4 = velocity + h * v3;
            v4 = force(problem, z, p + h * p3) - opt.eta * p4;
            p = p + h * (p1 + 2*p2 + 2*p3 + p4) / 6;
            velocity = velocity + h * (v1 + 2*v2 + 2*v3 + v4) / 6;
        case 'soar'
            t = opt.initial_time + (i-1) * h;
            damping = (1 + 2*opt.s) / t;
            half_velocity = (velocity + h * force(problem, z, p) / 2) / (1 + h*damping/2);
            p = p + h * half_velocity;
            next_damping = (1 + 2*opt.s) / (t + h);
            velocity = half_velocity + h * (force(problem, z, p) - next_damping*half_velocity) / 2;
        case 'nesterov'
            if i == 1
                extrapolated = p;
            else
                extrapolated = p + (i-2)/(i-1+opt.omega) * (p-previous);
            end
            previous = p;
            p = extrapolated + h * force(problem, z, extrapolated);
        case 'far'
            theta = opt.theta;
            if i > size(history, 2)
                history(:, min(limit, max(i, 2*max(128, size(history, 2))))) = 0;
            end
            history(:, i) = z - problem.F*p;
            lag = i-1:-1:0;
            b = h^theta/theta * ((lag+1).^theta - lag.^theta);
            predicted = problem.G * (history(:, 1:i) * b.') / gamma(theta);
            a = zeros(1, i);
            a(1) = (i-1)^(theta+1) - (i-1-theta)*i^theta;
            if i > 1
                lag = i-2:-1:0;
                a(2:end) = (lag+2).^(theta+1) + lag.^(theta+1) - 2*(lag+1).^(theta+1);
            end
            scale = h^theta / (theta*(theta+1));
            p = scale * (force(problem, z, predicted) + problem.G*(history(:, 1:i)*a.')) / gamma(theta);
        case 'ark'
            p = p + h * velocity;
            t = h*i;
            alpha = 1 + (t^(-opt.kappa)-opt.kappa)/i;
            coefficient = h*t^opt.kappa;
            rhs = velocity + h * force(problem, z, p);
            if ~isfinite(alpha) || alpha <= 0 || ~isfinite(coefficient)
                error('pde_fit:ARKCoefficient', 'AR-kappa has a singular or nonfinite step coefficient.');
            end
            % Woodbury solves the implicit source update in boundary space.
            H = alpha*eye(size(problem.F, 1)) + coefficient*normal_boundary;
            if rcond(H) < eps
                error('pde_fit:ARKSolve', 'The AR-kappa step matrix is numerically singular.');
            end
            velocity = (rhs - coefficient*problem.G*(H \ (problem.F*rhs))) / alpha;
    end
    if any(~isfinite(p)) || any(~isfinite(velocity))
        error('pde_fit:NonFiniteIterate', 'The iteration produced a nonfinite state at step %d.', i);
    end
    k = i;
end
if isempty(fixed_steps)
    cap_hit = k == limit && boundary_norm(problem, z - problem.F*p) > threshold;
end
end

function value = force(problem, z, p)
value = problem.G * (z - problem.F * p);
end

function value = boundary_norm(problem, r)
energy = r.' * problem.W * r;
if ~isfinite(energy) || energy < -100*eps*max(1, norm(r)^2 * norm(problem.W, 'fro'))
    error('pde_fit:ResidualEnergy', 'The boundary residual energy is invalid.');
end
value = sqrt(max(0, energy));
end
