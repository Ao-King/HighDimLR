function result = pde_fit_batched(problem, Z, method, options, noise_norm, fixed_steps)
%PDE_FIT_BATCHED Double-precision boundary-coordinate solver.
% Each column is an independent pde_fit call. Z selects CPU or GPU execution.
% The source state p=G*q is represented by q, so its force is Z-(F*G)*q.
% This changes floating-point association, not the mathematical iteration.
% Stopped columns remain frozen; debiasing replays each column's own step count.
% Source/residual outputs remain on Z's device. Step and cap metadata are CPU.
if nargin < 4 || isempty(options), options = struct; end
if nargin < 5 || isempty(noise_norm), noise_norm = 0; end
if nargin < 6, fixed_steps = []; end
method = lower(char(method));
if ~ismember(method, {'landweber','showalter','hbf','soar','nesterov','far','ark'})
    error('pde_fit_batched:Method', 'Unknown PDE method.');
end
options = complete_options(options);
if ~(isnumeric(Z) || isa(Z,'gpuArray')) || ~isreal(Z) || ~ismatrix(Z)
    error('pde_fit_batched:Observation', 'Z must be a real numeric matrix.');
end
Z = double(Z);
if host(any(~isfinite(Z),'all'))
    error('pde_fit_batched:Observation', 'Z must be finite.');
end
[b,m] = size(Z);
if b ~= size(problem.F,1) || m < 1
    error('pde_fit_batched:ObservationSize', 'Z must have one row per boundary node and at least one column.');
end
noise_norm = sample_row(noise_norm,m,'noise_norm',false);
if ~isempty(fixed_steps)
    fixed_steps = sample_row(fixed_steps,m,'fixed_steps',true);
end
P = struct('F',on_device(problem.F,Z),'G',on_device(problem.G,Z), ...
    'W',on_device(problem.W,Z),'d0',double(host(problem.d0)));
if isfield(problem,'normal_boundary')
    P.A = on_device(problem.normal_boundary,Z);
else
    P.A = P.F*P.G;
end
P.W_fro = host(norm(P.W,'fro'));
if strcmp(method,'ark')
    % rcond is checked on the host to retain the production singularity guard
    % without requiring GPU rcond support. Its cost is included in timing.
    P.A_host = host(P.A);
end
[q,steps,cap_hit] = run_kernel(P,Z,method,options,noise_norm,fixed_steps);
original = P.G*q;
check_source(original);
residual_original = Z-P.F*original;
norm_original = boundary_norm(P,residual_original);
debiased = [];
residual_debiased = [];
norm_debiased = nan(1,m,'like',Z);
if options.compute_debiased
    fitted_response = P.F*original;
    q_refit = run_kernel(P,fitted_response,method,options,zeros(1,m),steps);
    debiased = 2*original-P.G*q_refit;
    check_source(debiased);
    residual_debiased = Z-P.F*debiased;
    norm_debiased = boundary_norm(P,residual_debiased);
end
result = struct('original',original,'debiased',debiased,'steps',steps, ...
    'residual_original',residual_original,'residual_debiased',residual_debiased, ...
    'residual_norm',norm_original,'residual_norm_original',norm_original, ...
    'residual_norm_debiased',norm_debiased, ...
    'sigma2_original',norm_original.^2/P.d0, ...
    'sigma2_debiased',norm_debiased.^2/P.d0,'cap_hit',cap_hit, ...
    'method',method,'options',options,'noise_norm',noise_norm);
end

function [q,steps,cap_hit] = run_kernel(P,Z,method,opt,noise_norm,fixed_steps)
[b,m] = size(Z);
q = zeros(b,m,'like',Z);
previous = q;
velocity = q;
history = zeros(b,m,0,'like',Z);
steps = zeros(1,m);
cap_hit = false(1,m);
h = opt.step;
threshold = opt.tau*noise_norm;
is_fixed = ~isempty(fixed_steps);
if is_fixed
    limits = fixed_steps;
else
    limits = repmat(opt.max_steps,1,m);
end
limit = max(limits);
active = true(1,m);
for completed = 0:limit
    % One gathered packet checks all columns and their residual energies.
    % Checking immediately after an update also covers the final capped step.
    R = Z-P.A*q;
    energy = sum(R.*(P.W*R),1);
    packet = host([energy;sum(R.^2,1); ...
        any(~isfinite(q),1);any(~isfinite(velocity),1)]);
    if any(packet(3:4,:),'all')
        error('pde_fit_batched:NonFiniteIterate', ...
            'The iteration produced a nonfinite state at step %d.',completed);
    end
    validate_energy(packet(1,:),packet(2,:),P.W_fro);
    residual_norm_host = sqrt(max(0,packet(1,:)));
    active = active & steps < limits;
    if ~is_fixed
        active = active & residual_norm_host > threshold;
    end
    if ~any(active), break; end
    i = completed+1;
    qa = q(:,active);
    va = velocity(:,active);
    za = Z(:,active);
    force = @(state) za-P.A*state;
    switch method
        case 'landweber'
            next_q = qa+h*force(qa);
        case 'showalter'
            k1 = force(qa);
            k2 = force(qa+h*k1/2);
            k3 = force(qa+h*k2/2);
            k4 = force(qa+h*k3);
            next_q = qa+h*(k1+2*k2+2*k3+k4)/6;
        case 'hbf'
            p1 = va;
            v1 = force(qa)-opt.eta*va;
            p2 = va+h*v1/2;
            v2 = force(qa+h*p1/2)-opt.eta*p2;
            p3 = va+h*v2/2;
            v3 = force(qa+h*p2/2)-opt.eta*p3;
            p4 = va+h*v3;
            v4 = force(qa+h*p3)-opt.eta*p4;
            next_q = qa+h*(p1+2*p2+2*p3+p4)/6;
            velocity(:,active) = va+h*(v1+2*v2+2*v3+v4)/6;
        case 'soar'
            t = opt.initial_time+(i-1)*h;
            damping = (1+2*opt.s)/t;
            half_velocity = (va+h*force(qa)/2)/(1+h*damping/2);
            next_q = qa+h*half_velocity;
            next_damping = (1+2*opt.s)/(t+h);
            velocity(:,active) = half_velocity+h*(force(next_q)-next_damping*half_velocity)/2;
        case 'nesterov'
            if i == 1
                extrapolated = qa;
            else
                extrapolated = qa+(i-2)/(i-1+opt.omega)*(qa-previous(:,active));
            end
            previous(:,active) = qa;
            next_q = extrapolated+h*force(extrapolated);
        case 'far'
            theta = opt.theta;
            capacity = size(history,3);
            if i > capacity
                new_capacity = min(limit,max(i,2*max(128,capacity)));
                history(b,m,new_capacity) = 0;
            end
            history(:,active,i) = R(:,active);
            lag = i-1:-1:0;
            weights_b = h^theta/theta*((lag+1).^theta-lag.^theta);
            flat_history = reshape(history(:,:,1:i),b*m,i);
            predicted = reshape(flat_history*on_device(weights_b.',Z),b,m)/gamma(theta);
            weights_a = zeros(1,i);
            weights_a(1) = (i-1)^(theta+1)-(i-1-theta)*i^theta;
            if i > 1
                lag = i-2:-1:0;
                weights_a(2:end) = (lag+2).^(theta+1)+lag.^(theta+1)-2*(lag+1).^(theta+1);
            end
            accumulated = reshape(flat_history*on_device(weights_a.',Z),b,m);
            scale = h^theta/(theta*(theta+1));
            next_q = scale*(force(predicted(:,active))+accumulated(:,active))/gamma(theta);
        case 'ark'
            next_q = qa+h*va;
            t = h*i;
            alpha = 1+(t^(-opt.kappa)-opt.kappa)/i;
            coefficient = h*t^opt.kappa;
            rhs = va+h*force(next_q);
            if ~isfinite(alpha) || alpha <= 0 || ~isfinite(coefficient)
                error('pde_fit_batched:ARKCoefficient', ...
                    'AR-kappa has a singular or nonfinite step coefficient.');
            end
            H_host = alpha*eye(b)+coefficient*P.A_host;
            if rcond(H_host) < eps
                error('pde_fit_batched:ARKSolve', ...
                    'The AR-kappa step matrix is numerically singular.');
            end
            H = on_device(H_host,Z);
            % Keep the same Woodbury expression as pde_fit, represented in q.
            velocity(:,active) = (rhs-coefficient*(H\(P.A*rhs)))/alpha;
    end
    q(:,active) = next_q;
    steps(active) = i;
end
if ~is_fixed
    cap_hit = steps == limits & residual_norm_host > threshold;
end
end

function value = boundary_norm(P,R)
energy = sum(R.*(P.W*R),1);
packet = host([energy;sum(R.^2,1)]);
validate_energy(packet(1,:),packet(2,:),P.W_fro);
value = sqrt(max(0,energy));
end

function validate_energy(energy,residual_squared_norm,W_fro)
if any(~isfinite(energy)) || any(energy < -100*eps*max(1,residual_squared_norm*W_fro))
    error('pde_fit_batched:ResidualEnergy', 'The boundary residual energy is invalid.');
end
end

function check_source(source)
if host(any(~isfinite(source),'all'))
    error('pde_fit_batched:NonFiniteSource', 'Source reconstruction produced a nonfinite value.');
end
end

function row = sample_row(value,m,name,integer_required)
row = double(host(value));
validateattributes(row,{'numeric'},{'vector','real','finite','nonnegative'},mfilename,name);
if isscalar(row), row = repmat(row,1,m); else, row = reshape(row,1,[]); end
if numel(row) ~= m || (integer_required && any(row ~= fix(row)))
    error('pde_fit_batched:SampleParameter', '%s must be scalar or contain one valid value per sample.',name);
end
end

function value = on_device(value,prototype)
if isa(prototype,'gpuArray')
    if ~isa(value,'gpuArray'), value = gpuArray(value); end
else
    value = host(value);
end
value = double(value);
end

function value = host(value)
if isa(value,'gpuArray'), value = gather(value); end
end

function options = complete_options(options)
if ~isstruct(options) || ~isscalar(options)
    error('pde_fit_batched:Options', 'options must be a scalar structure.');
end
if isfield(options,'alpha') && ~isfield(options,'omega'), options.omega = options.alpha; end
defaults = struct('step',1,'tau',1.1,'max_steps',100000,'eta',1, ...
    's',1.5,'initial_time',1,'omega',3,'kappa',1,'theta',1.25,'compute_debiased',true);
fields = fieldnames(defaults);
for j = 1:numel(fields)
    key = fields{j};
    if ~isfield(options,key), options.(key) = defaults.(key); end
end
positive = {'step','tau','initial_time'};
for j = 1:numel(positive)
    validateattributes(options.(positive{j}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(options.max_steps,{'numeric'},{'scalar','integer','finite','nonnegative'});
validateattributes(options.eta,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(options.s,{'numeric'},{'scalar','real','finite','>',-0.5});
validateattributes(options.omega,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(options.kappa,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(options.theta,{'numeric'},{'scalar','real','finite','>',0,'<',2});
validateattributes(options.compute_debiased,{'logical','numeric'},{'scalar','real','finite'});
options.compute_debiased = logical(options.compute_debiased);
end
