function result = pde_legacy_solver(method, rt, tau, step, noise_norm, gD_obs, gN_obs, parameters)
%PDE_LEGACY_SOLVER Adapt historical solver inputs to the shared boundary map.
started = cputime;
problem = pde_problem(rt);
if size(gD_obs, 2) < 4 || size(gD_obs, 1) ~= numel(gN_obs)
    error('pde_legacy_solver:Observations', 'Boundary observations have incompatible dimensions.');
end
boundary_values = zeros(numel(problem.boundary), 1);
neumann_values = boundary_values;
for j = 1:numel(problem.boundary)
    xy = problem.mesh.node(problem.boundary(j), :);
    distances = sum(bsxfun(@minus, gD_obs(:, 1:2), xy).^2, 2);
    [distance, index] = min(distances);
    if distance > 1e-16
        error('pde_legacy_solver:MeshMismatch', 'Boundary observation coordinates do not match the solver mesh.');
    end
    boundary_values(j) = gD_obs(index, 4);
    neumann_values(j) = gN_obs(index);
end
parameters.step = step;
parameters.tau = tau;
z = boundary_values - problem.B * neumann_values;
result = pde_fit(problem, z, method, parameters, noise_norm);
if result.cap_hit
    warning('pde_legacy_solver:CapHit', 'The iteration cap was reached before the discrepancy criterion.');
end
result.g1 = zeros(size(problem.CM, 1), 1);
result.g2 = result.g1;
result.g1(problem.boundary) = boundary_values;
result.g2(problem.boundary) = neumann_values;
error_original = problem.source_mass_root * (result.original - problem.pe);
result.L2_original = norm(error_original);
result.Linf_original = norm(error_original, inf);
result.L2_debiased = NaN;
result.Linf_debiased = NaN;
if ~isempty(result.debiased)
    error_debiased = problem.source_mass_root * (result.debiased - problem.pe);
    result.L2_debiased = norm(error_debiased);
    result.Linf_debiased = norm(error_debiased, inf);
end
result.cpu_time = cputime - started;
end
