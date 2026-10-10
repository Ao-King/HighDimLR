function [gD_obs,gN_obs,u,norm_e,details] = Observe(rt,delta,source,stream)
%OBSERVE Same-noise boundary observations, optionally from a fitted source.
problem = pde_problem(rt);
if nargin < 3 || isempty(source)
    source = problem.pe;
end
if nargin < 4
    stream = [];
end
validateattributes(delta, {'numeric'}, {'real','finite','scalar','nonnegative'});
validateattributes(source, {'numeric'}, {'real','finite','column', ...
    'numel',size(problem.M0,1)});
if isempty(stream)
    e = delta*randn(numel(problem.boundary),1);
else
    e = delta*randn(stream,numel(problem.boundary),1);
end
u = problem.solve(problem.M0*source);
exact_boundary = u(problem.boundary);
gD_obs = [problem.mesh.node(problem.boundary,:), exact_boundary, exact_boundary+e];
gN_obs = e;
effective_noise = problem.Q*e;
norm_e = sqrt(max(0,effective_noise'*problem.W*effective_noise));
details = struct('node_noise',e,'effective_noise',effective_noise, ...
    'z',exact_boundary+effective_noise,'source',source);
end
