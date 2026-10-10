function problem = pde_problem(rt)
%PDE_PROBLEM Assemble a common mesh, forward map, and boundary noise metric.
persistent cachedRefinement cachedProblem
assemblyVersion = 'closed-source-region-v2';
validateattributes(rt, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'});
if ~isempty(cachedProblem) && isequal(rt, cachedRefinement) && ...
        isfield(cachedProblem, 'assembly_version') && ...
        strcmp(cachedProblem.assembly_version, assemblyVersion)
    problem = cachedProblem;
    return
end
previousRandomState = rng;
restoreRandomState = onCleanup(@() rng(previousRandomState)); %#ok<NASGU>
rng(0, 'twister');
mesh = Mesh(rt);
[K, M, C, M0] = Stiff(rt, mesh);
CM = K + M;
factor = decomposition(CM, 'chol');
solve = @(rhs) factor \ rhs;
n = size(CM, 1);
boundary = mesh.Dirichlet(:, 1);
if numel(unique(boundary)) ~= numel(boundary)
    error('pde_problem:Boundary', 'Boundary nodes must occur exactly once.');
end
node0 = find(any(M0 ~= 0, 2));
W = full(C(boundary, boundary));
G = full(solve(C(:, boundary)));
B = G(boundary, :);
Q = eye(numel(boundary)) - B;
F = zeros(numel(boundary), n);
sourceSolutions = solve(M0(:, node0));
F(:, node0) = sourceSolutions(boundary, :);
d0 = trace(Q.' * W * Q);
if ~isfinite(d0) || d0 <= 0
    error('pde_problem:NoiseScale', 'The effective boundary noise energy must be positive.');
end
pe = zeros(n, 1);
pe(node0) = 1 + mesh.node(node0, 1) + mesh.node(node0, 2);
problem = struct('rt', rt, 'mesh', mesh, 'CM', CM, 'C', C, 'M0', M0, ...
    'pe', pe, 'boundary', boundary, 'node0', node0, 'W', W, 'B', B, ...
    'Q', Q, 'd0', d0, 'F', F, 'G', G, 'solve', solve, ...
    'source_mass_root', source_mass_sqrt(M0));
problem.normal_boundary = F * G;
problem.noise_model = 'shared-dirichlet-neumann';
problem.mesh_seed = 0;
problem.assembly_version = assemblyVersion;
cachedRefinement = rt;
cachedProblem = problem;
end
