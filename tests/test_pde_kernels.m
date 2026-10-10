function tests = test_pde_kernels
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.oldPath = path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.problem = pde_problem(0);
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testZeroStepsAndAllMethods(testCase)
problem = scalar_problem;
methods = method_names;
for j = 1:numel(methods)
    r = pde_fit(problem, 1, methods{j}, struct('step', 0.1), 1);
    verifyEqual(testCase, r.steps, 0);
    verifyEqual(testCase, r.original, 0);
    verifyEqual(testCase, r.debiased, 0);
    verifyFalse(testCase, r.cap_hit);
end
end

function testFixedStepDoubleFitAgainstExplicitLinearMap(testCase)
problem = scalar_problem;
problem.F = [0.8, 0.2; -0.1, 0.6];
problem.G = problem.F.';
problem.W = eye(2);
problem.d0 = 2;
z = [1; -0.4];
methods = method_names;
for j = 1:numel(methods)
    for k = [0, 1, 3, 6]
        opt = struct('step', 0.1, 'compute_debiased', false, 'max_steps', 1);
        L = zeros(2);
        for col = 1:2
            basis = zeros(2, 1);
            basis(col) = 1;
            fit = pde_fit(problem, basis, methods{j}, opt, 1e6, k);
            L(:, col) = fit.original;
            verifyEqual(testCase, fit.steps, k);
        end
        opt.compute_debiased = true;
        fit = pde_fit(problem, z, methods{j}, opt, 1e6, k);
        verifyEqual(testCase, fit.original, L*z, 'AbsTol', 2e-12);
        verifyEqual(testCase, fit.debiased, (2*L-L*problem.F*L)*z, 'AbsTol', 2e-12);
        verifyFalse(testCase, fit.cap_hit);
    end
end
end

function testScalarLandweberShowalterNesterov(testCase)
problem = scalar_problem;
opt = struct('step', 0.1, 'omega', 3);
fit = pde_fit(problem, 1, 'Landweber', opt, 0, 3);
verifyEqual(testCase, fit.original, 1-0.9^3, 'AbsTol', 1e-14);
verifyEqual(testCase, fit.debiased, 1-0.9^6, 'AbsTol', 1e-14);
fit = pde_fit(problem, 1, 'Showalter', opt, 0, 3);
rk = 1-0.1+0.1^2/2-0.1^3/6+0.1^4/24;
verifyEqual(testCase, fit.original, 1-rk^3, 'AbsTol', 1e-14);
fit = pde_fit(problem, 1, 'Nesterov', opt, 0, 3);
verifyEqual(testCase, fit.original, 0.2872, 'AbsTol', 1e-14);
end

function testHBFRungeKuttaMatrixPolynomial(testCase)
problem = scalar_problem;
h = 0.1;
eta = 0.7;
A = [0, 1, 0; -1, -eta, 1; 0, 0, 0];
R = eye(3)+h*A+(h*A)^2/2+(h*A)^3/6+(h*A)^4/24;
expected = R^4 * [0; 0; 1];
fit = pde_fit(problem, 1, 'HBF', struct('step', h, 'eta', eta), 0, 4);
verifyEqual(testCase, fit.original, expected(1), 'AbsTol', 1e-14);
end

function testSOARShiftedVerlet(testCase)
problem = scalar_problem;
h = 0.1;
p = 0;
v = 0;
for i = 1:4
    t = 1+(i-1)*h;
    half = (v+h*(1-p)/2)/(1+2*h/t);
    p = p+h*half;
    v = half+h*(1-p-4*half/(t+h))/2;
end
fit = pde_fit(problem, 1, 'SOAR', struct('step', h, 's', 1.5, 'initial_time', 1), 0, 4);
verifyEqual(testCase, fit.original, p, 'AbsTol', 1e-14);
end

function testFARDirectHistoryFormula(testCase)
problem = scalar_problem;
problem.F = 0.5;
theta = 1.25;
h = 0.025;
p = zeros(1, 5);
for k = 0:3
    predicted = 0;
    accumulated = 0;
    for j = 0:k
        b = h^theta/theta*((k-j+1)^theta-(k-j)^theta);
        if j == 0
            d = k^(theta+1)-(k-theta)*(k+1)^theta;
        else
            d = (k-j+2)^(theta+1)+(k-j)^(theta+1)-2*(k-j+1)^(theta+1);
        end
        predicted = predicted+b*(1-0.5*p(j+1))/gamma(theta);
        accumulated = accumulated+d*(1-0.5*p(j+1));
    end
    p(k+2) = h^theta/(theta*(theta+1)*gamma(theta))*(1-0.5*predicted+accumulated);
    fit = pde_fit(problem, 1, 'FAR', struct('step', h, 'theta', theta), 0, k+1);
    verifyEqual(testCase, fit.original, p(k+2), 'AbsTol', 1e-14);
end
verifyEqual(testCase, p(3), 0.0207803458632375, 'AbsTol', 1e-14);
end

function testFARBeyondLegacyHistoryLimit(testCase)
problem = scalar_problem;
problem.G = 0;
opt = struct('step', 0.01, 'compute_debiased', false);
fit = pde_fit(problem, 1, 'FAR', opt, 0, 3002);
verifyEqual(testCase, fit.steps, 3002);
verifyEqual(testCase, fit.original, 0);
end

function testARKWoodburyMatchesFullImplicitSolve(testCase)
problem = scalar_problem;
problem.F = [0.7, 0.1, 0; -0.2, 0.5, 0];
problem.G = problem.F.';
problem.W = eye(2);
problem.d0 = 2;
z = [1; -0.3];
h = 0.725;
kappa = 1;
p = zeros(3, 1);
v = p;
for i = 1:5
    p = p+h*v;
    t = h*i;
    alpha = 1+(t^(-kappa)-kappa)/i;
    verifyGreaterThan(testCase, alpha, 0);
    v = (alpha*eye(3)+h*t^kappa*problem.G*problem.F) \ (v+h*problem.G*(z-problem.F*p));
end
fit = pde_fit(problem, z, 'ARk', struct('step', h, 'kappa', kappa), 0, 5);
verifyEqual(testCase, fit.original, p, 'AbsTol', 1e-13);
verifyEqual(testCase, fit.original(3), 0);
end

function testStoppingUsesReturnedSource(testCase)
problem = scalar_problem;
fit = pde_fit(problem, 1, 'Nesterov', struct('step', 0.1, 'tau', 1, 'max_steps', 20), 0.85);
verifyEqual(testCase, fit.steps, 2);
verifyEqual(testCase, fit.original, 0.19, 'AbsTol', 1e-14);
verifyEqual(testCase, fit.residual_norm_original, 0.81, 'AbsTol', 1e-14);
end

function testCapsAndInvalidStates(testCase)
problem = scalar_problem;
fit = pde_fit(problem, 1, 'Landweber', struct('step', 0.1, 'max_steps', 1), 0);
verifyTrue(testCase, fit.cap_hit);
verifyEqual(testCase, fit.steps, 1);
fit = pde_fit(problem, 1, 'Landweber', struct('max_steps', 0), 0);
verifyTrue(testCase, fit.cap_hit);
verifyEqual(testCase, fit.original, 0);
verifyError(testCase, @() pde_fit(problem, NaN, 'Landweber', struct, 0), 'MATLAB:expectedFinite');
verifyError(testCase, @() pde_fit(problem, 1, 'unknown', struct, 0), 'pde_fit:Method');
verifyError(testCase, @() pde_fit(problem, 1, 'ARk', struct('step', 1, 'kappa', 2), 0, 1), 'pde_fit:ARKCoefficient');
problem.G = realmax;
verifyError(testCase, @() pde_fit(problem, 1, 'Landweber', struct('step', realmax), 0, 1), 'pde_fit:NonFiniteIterate');
end

function testActualProblemForwardGradientAndEquivariance(testCase)
problem = testCase.TestData.problem;
n = size(problem.CM, 1);
f = sin((1:n).');
g = cos((1:numel(problem.boundary)).');
u = problem.solve(problem.M0*f+problem.C(:, problem.boundary)*g);
verifyEqual(testCase, u(problem.boundary), problem.F*f+problem.B*g, 'AbsTol', 1e-11);
r = sin(2*(1:numel(problem.boundary)).');
verifyEqual(testCase, f.'*problem.M0*(problem.G*r), (problem.F*f).'*problem.W*r, 'AbsTol', 1e-11);
h = cos(3*(1:numel(problem.boundary)).');
gD = problem.F*f+problem.B*g;
z1 = gD-problem.B*g;
z2 = (gD+problem.B*h)-problem.B*(g+h);
opt = struct('step', 0.1);
methods = method_names;
for j = 1:numel(methods)
    fit1 = pde_fit(problem, z1, methods{j}, opt, 0, 3);
    fit2 = pde_fit(problem, z2, methods{j}, opt, 0, 3);
    verifyEqual(testCase, fit1.original, fit2.original, 'AbsTol', 1e-11);
    verifyEqual(testCase, fit1.debiased, fit2.debiased, 'AbsTol', 1e-11);
end
outside = setdiff(1:n, problem.node0);
f(:) = 0;
f(outside) = 1;
verifyEqual(testCase, norm(problem.F*f), 0);
verifyGreaterThan(testCase, problem.d0, 0);
end

function testAllLegacyOutputContractsAtZeroSteps(testCase)
problem = testCase.TestData.problem;
rt = problem.rt;
gD = [problem.mesh.node(problem.boundary, :), problem.F*problem.pe, problem.F*problem.pe];
gN = zeros(numel(problem.boundary), 1);
common = {rt, 1.1};
combined = {
    'Solver_Landweber_1', [common, {0.1, 1e6, gD, gN}];
    'Solver_Showalter_1', [common, {0.1, 1e6, gD, gN}];
    'Solver_HBF_2', [common, {0.1, 1, 1e6, gD, gN}];
    'Solver_SOAR_2', [common, {0.1, 1e6, gD, gN}];
    'Solver_Nesterov_2', [common, {3, 0.1, 1e6, gD, gN}];
    'Solver_ARk_2', [common, {0.1, 1, 1e6, gD, gN}];
    'Solver_FAR_2', [common, {1.25, 0.1, 1e6, gD, gN, [], [], []}]
    };
for j = 1:size(combined, 1)
    args = combined{j, 2};
    [beta, beta1, L2, L21, inferr, inferr1, residual, residual1] = feval(combined{j, 1}, args{:});
    verifyEqual(testCase, norm(beta), 0);
    verifyEqual(testCase, norm(beta1), 0);
    verifyEqual(testCase, [L2, inferr, residual], [L21, inferr1, residual1], 'AbsTol', 1e-13);
end
ordinary = {
    'Solver_HBF_1', [common, {0.1, 1, 1e6, gD, gN}];
    'Solver_SOAR_1', [common, {0.1, 1e6, gD, gN}];
    'Solver_Nesterov_1', [common, {3, 0.1, 1e6, gD, gN}];
    'Solver_ARk_1', [common, {0.1, 1, 1e6, gD, gN}];
    'Solver_FAR_1', [common, {1.25, 0.1, 1e6, gD, gN, [], [], []}]
    };
for j = 1:size(ordinary, 1)
    args = ordinary{j, 2};
    [ph, L2, inferr, fourth, k, g1, g2] = feval(ordinary{j, 1}, args{:});
    verifyEqual(testCase, norm(ph), 0);
    verifyEqual(testCase, k, 0);
    verifyGreaterThanOrEqual(testCase, [L2, inferr, fourth], [0, 0, 0]);
    verifyEqual(testCase, g1(problem.boundary), gD(:, 4), 'AbsTol', 1e-13);
    verifyEqual(testCase, norm(g2), 0);
end
end

function testProblemAssemblyPreservesRandomState(testCase)
before = rng;
problem = pde_problem(1);
verifyEqual(testCase, rng, before);
verifyEqual(testCase, problem.mesh_seed, 0);
end

function problem = scalar_problem
problem = struct('F', 1, 'G', 1, 'W', 1, 'd0', 1);
end

function names = method_names
names = {'Landweber', 'Showalter', 'HBF', 'SOAR', 'Nesterov', 'FAR', 'ARk'};
end
