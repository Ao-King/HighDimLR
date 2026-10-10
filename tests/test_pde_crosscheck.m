function tests = test_pde_crosscheck
% Independent checks of matrix evolution, first-hit stopping, and noise scale.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.oldPath = path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testMatrixRKPolynomials(testCase)
F = [0.5,-0.2,0.1; 0.1,0.4,-0.3];
W = [1,0.1;0.1,2];
G = diag([1,0.5,2])*F.'*W;
problem = struct('F',F,'G',G,'W',W,'d0',trace(W));
z = [0.8;-0.2];
h = 0.12;
steps = 5;
n = size(F,2);
A = [-G*F,G*z;zeros(1,n+1)];
R = rk4_polynomial(h*A);
expected = R^steps*[zeros(n,1);1];
fit = pde_fit(problem,z,'Showalter',struct('step',h),0,steps);
verifyEqual(testCase,fit.original,expected(1:n),'AbsTol',1e-13);
eta = 0.7;
A = [zeros(n),eye(n),zeros(n,1); ...
    -G*F,-eta*eye(n),G*z;zeros(1,2*n+1)];
R = rk4_polynomial(h*A);
expected = R^steps*[zeros(2*n,1);1];
fit = pde_fit(problem,z,'HBF',struct('step',h,'eta',eta),0,steps);
verifyEqual(testCase,fit.original,expected(1:n),'AbsTol',1e-13);
end

function testEveryMethodUsesFirstHitOfReturnedState(testCase)
problem = struct('F',1,'G',1,'W',1,'d0',1);
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
options = struct('step',0.1,'tau',1.1,'max_steps',200,'compute_debiased',false);
for j = 1:numel(methods)
    selected = pde_fit(problem,1,methods{j},options,0.2);
    verifyFalse(testCase,selected.cap_hit);
    verifyLessThanOrEqual(testCase,selected.residual_norm_original,0.22+1e-14);
    for k = 0:selected.steps-1
        earlier = pde_fit(problem,1,methods{j},options,0.2,k);
        verifyGreaterThan(testCase,earlier.residual_norm_original,0.22);
    end
end
end

function testActualMeshVarianceAndNoiseEnergy(testCase)
problem = pde_problem(0);
noise = 0.002*cos((1:numel(problem.boundary)).');
z = problem.F*problem.pe+problem.Q*noise;
fit = pde_fit(problem,z,'Nesterov',struct('step',0.1),0,4);
r_original = z-problem.F*fit.original;
r_debiased = z-problem.F*fit.debiased;
verifyEqual(testCase,fit.sigma2_original, ...
    (r_original.'*problem.W*r_original)/problem.d0,'RelTol',1e-12);
verifyEqual(testCase,fit.sigma2_debiased, ...
    (r_debiased.'*problem.W*r_debiased)/problem.d0,'RelTol',1e-12);
total = 0;
for j = 1:numel(problem.boundary)
    q = problem.Q(:,j);
    total = total+q.'*problem.W*q;
end
verifyEqual(testCase,problem.d0,total,'RelTol',1e-13);
verifyGreaterThan(testCase,norm(problem.Q-eye(size(problem.Q)),'fro'),0.01);
end

function R = rk4_polynomial(A)
R = eye(size(A))+A+A^2/2+A^3/6+A^4/24;
end
