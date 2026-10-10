function tests = test_pde_source_region
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.oldPath = path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testSourceAreaAndLinearMomentsAtEveryRefinement(testCase)
for rt = 0:2
    problem = pde_problem(rt);
    one = ones(size(problem.pe));
    x = problem.mesh.node(:,1);
    y = problem.mesh.node(:,2);
    verifyEqual(testCase,one.'*problem.M0*one,1,'AbsTol',2e-13);
    verifyEqual(testCase,one.'*problem.M0*problem.pe,1,'AbsTol',2e-13);
    verifyEqual(testCase,x.'*problem.M0*problem.pe,1/12,'AbsTol',2e-13);
    verifyEqual(testCase,y.'*problem.M0*problem.pe,1/12,'AbsTol',2e-13);
    verifyEqual(testCase,problem.pe.'*problem.M0*problem.pe,7/6,'AbsTol',2e-13);
    verifyEqual(testCase,norm(problem.source_mass_root*problem.pe)^2, ...
        7/6,'AbsTol',2e-12);
end
end

function testInterfaceNodesAndForwardLoadAreConsistent(testCase)
problem = pde_problem(2);
nodes = problem.mesh.node;
inside = all(abs(nodes) <= 0.5+1e-13,2);
verifyEqual(testCase,problem.node0,find(inside));
verifyEqual(testCase,numel(problem.node0),753);
interface = inside & any(abs(abs(nodes)-0.5) <= 1e-13,2);
verifyGreaterThan(testCase,nnz(interface),0);
verifyGreaterThan(testCase,full(min(diag(problem.M0(interface,interface)))),0);
verifyEqual(testCase,problem.pe(inside),1+nodes(inside,1)+nodes(inside,2));
verifyEqual(testCase,problem.pe(~inside),zeros(nnz(~inside),1));
u = problem.solve(problem.M0*problem.pe);
verifyEqual(testCase,problem.F*problem.pe,u(problem.boundary),'AbsTol',2e-13);
verifyLessThan(testCase,norm(problem.M0*problem.G-problem.F.'*problem.W,'fro'),1e-12);
end

function testClosedSquareIncludesBoundaryCellsButNotExteriorCells(testCase)
mesh = struct('node',[-0.5,-0.5;0.5,-0.5;0.5,0.5;-0.5,0.5;1,0], ...
    'elem',[1,2,3;1,3,4;2,5,3],'Dirichlet',zeros(0,2));
[~,~,~,M0] = Stiff(0,mesh);
verifyEqual(testCase,full(sum(M0(:))),1,'AbsTol',1e-14);
verifyEqual(testCase,full(M0(5,:)),zeros(1,5));
mesh.node(2,1) = mesh.node(2,1)+eps;
[~,~,~,M0] = Stiff(0,mesh);
verifyEqual(testCase,full(sum(M0(:))),1,'AbsTol',1e-14);
end

function testCurrentExperimentVersionRejectsOldSourceRegion(testCase)
config = pde_experiment_config('Landweber');
verifyEqual(testCase,config.version,'pde-calibrated-bootstrap-v2');
verifyTrue(testCase,endsWith(config.output_root,'paper-table6'));
problem = pde_problem(2);
verifyEqual(testCase,problem.assembly_version,'closed-source-region-v2');
cached = pde_problem(2);
verifyEqual(testCase,cached.assembly_version,problem.assembly_version);
verifyEqual(testCase,cached.M0,problem.M0);
end

function testDefaultNesterovIsStableAndReachesDiscrepancy(testCase)
problem = pde_problem(2);
config = pde_experiment_config('Nesterov');
lambda = max(real(eig(problem.normal_boundary)));
verifyEqual(testCase,config.options.step,1);
verifyEqual(testCase,config.options.omega,5);
verifyLessThanOrEqual(testCase,config.options.step*lambda,1);
noise = problem.Q*(config.sigma*pde_noise(config.seed,1,0,0,numel(problem.boundary)));
noise_norm = sqrt(noise.'*problem.W*noise);
z = problem.F*problem.pe+noise;
fit = pde_fit(problem,z,config.method,config.options,noise_norm);
verifyFalse(testCase,fit.cap_hit);
verifyGreaterThan(testCase,fit.steps,0);
verifyLessThanOrEqual(testCase,fit.residual_norm_original,config.options.tau*noise_norm);
verifyTrue(testCase,all(isfinite(fit.debiased)));
earlier = pde_fit(problem,z,config.method,config.options,noise_norm,fit.steps-1);
verifyGreaterThan(testCase,earlier.residual_norm_original,config.options.tau*noise_norm);
end
