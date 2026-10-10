% Legacy export; use the same mesh and source region as the experiment runner.
rt = 2;
problem = pde_problem(rt);
mesh = problem.mesh;
n = size(mesh.node,1);
CM = problem.CM;
C = problem.C;
M0 = problem.M0;
node0 = problem.node0;
pe = problem.pe;
p2norm = sqrt(pe'*M0*pe);

save infor CM C M0 p2norm node0 pe mesh rt
