% FAR Table 6 experiment; set pde_overrides for a smaller run.
if ~exist('pde_overrides','var')
    pde_overrides = struct();
end
pde_result = run_paper_table6('FAR',pde_overrides);
source_error_summary = pde_result.source_error_summary;
M0 = pde_result.figure_data.M0;
mesh = pde_result.figure_data.mesh;
pe = pde_result.figure_data.pe;
node0 = pde_result.figure_data.node0;
n = pde_result.figure_data.n;
disp(pde_result.summary);
