function [ph,L2RelErr,LinfErr,residue,k,g1,g2] = Solver_FAR_1(rt,tau,theta,dt,norm_e,gD_obs,gN_obs,bjk,ajk1,ajk2)
% Compatibility wrapper. norm_e is the effective boundary W-norm of noise.
% Legacy weight arrays are accepted; the shared kernel derives correct weights.
r = pde_legacy_solver('FAR',rt,tau,dt,norm_e,gD_obs,gN_obs,struct('theta',theta,'compute_debiased',false));
ph = r.original;
L2RelErr = r.L2_original;
LinfErr = r.Linf_original;
residue = r.residual_norm_original;
k = r.steps;
g1 = r.g1;
g2 = r.g2;
end
