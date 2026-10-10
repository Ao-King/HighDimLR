function [ph,L2RelErr,LinfErr,ct,k,g1,g2] = Solver_ARk_1(rt,tau,dt,kappa,norm_e,gD_obs,gN_obs)
% Compatibility wrapper. norm_e is the effective boundary W-norm of noise.
r = pde_legacy_solver('ARk',rt,tau,dt,norm_e,gD_obs,gN_obs,struct('kappa',kappa,'compute_debiased',false));
ph = r.original;
L2RelErr = r.L2_original;
LinfErr = r.Linf_original;
ct = r.cpu_time;
k = r.steps;
g1 = r.g1;
g2 = r.g2;
end
