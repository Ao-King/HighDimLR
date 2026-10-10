function [beta,beta1,L2err,L2err1,LinfErr,LinfErr1,Residue,Residue1] = Solver_HBF_2(rt,tau,dt,eta,norm_e,gD_obs,gN_obs)
% Compatibility wrapper. norm_e is the effective boundary W-norm of noise.
r = pde_legacy_solver('HBF',rt,tau,dt,norm_e,gD_obs,gN_obs,struct('eta',eta));
beta = r.debiased;
beta1 = r.original;
L2err = r.L2_debiased;
L2err1 = r.L2_original;
LinfErr = r.Linf_debiased;
LinfErr1 = r.Linf_original;
Residue = r.residual_norm_debiased;
Residue1 = r.residual_norm_original;
end
