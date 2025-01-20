parpool; % 启动默认的并行池
delta = 0.002;
rt    = 2;
N     = 1000;
B     = 500;
tau   = 1.1;
omega = 2.125;
% omega = 2.5 100 I1=95;
% omega = 2.25 98 I1=3; 996 17
% omega = 2.125 98 I1=0; 996 0
% omega = 2.375 100 I1=42
I      = 0;
I1     = 0;
Inf    = 0;
Inf1   = 0;
Err    = 0;
Err1   = 0;
Sigma  = 0;
Sigma1 = 0; 

for i = 1:N
    i
    rng(i)
    inf1  = [];
    inf   = [];
    [gD_obs,gN_obs,u,norm_e] = Observe(rt,delta);
    [beta,beta1,L2err,L2err1,Linferr,Linferr1,Residue,Residue1] = Solver_Landweber_1(rt,tau,omega,norm_e,gD_obs,gN_obs);
    parfor j = 1:B
        [LInfErr1,LInfErr]=compute_Landweber(j,M0,rt,Residue,Residue1,tau,omega,beta,beta1);
        inf1 = [inf1,LInfErr1];
        inf  = [inf,LInfErr];
    end
    q1_95  = quantile(inf1, 0.95);
    q_95   = quantile(inf, 0.95);
    I1     = I1 + 1 * (Linferr1 <= q1_95);
    I      = I + 1 * (Linferr <= q_95);
    Inf1   = Inf1 + Linferr1;
    Inf    = Inf + Linferr;
    Err1   = Err1 + L2err1;
    Err    = Err + L2err;
    Sigma  = Sigma + abs(Residue^2-0.002^2);
    Sigma1 = Sigma1 + abs(Residue1^2-0.002^2);
end

averInf    = Inf/N;
averInf1   = Inf1/N;
averErr    = Err/N;
averErr1   = Err1/N;
averSigma  = Sigma/N;
averSigma1 = Sigma1/N;

function [LInfErr1,LInfErr]=compute_Landweber(seed,M0,rt,Residue,Residue1,tau,omega,beta,beta1)
rng(seed)
[U, S, V] = svds(M0);
S_sqrt  = sqrtm(S);
[gD_obs,gN_obs,~,norm_e] = Observe(rt,Residue);
[gD_obs1,gN_obs1,~,norm_e1] = Observe(rt,Residue1);
[~,beta1star,~,~,~,~,~,~] = Solver_Landweber_1(rt,tau,omega,norm_e1,gD_obs1,gN_obs1);
Ph1        = U * S_sqrt * V'*(beta1star-beta1);
LInfErr1   = max(abs(Ph1));
[betastar,~,~,~,~,~,~,~] = Solver_Landweber_1(rt,tau,omega,norm_e,gD_obs,gN_obs);
Ph         = U * S_sqrt * V'*(betastar-beta);
LInfErr    = max(abs(Ph));
end
