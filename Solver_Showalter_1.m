function [beta,beta1,L2err,L2err1,LinfErr,LinfErr1,Residue,Residue1] = Solver_Showalter_1(rt,tau,dt,norm_e,gD_obs,gN_obs)
% Solver compute approximate solution and its corresponding error,
% mean value as well as its variance
% USAGE
%   [fh,L2err,ct,iternum]=Solver(rt,vap)
%
% Input
%   rt: The times of the refinement of triangulation
%   m: Reflect different noise level
%   tau: A parameter in discrepancy function
%   dt: Time step size
%   s: constant in damping parameter (1+2s)/t
%   alpha: Step size

% output
%   ph: Reconstructed solution
%   L2err: L^2-norm error in ph
%   ct: computional time of cpu cost
%   k(iternum): iterative number

% R.F. Gong 2-20-2019

%--------------------------------------------------------------------------
% Construct finite element mesh 
%--------------------------------------------------------------------------
ct = cputime;  
mesh = Mesh(rt);
n = size(mesh.node,1);

load infor.mat CM C M0 p2norm node0 pe

% right-hand-side for Dirichlet and Neumann boundary condition
g1 = sparse(n,1); % initialize observation g1 on current mesh
g1_exact = sparse(n,1);
g2 = sparse(n,1); % initialize observation g1 on current mesh
g2_exact = sparse(n,1);
for k = 1:n
    if (ismember(k,mesh.Dirichlet(:,1))) % if kth node is a boundary one 
       x = mesh.node(k,1)-gD_obs(:,1);
       y = mesh.node(k,2)-gD_obs(:,2);
       d2 = x.^2 + y.^2;
       [~,ix] = sort(d2);
       g1(k) = gD_obs(ix(1),4); 
       g1_exact(k) = gD_obs(ix(1),3);
       g2(k) = gN_obs(ix(1));
    end 
end

%--------------------------------------------------------------------------
% Set 
%--------------------------------------------------------------------------

tolerr = 1e-6;
maxiter = 200000;
b = C*g1;
evalc('[B] = bicg(CM,b,tolerr,maxiter);');

Cg2 = C*g2;
b = Cg2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');  % [u] = bicg(CM,b,tolerr,maxiter);
diff = u-g1;
residue = sqrt(diff'*C*diff);

p = sparse(n,1);
k = 0;

while residue > tau*norm_e && k < 1e5
    %--------------------------------------------------------------------------
    % Showalte
    %--------------------------------------------------------------------------
    k = k+1;
    b = M0*p+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    k1 = -w+B;
    
    b = M0*(p+dt*k1/2)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    k2 = -w+B;

    b = M0*(p+dt*k2/2)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    k3 = -w+B;
    
    b = M0*(p+dt*k3)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    k4 = -w+B;

    p = p+dt*(k1+2*k2+2*k3+k4)/6;

    b = M0*p+C*g2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    diff = u-g1;
    residue = sqrt(diff'*C*diff);
end
beta1 = p;
Residue1 = residue;
L2err1 = sqrt((beta1-pe)'*M0*(beta1-pe));
[U, S, V] = svds(M0);
S_sqrt  = sqrtm(S);
Ph1       = U * S_sqrt * V'*(beta1-pe);
LinfErr1  = max(abs(Ph1));
i = 0;

while i <= k
%--------------------------------------------------------------------------
% Landweber_debias
%--------------------------------------------------------------------------
    i = i+1;
    b = M0*p+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    k1 = -w+B;
    
    b = M0*(p+dt*k1/2)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    k2 = -w+B;

    b = M0*(p+dt*k2/2)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    k3 = -w+B;
    
    b = M0*(p+dt*k3)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    k4 = -w+B;

    p = p+dt*(k1+2*k2+2*k3+k4)/6;

end
beta = p;
b = M0*beta+C*g2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');
Residue  = sqrt((u-g1)'*C*(u-g1));
L2err    = sqrt((beta-pe)'*M0*(beta-pe));
[U, S, V] = svds(M0);
S_sqrt  = sqrtm(S);
Ph1       = U * S_sqrt * V'*(beta-pe);
LinfErr  = max(abs(Ph1));
ct = cputime -ct;
   

%--------------------------------------------------------------------------
% End of function AcceleratedBLTSolver1
%--------------------------------------------------------------------------