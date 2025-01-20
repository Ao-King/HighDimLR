function [ph,L2RelErr,LinfErr,ct,k,g1,g2]=Solver_HBF_1(rt,tau,dt,eta,norm_e,gD_obs,gN_obs) % (2,6,1.01,3,2.5)
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

% right-hand-side for Neumann boundary condition
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
Cg2 = C*g2;
b = Cg2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');  % [u] = bicg(CM,b,tolerr,maxiter);
diff = u-g1;
residue = sqrt(diff'*C*diff);

b = C*g1;
evalc('[B] = bicg(CM,b,tolerr,maxiter);');% Xn'Yn

p = sparse(n,1);
q = sparse(n,1);
k = 0;

while residue > tau*norm_e && k < 1e5
    %--------------------------------------------------------------------------
    % Second order method----Runge-Kutta method
    %--------------------------------------------------------------------------
    k = k+1;
    
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    K11 = q;
    K12 = -w-eta*q+B;
    
    b = M0*(p+dt*K11/2)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    K21 = q+dt*K12/2;
    K22 = -w-eta*K21+B;

    b = M0*(p+dt*K21/2)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    K31 = q+dt*K22/2;
    K32 = -w-eta*K31+B;

    b = M0*(p+dt*K31)+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*u;
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    K41 = q+dt*K32;
    K42 = -w-eta*K41+B;

    p = p+dt*(K11+2*K21+2*K31+K41)/6;
    q = q+dt*(K12+2*K22+2*K32+K42)/6;

    b = M0*p+Cg2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    diff = u-g1;
    residue = sqrt(diff'*C*diff);
end
ph = p;
L2RelErr = sqrt((ph-pe)'*M0*(ph-pe));
[U, S, V] = svds(M0);
S_sqrt  = sqrtm(S);
Ph1       = U * S_sqrt * V'*(ph-pe);
LinfErr  = max(abs(Ph1));
% M0_sqrt = sqrtm(M0);
% Ph1     = M0_sqrt*(ph-pe);
% LinfErr  = max(abs(Ph1));

% x = 0:1:k-1;
% x = x';
% plot(x,L2RelErr)
% axis([1, k 0 max(L2RelErr)])
% node = (1:n)';
% node1 = setdiff(node,node0);
% ph(node1) = 0;
% ph = full(ph);
% pe = full(pe);
% 
% subplot(1,2,1)
% trisurf(mesh.elem, mesh.node(:,1), mesh.node(:,2), pe, ...
%        'FaceColor', 'interp', 'EdgeColor', 'interp');
% title('The exact source','FontSize', 10)
% view(2); axis equal
% colorbar
% axis([-1 1 -1 1])
% xlabel('x')
% ylabel('y')
% set(gca,'CLim',[0,2])
% 
% subplot(1,2,2)
% trisurf(mesh.elem, mesh.node(:,1), mesh.node(:,2), ph, ...
%         'FaceColor', 'interp', 'EdgeColor', 'interp');
% title('The approximate solution','FontSize', 10)
% view(2); axis equal
% colorbar
% axis([-1 1 -1 1])
% xlabel('x')
% ylabel('y')
% set(gca,'CLim',[0,2])
ct = cputime -ct;
   

%--------------------------------------------------------------------------
% End of function AcceleratedBLTSolver1
%--------------------------------------------------------------------------