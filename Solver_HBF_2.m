function [beta,beta1,L2err,L2err1,LinfErr,LinfErr1,Residue,Residue1] = Solver_HBF_2(rt,tau,dt,eta,norm_e,gD_obs,gN_obs) % (2,1.01,3,2.5,0.02)
% Second order aymptotic
% delta = 0.002; dt = 0.9; tau = 1.1; eta = 0.2;

[ph,L2err1,LinfErr1,~,k,g1,g2]=Solver_HBF_1(rt,tau,dt,eta,norm_e,gD_obs,gN_obs);
beta1 = ph;
mesh = Mesh(rt);
n = size(mesh.node,1);

load infor.mat CM C M0 node0 pe

p = sparse(n,1);
q = sparse(n,1);
i = 0;
tolerr = 1e-6;
maxiter = 200000;
Cg2 = C*g2;
b = M0*ph+Cg2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);'); % [u] = bicg(CM,b,tolerr,maxiter);
Residue1 = sqrt((u-g1)'*C*(u-g1));

b = C*u;
evalc('[B] = bicg(CM,b,tolerr,maxiter);');% Xn'Xnbeta_k

while i <= k
%--------------------------------------------------------------------------
    % Second order method----Runge-Kutta method
    %--------------------------------------------------------------------------
    i = i+1;
    
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
end
beta2 = p;
beta = 2*beta1-beta2;
b = M0*beta+C*g2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');
Residue  = sqrt((u-g1)'*C*(u-g1));
L2err    = sqrt((beta-pe)'*M0*(beta-pe));
[U, S, V] = svds(M0);
S_sqrt  = sqrtm(S);
Ph1       = U * S_sqrt * V'*(beta-pe);
LinfErr  = max(abs(Ph1));

% figure(2)
% node = (1:n)';
% node1 = setdiff(node,node0);
% beta(node1) = 0;
% beta = full(beta);
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
% trisurf(mesh.elem, mesh.node(:,1), mesh.node(:,2), beta, ...
%         'FaceColor', 'interp', 'EdgeColor', 'interp');
% title('The approximate solution','FontSize', 10)
% view(2); axis equal
% colorbar
% axis([-1 1 -1 1])
% xlabel('x')
% ylabel('y')
% set(gca,'CLim',[0,2])
