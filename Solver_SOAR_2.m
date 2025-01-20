function [beta,beta1,L2err,L2err1,LinfErr,LinfErr1,Residue,Residue1]=Solver_SOAR_2(rt,tau,dt,norm_e,gD_obs,gN_obs) % (2,1.01,1,0.02,1)
% SOAR
%rt = 2; dt = 1;
[ph,L2err1,LinfErr1,~,k,g1,g2]=Solver_SOAR_1(rt,tau,dt,norm_e,gD_obs,gN_obs); % ph等价betak 
beta1 = ph;
mesh = Mesh(rt);
n = size(mesh.node,1);

load infor.mat CM C M0 node0 pe

i = 1;
b = M0*ph+C*g2;
tolerr = 1e-6;
maxiter = 200000;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');
Residue1 = sqrt((u-g1)'*C*(u-g1));
y = u; % y等价于Yn

t = 1;
q = sparse(n,1);
p = sparse(n,1);
eta = 4/t;
b = -C*y;
evalc('[w] = bicg(CM,b,tolerr,maxiter);');

while i <= k
    %--------------------------------------------------------------------------
    % SOAR
    %--------------------------------------------------------------------------
    q = (q-dt*w/2)/(1+dt*eta/2);
    p = p+dt*q;

    b = M0*p+C*g2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*(u-y);
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');

    t = t+dt;
    eta = 4/t;
    q = q-dt*(eta*q+w)/2;

    i = i+1;
end
beta2 = p;
beta = 2*beta1-beta2;
b = M0*beta+C*g2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');
Residue = sqrt((u-g1)'*C*(u-g1));
L2err = sqrt((beta-pe)'*M0*(beta-pe));
[U, S, V] = svds(M0);
S_sqrt  = sqrtm(S);
Ph1       = U * S_sqrt * V'*(beta-pe);
LinfErr  = max(abs(Ph1));

% figure(2)
% 
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