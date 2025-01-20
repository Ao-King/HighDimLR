function [beta,beta1,L2err,L2err1,LinfErr,LinfErr1,Residue,Residue1] = Solver_Nesterov_2(rt,tau,alpha,omega,norm_e,gD_obs,gN_obs) % (2,1.01,3,2.5,0.02)
% Nesterov

[ph,L2err1,LinfErr1,~,k,g1,g2]=Solver_Nesterov_1(rt,tau,alpha,omega,norm_e,gD_obs,gN_obs);
beta1 = ph;
mesh = Mesh(rt);
n = size(mesh.node,1);

load infor.mat CM C M0 node0 pe

p0 = sparse(n,1);
p1 = sparse(n,1);
i = 0;
b = M0*ph+C*g2;
tolerr = 1e-6;
maxiter = 200000;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');
% [u] = bicg(CM,b,tolerr,maxiter);
Residue1 = sqrt((u-g1)'*C*(u-g1));
y = u; % y等价于Yn
while i <= k
%--------------------------------------------------------------------------
% Nesterov
%--------------------------------------------------------------------------
    i = i+1;
    constant = (i-1)/(i+alpha-1);
    z = p1+constant*(p1-p0);

    b = M0*z+C*g2;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    % [u] = bicg(CM,b,tolerr,maxiter);
    b = C*(u-y);
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    % [w] = bicg(CM,b,tolerr,maxiter);

    p = z-omega*w;
    p0 = p1; p1 = p;
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
