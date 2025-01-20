function [beta,beta1,L2err,L2err1,LinfErr,LinfErr1,Residue,Residue1] = Solver_ARk_2(rt,tau,dt,kappa,norm_e,gD_obs,gN_obs) 
% delta = 0.002; rt = 2; tau = 2; dt = 1.6; kappa = 0.1;

[ph,L2err1,LinfErr1,~,k,g1,g2]=Solver_ARk_1(rt,tau,dt,kappa,norm_e,gD_obs,gN_obs);
beta1 = ph;
mesh = Mesh(rt);
n = size(mesh.node,1);

load infor.mat CM C M0 node0 pe

p = sparse(n,1);
v = sparse(n,1);
tolerr = 1e-6;
maxiter = 200000;
Cg2 = C*g2;
b = M0*ph+Cg2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);'); % [u] = bicg(CM,b,tolerr,maxiter);
Residue1 = sqrt((u-g1)'*C*(u-g1));
y = u; % y等价于Yn
i = 0;
while i <= k
    %--------------------------------------------------------------------------
    % AR^Kappa method
    %--------------------------------------------------------------------------
    i = i+1;
    p = p+dt*v;
    i_kappa = i^(-kappa);
    tk = dt*i_kappa;
    tk_kappa = tk^kappa;
    alpha = 1+((tk_kappa)^(-1)-kappa)/k;
    b = M0*p+Cg2; 
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    diff = y-u; 
    b = C*diff;
    evalc('[g] = bicg(CM,b,tolerr,maxiter);');
    Y = C*(v+dt*g)./(dt*tk_kappa);

    b = M0*v+Cg2; 
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');

    b = C*(u-Y);
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');

    v = -dt*tk_kappa*w./alpha;
end
beta2 = p;
beta = 2*beta1-beta2;
b = M0*beta+Cg2;
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
