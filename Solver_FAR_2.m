function [beta,beta1,L2err,L2err1,LinfErr,LinfErr1,Residue,Residue1] = Solver_FAR_2(rt,tau,theta,dt,norm_e,gD_obs,gN_obs,bjk,ajk1,ajk2)
% delta = 0.002; rt = 2; tau = 2; dt = 0.1; theta = 1.5; 

[ph,L2err1,LinfErr1,Residue1,k,g1,g2] = Solver_FAR_1(rt,tau,theta,dt,norm_e,gD_obs,gN_obs,bjk,ajk1,ajk2);
beta1 = ph;
mesh = Mesh(rt);
n = size(mesh.node,1);

load infor.mat CM C M0 node0 pe

p = sparse(n,1);
tolerr = 1e-6;
maxiter = 200000;
Cg2 = C*g2;
b = M0*ph+Cg2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');
y = u; % y等价于Yn

b = Cg2;
evalc('[u] = bicg(CM,b,tolerr,maxiter);');

gradient = [];
gradient1 = sparse(n,1);
gradient2 = sparse(n,1);
cb = dt^theta/theta;
ca = dt^theta/(theta*(theta+1));
gamma_theta = gamma(theta);
i = 0;

% % 初始化参数
% N = 1000; 
% 
% % 计算 b_{j,k+1}
% bjk = cb.*((1:N+1).^theta-(0:N).^theta);
% 
% % 计算 d_{j,k+1} 
% ajk1 = ca.*((0:N).^(theta+1)-((0:N) - theta).*(1:N+1).^theta); % j = 0
% ajk2 = ca.*((3:N+1).^(theta+1) + (1:N-1).^(theta+1) ...
%              - 2.*(2:N).^(theta+1)); % 1 <= j <= N-1


while i <= k
    %--------------------------------------------------------------------------
    % Fractional asymptotical regularization
    %--------------------------------------------------------------------------
    i = i+1;
    
    b = C*(y-u);
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    gradient = [gradient,w];
    L = size(gradient,2);
    gradient1 = sum(bjk(1:L).*gradient,2);
    gradient2 = ajk1(i).*gradient(:,1)+sum(ajk2(1:L).*gradient,2);
    pp = gradient1/gamma_theta;
    b = Cg2+M0*pp;
    evalc('[u] = bicg(CM,b,tolerr,maxiter);');
    b = C*(g1-u);
    evalc('[w] = bicg(CM,b,tolerr,maxiter);');
    p = (ca*w+gradient2)/gamma_theta;
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

