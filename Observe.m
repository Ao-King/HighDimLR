 function [gD_obs,gN_obs,u,norm_e] = Observe(rt,delta)
% Observe compute observation from exact solution obtained by solving governning
% equation in a rather small meshsize domain
% USAGE
% [gD_obs,u] = Observe(rt)
%
% Input
%   rt: the times of the refinement of triangulation
% output
% gD_obs: Dirichlet observation on boundary \Gamma=\partial\Omega
% u: the finite element solution corresponding to the exact source function

% R.F. Gong 8-23-2018

%--------------------------------------------------------------------------
% Construct finite element mesh 
%--------------------------------------------------------------------------
mesh = Mesh(rt);
n = size(mesh.node,1);
        
%--------------------------------------------------------------------------
% Exact light source function p
%--------------------------------------------------------------------------
% find out elements in \Omega_0: -0.5<x,y<0.5 where f = 1+x+y             
indexe = find(mesh.node(mesh.elem(:,1),1)>-0.5 & mesh.node(mesh.elem(:,1),1)<0.5 ...
            & mesh.node(mesh.elem(:,1),2)>-0.5 & mesh.node(mesh.elem(:,1),2)<0.5 ...
            & mesh.node(mesh.elem(:,2),1)>-0.5 & mesh.node(mesh.elem(:,2),1)<0.5 ...
            & mesh.node(mesh.elem(:,2),2)>-0.5 & mesh.node(mesh.elem(:,2),2)<0.5 ...
            & mesh.node(mesh.elem(:,3),1)>-0.5 & mesh.node(mesh.elem(:,3),1)<0.5 ...
            & mesh.node(mesh.elem(:,3),2)>-0.5 & mesh.node(mesh.elem(:,3),2)<0.5);
node = [mesh.elem(indexe,1);mesh.elem(indexe,2);mesh.elem(indexe,3)];         
node = unique(node);
p = sparse(n,1);
% Set true light source p in \Omega_0: p = 1+x+y
p(node)=SetSource(mesh.node(node,1),mesh.node(node,2));

%--------------------------------------------------------------------------
% Solve the governing PDE for the exact light source function f
%--------------------------------------------------------------------------

%Establish linear system
[K,M,~,M0] = Stiff(rt);%Assemble coefficient matrix

% Assemble right-hand-side
b = M0*p;

CM = K+M;

% Solve the governing PDE
tolerr = 1e-8;
maxiter = 20000;
[u,~,~,~] = bicg(CM,b,tolerr,maxiter);


% trisurf(mesh.elem, mesh.node(:,1), mesh.node(:,2), u , ...
% 'FaceColor', 'interp', 'EdgeColor', 'interp');
% view(2)
%--------------------------------------------------------------------------
% Compute observation g1 and g2 on the boundary
%--------------------------------------------------------------------------
gD_obs(:,1) = mesh.node(mesh.Dirichlet(:,1),1);
gD_obs(:,2) = mesh.node(mesh.Dirichlet(:,1),2);
gD_obs(:,3) = u(mesh.Dirichlet(:,1));

noise = randn(size(gD_obs,1),1);
e = delta*noise;
norm_e = norm(e);
% norm_noise = norm(noise);
% e = delta*noise/norm_noise;
gD_obs(:,4) = gD_obs(:,3)+e;
gN_obs = e;

% save gD_obs in file Obs1_5.mat, where number 5 denotes 
% the number of refinement of mesh 
function y = SetSource(x,y)% Source function p
        y = 1+x+y;

%--------------------------------------------------------------------------
% End of function Observe
%--------------------------------------------------------------------------