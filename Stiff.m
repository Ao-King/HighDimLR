function [K,M,C,M0] = Stiff(rt)
% Stiff assemble stiff and mass matrices
% USAGE
%   [K,M,C,M0] = Stiff(rt)
%
% Input
%   rt: The times of the refinement of triangulation
% Output
%   K: Stiff matrix by FEM
%   M: Mass matrix
%   C: Stiff corresponding to Neumann boundary condition
%   M0: Mass matrix in \Omega_0

%   R.F. Gong 1-28-2019

%--------------------------------------------------------------------------
% Initialize the data 
%--------------------------------------------------------------------------
mesh = Mesh(rt);
n = size(mesh.node,1);
K = sparse(n,n); % stiff matrix
C = sparse(n,n); % stiff corresponding to Neumann boundary condition

%--------------------------------------------------------------------------
% Compute vedge: edge as a vector, area and midpoint of each element
%--------------------------------------------------------------------------
ve(:,:,1) = mesh.node(mesh.elem(:,3),:)-mesh.node(mesh.elem(:,2),:);
ve(:,:,2) = mesh.node(mesh.elem(:,1),:)-mesh.node(mesh.elem(:,3),:);
ve(:,:,3) = mesh.node(mesh.elem(:,2),:)-mesh.node(mesh.elem(:,1),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));

%--------------------------------------------------------------------------
% Assemble Stiffness matrix A, Mass matrix M and then coefficient matrix K
%--------------------------------------------------------------------------
% Assemble Stiffness matrix A
for i = 1:3
    for j = 1:3
        % Aij = \int_T grad u_i * grad u_j
        Kij = (ve(:,1,i).*ve(:,1,j)+ve(:,2,i).*ve(:,2,j))./(4*area); %local stiff matrix
        K = K + sparse(mesh.elem(:,i),mesh.elem(:,j),Kij,n,n);
    end
end
% Assemble Mass matrix M
M = sparse(mesh.elem(:,[1,1,1,2,2,2,3,3,3]), ...
           mesh.elem(:,[1,2,3,1,2,3,1,2,3]), ...
                  area*[2,1,1,1,2,1,1,1,2]/12, n, n) ;


%--------------------------------------------------------------------------
% Assemble Mass matrix M0 associated with \Omega_0
%--------------------------------------------------------------------------              
% find out elements in \Omega_0: -0.5<x,y<0.5 where f = 1+x+y             
indexe = find(mesh.node(mesh.elem(:,1),1)>-0.5 & mesh.node(mesh.elem(:,1),1)<0.5 ...
            & mesh.node(mesh.elem(:,1),2)>-0.5 & mesh.node(mesh.elem(:,1),2)<0.5 ...
            & mesh.node(mesh.elem(:,2),1)>-0.5 & mesh.node(mesh.elem(:,2),1)<0.5 ...
            & mesh.node(mesh.elem(:,2),2)>-0.5 & mesh.node(mesh.elem(:,2),2)<0.5 ...
            & mesh.node(mesh.elem(:,3),1)>-0.5 & mesh.node(mesh.elem(:,3),1)<0.5 ...
            & mesh.node(mesh.elem(:,3),2)>-0.5 & mesh.node(mesh.elem(:,3),2)<0.5);
% trisurf(mesh.elem, mesh.node(:,1), mesh.node(:,2), zeros(size(mesh.node,1),1));
% view(2), hold on
% node = [mesh.elem(indexe,1);mesh.elem(indexe,2);mesh.elem(indexe,3)];         
% node = unique(node);         
% plot(mesh.node(node,1), mesh.node(node,2),'y.', 'MarkerSize',10)    

M0 = sparse(mesh.elem(indexe,[1,1,1,2,2,2,3,3,3]), ...
           mesh.elem(indexe,[1,2,3,1,2,3,1,2,3]), ...
                  area(indexe)*[2,1,1,1,2,1,1,1,2]/12, n, n);          

%--------------------------------------------------------------------------
% Neumann boundary conditions
%--------------------------------------------------------------------------
 if (~isempty(mesh.Dirichlet))
%     % boundary segment as a vector
     vs(:,:) = mesh.node(mesh.Dirichlet(:,2),:) - ...
               mesh.node(mesh.Dirichlet(:,1),:);
%     % the lenght of boundary segment
     ls = sqrt(sum(vs.^2,2));
     C = C + sparse(mesh.Dirichlet(:,[1 1 2 2]),mesh.Dirichlet(:,[1 2 1 2]),...
                   ls(:)*[2,1,1,2]/6,n,n);
 end
              
%--------------------------------------------------------------------------
% End of function Stiff
%--------------------------------------------------------------------------