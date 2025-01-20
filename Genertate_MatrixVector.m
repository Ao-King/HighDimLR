% Save some matrix, vector 
rt = 2; 
mesh = Mesh(rt);
n = size(mesh.node,1);

[K,M,C,M0] = Stiff(rt); % Assemble coefficient matrix
CM = K+M;               % Coefficient Matrix

%--------------------------------------------------------------------------
% find out elements in \Omega_0: -0.5<x,y<0.5 where f = 1+x+y
%--------------------------------------------------------------------------
indexe = find(mesh.node(mesh.elem(:,1),1)>-0.5 & mesh.node(mesh.elem(:,1),1)<0.5 ...
            & mesh.node(mesh.elem(:,1),2)>-0.5 & mesh.node(mesh.elem(:,1),2)<0.5 ...
            & mesh.node(mesh.elem(:,2),1)>-0.5 & mesh.node(mesh.elem(:,2),1)<0.5 ...
            & mesh.node(mesh.elem(:,2),2)>-0.5 & mesh.node(mesh.elem(:,2),2)<0.5 ...
            & mesh.node(mesh.elem(:,3),1)>-0.5 & mesh.node(mesh.elem(:,3),1)<0.5 ...
            & mesh.node(mesh.elem(:,3),2)>-0.5 & mesh.node(mesh.elem(:,3),2)<0.5);
node0 = [mesh.elem(indexe,1);mesh.elem(indexe,2);mesh.elem(indexe,3)];         
node0 = unique(node0);
node = (1:n)';
pe = sparse(n,1);
pe(node0) = SetSource(mesh.node(node0,1),mesh.node(node0,2));
p2norm = sqrt(pe'*M0*pe);

save infor CM C M0 p2norm node0 pe 

function z = SetSource(x,y)% Source function
         z = 1+x+y;
end