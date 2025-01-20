function [mesh] = Mesh(rt)
% Mesh creates the mesh on which the finite element method is implemented
% 
% USAGE
%    [mesh] = Mesh(t)
%
% INPUT
%    rt: the times of the refinement of triangulation
% OUTPUT
%    mesh: mesh structure including information on the nodes, elements,
%    boundary on which Dirichlet and Neumann condtions imposed, and the
%    type of the nodes

%    R.F. Gong 1-28-2019
%--------------------------------------------------------------------------
% Initialize the mesh and refine
%--------------------------------------------------------------------------
[p,e,t] = initmesh('Geom');
p = jigglemesh(p,e,t);
for k = 1:rt
    [p,e,t] = refinemesh('Geom',p,e,t);
    p = jigglemesh(p,e,t);
end
node = p';
elem = t';elem(:,4) = [];

%
%--------------------------------------------------------------------------
% Find out points on the boundary and arrange them anticlockwise
% the domain is \Omega=\{(x,y)\in \Omega\mid x^2+y^2\leq 1\}
%--------------------------------------------------------------------------
ix = find(sum(node.^2,2)>=1-10^(-5));
[temp0,ix0] = sort(node(ix,1));  
ix = ix(ix0);
ix0 = find(node(ix,2)<=0);
ix1 = ix(ix0);
ix2 = ix; 
ix2(ix0)=[];
ix0 = [];
for k=1:length(ix2)                   
    ix0(k) = ix2(length(ix2)-k+1);
end
ix2 = ix0';
ix = [ix1; ix2];         

for k = 1:(length(ix)-1)
    bedge(k,:) = [ix(k) ix(k+1)];   
end
bedge(length(ix),:)=[ix(length(ix)) ix(1)];

%-------------------------------------------------------------------------
% specify Dirichlet/Neumann boundary \Gamma=\partial\Omega
 Dirichlet = bedge;
% Neumann = bedge;

%--------------------------------------------------------------------------
mesh = meshstruct(node,elem,Dirichlet);

%--------------------------------------------------------------------------
% End of function Mesh
%--------------------------------------------------------------------------

