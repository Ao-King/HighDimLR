function mesh = meshstruct(node,elem,Dirichlet)
% MESHSTRUCT1 generates initial mesh data structure.
% 
% USAGE
%    mesh = meshstruct1(node,elem,Neumann)
%
% INPUT 
%       node:  array of nodes (required)
%       elem:  array of element indices (required)
%    Neumann:  Neumann boundary edges (required)
%
% OUTPUT
%    mesh:  current mesh
%

% R.F. Gong 1-28-2019

%-------------------------------------------------------------------------
% Generate mesh data structure
%--------------------------------------------------------------------------
mesh = struct('node',node, 'elem',elem, 'Dirichlet',Dirichlet);

      
%--------------------------------------------------------------------------
% End of function MESHSTRUCT1
%--------------------------------------------------------------------------