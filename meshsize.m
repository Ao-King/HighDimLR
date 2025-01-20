function [Maxmesh,Minmesh] = meshsize(rt)
% compute meshsize
% element boundary length

mesh = Mesh(rt);
ebl(:,1) = (mesh.node(mesh.elem(:,2),1)-mesh.node(mesh.elem(:,1),1)).^2 + ...
      (mesh.node(mesh.elem(:,2),2)-mesh.node(mesh.elem(:,1),2)).^2;
ebl(:,2) = (mesh.node(mesh.elem(:,3),1)-mesh.node(mesh.elem(:,1),1)).^2 + ...
      (mesh.node(mesh.elem(:,3),2)-mesh.node(mesh.elem(:,1),2)).^2;
ebl(:,3) = (mesh.node(mesh.elem(:,3),1)-mesh.node(mesh.elem(:,2),1)).^2 + ...
      (mesh.node(mesh.elem(:,3),2)-mesh.node(mesh.elem(:,2),2)).^2;
ebl = sqrt(ebl);
Maxmesh = max(max(ebl)); %meshsize
Minmesh = min(min(ebl)); %meshsize