% trimesh(mesh.elem(node0), mesh.node(node0,1), mesh.node(node0,2), full(pe(node0)));
% trimesh(mesh.elem, mesh.node(:,1), mesh.node(:,2), full(pe));
node = (1:n)';
node1 = setdiff(node,node0);
[U, S, V] = svds(M0);
S_sqrt  = sqrtm(S);
betanew    = U * S_sqrt * V'*beta1;
betau  = betanew + q1_95;
betau(node1) = 0;
betal  = betanew - q1_95;
betal(node1) = 0;
penew = U * S_sqrt * V'*pe;
figure;
colormap(jet);
hold on;
trimesh(mesh.elem, mesh.node(:,1), mesh.node(:,2), full(betau),'FaceAlpha', 0, 'EdgeColor', 'k');
trimesh(mesh.elem, mesh.node(:,1), mesh.node(:,2), full(betal),'FaceAlpha', 0, 'EdgeColor', 'k');
trisurf(mesh.elem, mesh.node(:,1), mesh.node(:,2), full(penew),'FaceAlpha', 1, 'FaceColor', 'interp', 'EdgeColor', 'interp');
view(3);
axis([-0.5 0.5 -0.5 0.5])
colorbar;
hold off;

% colorbar
% axis([-1 1 -1 1])
% xlabel('x')
% ylabel('y')
% set(gca,'CLim',[0,2])
