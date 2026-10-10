if ~exist('source_error_summary', 'var') || ~isstruct(source_error_summary) || ...
        ~all(isfield(source_error_summary, ...
        {'transform', 'mass_matrix', 'original_estimate', 'original_radius'})) || ...
        ~strcmp(source_error_summary.transform, 'full-source-mass-sqrt-v1')
    error('SourceError:LegacyQuantiles', ...
        'Rerun the updated bootstrap driver before plotting full-mass intervals.');
end
if ~isfield(source_error_summary, 'experiment_version') || ...
        ~strcmp(source_error_summary.experiment_version, 'pde-calibrated-bootstrap-v2')
    error('SourceError:LegacyExperiment', ...
        'Rerun the calibrated PDE experiment; a full-mass transform tag alone is insufficient.');
end
if ~isequal(source_error_summary.mass_matrix, M0)
    error('SourceError:MassMatrixMismatch', ...
        'The bootstrap summary and plot must use the same source mass matrix.');
end
if ~isfield(source_error_summary, 'run_status') || ...
        ~strcmp(source_error_summary.run_status, 'complete')
    error('SourceError:IncompleteExperiment', ...
        'Only a completed calibrated experiment can supply confidence-region plots.');
end
node = (1:n)';
node1 = setdiff(node,node0);
M0_sqrt = source_mass_sqrt(M0);
betanew    = M0_sqrt*source_error_summary.original_estimate;
betau  = betanew + source_error_summary.original_radius;
betau(node1) = 0;
betal  = betanew - source_error_summary.original_radius;
betal(node1) = 0;
penew = M0_sqrt*pe;
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

