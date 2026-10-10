function config = paper_table6_config(method,overrides)
%PAPER_TABLE6_CONFIG Parameters used for the seven PDE experiments in Table 6.
% Example: paper_table6_config('SOAR',struct('N',1,'B',8)).
if nargin < 2, overrides = struct(); end
config = pde_experiment_config(method,overrides);
config.execution_backend = 'cpu-double-batched';
config.far_backend = 'far-spectral';
end
