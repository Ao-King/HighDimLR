function render_pde(data_directory, output_directory)
% Reproduce the seven PDE confidence-region figures from saved plot data.
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, data_directory = fullfile(root,'data','plotting','pde'); end
if nargin < 2, output_directory = fullfile(root,'outputs','figures'); end
if ~isfolder(output_directory), mkdir(output_directory); end
methods = {'Landweber','Showalter','HBF','ARk','SOAR','Nesterov','FAR'};
for k = 1:numel(methods)
    method = methods{k};
    input_file = fullfile(data_directory,[method '_plot_data.mat']);
    assert(isfile(input_file), 'Missing plot data: %s',input_file);
    identity = load(input_file,'method');
    assert(strcmp(char(string(identity.method)),method), ...
        'The method field does not match %s.',input_file);
    render_pde_method(input_file,fullfile(output_directory,method));
end
end
