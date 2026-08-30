function setup_robustgasp()
%SETUP_ROBUSTGASP  Add the toolbox folders to the MATLAB/Octave path.
%
%   Run this once per session from anywhere:
%       run('/path/to/RobustGaSP-MATLAB/setup_robustgasp.m')

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, 'internal'));
addpath(fullfile(here, 'internal', 'mex'));
addpath(fullfile(here, 'testfun'));
addpath(fullfile(here, 'demos'));
addpath(fullfile(here, 'tests'));
% In Octave the Statistics and Optimization equivalents live in packages that
% must be loaded explicitly; in MATLAB the toolboxes are always on the path.
if exist('OCTAVE_VERSION', 'builtin') == 5
    for pk = {'statistics', 'optim'}
        try %#ok<TRYNC>
            evalc(['pkg load ' pk{1}]);
        end
    end
end

bits = {};
if rgasp_toolbox('fmincon', true), bits{end+1} = 'Optimization Toolbox'; end
if rgasp_toolbox('tinv', true),    bits{end+1} = 'Statistics Toolbox'; end
if rgasp_have_mex(true),           bits{end+1} = 'MEX acceleration'; end
if isempty(bits)
    fprintf(['RobustGaSP-MATLAB added to the path (%s).\n' ...
             'No optional components detected; built-in fallbacks will be used.\n'], here);
else
    fprintf('RobustGaSP-MATLAB added to the path (%s).\nUsing: %s.\n', ...
            here, strjoin(bits, ', '));
end
if ~rgasp_have_mex()
    fprintf('Optional: run rgasp_compile_mex() once for the C++ acceleration.\n');
end
end
