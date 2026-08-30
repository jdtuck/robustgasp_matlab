function opt = rgasp_options(n, p, varargin)
%RGASP_OPTIONS  Parse name/value options for RGASP and PPGASP.
%
%   Defaults follow the R package RobustGaSP (version 0.6.x):
%
%     trend            ones(n,1)              constant mean
%     zeroMean         false
%     nugget           0
%     nuggetEst        false
%     rangePar         []                     (fix the range parameters)
%     method           'post_mode'            'post_mode' | 'mmle' | 'mle'
%     a                0.2                    jointly robust prior
%     b                (a+p)/n^(1/p)          jointly robust prior
%     kernelType       'matern_5_2'           'matern_5_2'|'matern_3_2'|'pow_exp'
%     isotropic        false
%     alpha            1.9  (pow_exp only)
%     lowerBound       true
%     maxEval          max(30, 20+5p)
%     initialValues    []                     rows are starting log-parameters
%     numInitialValues 2
%     optimizer        'auto'                 'auto' uses fmincon (or fminunc
%                                             when unconstrained) from the
%                                             Optimization Toolbox and falls
%                                             back to the built-in projected
%                                             L-BFGS; may also be set to
%                                             'fmincon'|'fminunc'|'lbfgs'|
%                                             'neldermead'
%     optimTolFun      1e-8                   solver optimality tolerance
%     optimTolX        1e-10                  solver step tolerance
%     condNumUB        1e16
%     useMex           true                   set false to force the pure
%                                             MATLAB code paths
%     verbose          false

d = struct( ...
    'trend',            ones(n,1), ...
    'zeroMean',         false, ...
    'nugget',           0, ...
    'nuggetEst',        false, ...
    'rangePar',         [], ...
    'method',           'post_mode', ...
    'a',                0.2, ...
    'b',                [], ...
    'kernelType',       'matern_5_2', ...
    'isotropic',        false, ...
    'alpha',            1.9*ones(1,p), ...
    'lowerBound',       true, ...
    'maxEval',          max(30, 20 + 5*p), ...
    'initialValues',    [], ...
    'numInitialValues', 2, ...
    'optimizer',        'auto', ...
    'optimTolFun',      1e-8, ...
    'optimTolX',        1e-10, ...
    'condNumUB',        1e16, ...
    'useMex',           true, ...
    'verbose',          false);

if numel(varargin) == 1 && isstruct(varargin{1})
    userOpts = varargin{1};
    names = fieldnames(userOpts);
    for i = 1:numel(names)
        d = assign(d, names{i}, userOpts.(names{i}));
    end
else
    if mod(numel(varargin), 2) ~= 0
        error('rgasp:options', 'Options must be given as name/value pairs.');
    end
    for i = 1:2:numel(varargin)
        d = assign(d, varargin{i}, varargin{i+1});
    end
end

if ~any(strcmpi(d.method, {'post_mode','mmle','mle'}))
    error('rgasp:options', 'method must be post_mode, mmle or mle.');
end

% Resolve 'auto': prefer the Optimization Toolbox when it is installed.
if strcmpi(d.optimizer, 'auto')
    if d.lowerBound && rgasp_toolbox('fmincon')
        d.optimizer = 'fmincon';
    elseif ~d.lowerBound && rgasp_toolbox('fminunc')
        d.optimizer = 'fminunc';
    elseif rgasp_toolbox('fmincon')
        d.optimizer = 'fmincon';
    else
        d.optimizer = 'lbfgs';
    end
end
if any(strcmpi(d.optimizer, {'fmincon','fminunc'})) && ~rgasp_toolbox(lower(d.optimizer))
    warning('rgasp:optimizer', ...
        ['%s is not available (Optimization Toolbox); falling back to the ' ...
         'built-in projected L-BFGS.'], d.optimizer);
    d.optimizer = 'lbfgs';
end
if d.nuggetEst && d.nugget == 0
    % nugget is estimated: the supplied value is only a starting point
end
opt = d;
end

function d = assign(d, name, value)
fn = fieldnames(d);
idx = find(strcmpi(fn, name), 1);
if isempty(idx)
    error('rgasp:options', 'Unknown option ''%s''.', name);
end
d.(fn{idx}) = value;
end
