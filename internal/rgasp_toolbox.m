function tf = rgasp_toolbox(name, refresh)
%RGASP_TOOLBOX  Is a toolbox function available in this installation?
%
%   tf = RGASP_TOOLBOX('fmincon')
%   tf = RGASP_TOOLBOX('lhsdesign', true)     % ignore the cache
%
%   The toolbox is written to USE the Statistics and Optimization Toolboxes
%   when they are present: fmincon drives the marginal posterior mode search,
%   tinv/norminv give the predictive quantiles, lhsdesign generates designs,
%   mvnrnd/chi2rnd draw predictive sample paths, and pdist2 forms isotropic
%   distance matrices.  Pure-MATLAB equivalents are retained only so the code
%   still runs where a toolbox is missing.  RGASP_TOOLBOX_REPORT prints which
%   path is actually in use.

persistent names flags
if nargin < 2, refresh = false; end

% Global override (see RGASP_USE_TOOLBOXES): lets you check that a result does
% not depend on a toolbox by forcing every built-in fallback.
if ~rgasp_use_toolboxes()
    tf = false;
    return
end
if refresh
    names = {}; flags = [];
end
if isempty(names)
    names = {}; flags = logical([]);
end

i = find(strcmp(names, name), 1);
if ~isempty(i)
    tf = flags(i);
    return
end

tf = exist(name, 'file') > 0 || exist(name, 'builtin') > 0;
names{end+1} = name; %#ok<AGROW>
flags(end+1) = tf;   %#ok<AGROW>
end
