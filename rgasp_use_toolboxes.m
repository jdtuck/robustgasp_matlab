function state = rgasp_use_toolboxes(tf)
%RGASP_USE_TOOLBOXES  Enable or disable use of the optional toolboxes.
%
%   rgasp_use_toolboxes(false)   force every built-in fallback
%   rgasp_use_toolboxes(true)    use the toolboxes again (the default)
%   state = rgasp_use_toolboxes()
%
%   Useful for checking that a result does not depend on the Statistics or
%   Optimization Toolbox being installed -- the test suite runs itself both
%   ways.  It does not affect the MEX acceleration; see RGASP_HAVE_MEX.
%
%   See also RGASP_TOOLBOX_REPORT, RGASP_COMPILE_MEX.

persistent enabled
if isempty(enabled), enabled = true; end
if nargin >= 1 && ~isempty(tf)
    enabled = logical(tf);
    rgasp_toolbox('tinv', true);      % drop the availability cache
end
state = enabled;
end
