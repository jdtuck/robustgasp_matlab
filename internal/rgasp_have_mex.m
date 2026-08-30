function tf = rgasp_have_mex(refresh)
%RGASP_HAVE_MEX  True when the optional C++ acceleration is available.
%
%   tf = RGASP_HAVE_MEX()        cached check
%   tf = RGASP_HAVE_MEX(true)    re-check (call after RGASP_COMPILE_MEX)
%
%   Build it with RGASP_COMPILE_MEX.  Nothing in the toolbox requires it; the
%   pure-MATLAB paths give identical results.

persistent cached
if nargin >= 1 && refresh
    cached = [];
end
if isempty(cached)
    cached = (exist('rgasp_corr_mex', 'file') == 3) && ...
             (exist('rgasp_gradterms_mex', 'file') == 3);
end
tf = cached;
end
