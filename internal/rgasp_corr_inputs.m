function R = rgasp_corr_inputs(input1, input2, beta, kernel_type, alpha, isotropic)
%RGASP_CORR_INPUTS  Separable correlation matrix straight from the designs.
%
%   R = RGASP_CORR_INPUTS(input1, input2, beta, kernel_type, alpha, isotropic)
%
%   Equivalent to RGASP_CORR(RGASP_R0(input1, input2, isotropic), ...) but,
%   when the MEX acceleration is available, the p distance matrices are never
%   allocated and the product over dimensions is evaluated in a single pass.
%   That is the difference between p*n^2 and 0 extra doubles, which is what
%   makes n in the thousands practical.

if nargin < 6 || isempty(isotropic), isotropic = false; end
if isempty(input2), input2 = input1; end

if rgasp_have_mex()
    R = rgasp_corr_mex(input1, input2, beta(:), ...
                       rgasp_kernel_code(kernel_type), alpha(:), isotropic);
else
    R0 = rgasp_R0(input1, input2, isotropic);
    R  = rgasp_corr(R0, beta, kernel_type, alpha);
end
end
