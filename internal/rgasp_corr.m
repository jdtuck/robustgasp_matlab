function R = rgasp_corr(R0, beta, kernel_type, alpha)
%RGASP_CORR  Separable (product) correlation matrix R = R_1 o ... o R_p.
%
%   R = RGASP_CORR(R0, beta, kernel_type, alpha)
%
%   R0          1-by-p cell array of absolute-distance matrices, R0{l}(i,j) =
%               |x_il - x_jl| (or a single Euclidean-distance matrix when the
%               model is isotropic, in which case p = 1).
%   beta        p-vector of inverse range parameters, beta_l = 1/gamma_l.
%   kernel_type 1-by-p cell array of kernel names (one per input dimension).
%   alpha       p-vector of roughness parameters (used by 'pow_exp' only).
%
%   See also RGASP_CORR_1D.

p = numel(R0);
R = ones(size(R0{1}));
for l = 1:p
    R = R .* rgasp_corr_1d(R0{l}, beta(l), kernel_type{l}, alpha(l));
end
end
