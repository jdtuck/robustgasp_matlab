function c = rgasp_corr_1d(d, beta, kernel_type, alpha)
%RGASP_CORR_1D  One-dimensional correlation function evaluated on |x_i - x_j|.
%
%   c = RGASP_CORR_1D(d, beta, kernel_type, alpha) returns the correlation
%   matrix (elementwise) for the absolute-distance matrix d, using the
%   inverse-range parameterization beta = 1/gamma advocated in
%
%       Gu, Wang and Berger (2018), "Robust Gaussian Stochastic Process
%       Emulation", Annals of Statistics 46(6A), 3038-3066.
%
%   Supported kernels (Section 2.1 of the paper):
%       'matern_5_2' : (1 + u + u^2/3) exp(-u),  u = sqrt(5) beta d
%       'matern_3_2' : (1 + u) exp(-u),          u = sqrt(3) beta d
%       'pow_exp'    : exp(-(beta d)^alpha),     alpha in (0,2]
%
%   See also RGASP_DLOGCORR_1D, RGASP_CORR.

switch lower(kernel_type)
    case 'matern_5_2'
        u = sqrt(5) * beta * d;
        c = (1 + u + u.^2/3) .* exp(-u);
    case 'matern_3_2'
        u = sqrt(3) * beta * d;
        c = (1 + u) .* exp(-u);
    case 'pow_exp'
        c = exp(-(beta * d).^alpha);
    otherwise
        error('rgasp:kernel', ...
            'Unknown kernel_type ''%s''. Use matern_5_2, matern_3_2 or pow_exp.', ...
            kernel_type);
end
end
