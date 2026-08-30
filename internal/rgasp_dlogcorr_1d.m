function dlc = rgasp_dlogcorr_1d(d, beta, kernel_type, alpha)
%RGASP_DLOGCORR_1D  d log c_l(d) / d beta_l for the 1-d correlation functions.
%
%   Because the separable correlation is a Hadamard product,
%       R = R_1 o R_2 o ... o R_p,
%   the derivative with respect to beta_l is simply
%       dR/dbeta_l = R .* dlogc_l,
%   which avoids ever dividing by a (possibly underflowed) correlation entry.
%
%   Closed forms used here:
%       matern_5_2 : -sqrt(5) d * u (1+u) / (3 + 3u + u^2),  u = sqrt(5) beta d
%       matern_3_2 : -sqrt(3) d * u / (1+u),                 u = sqrt(3) beta d
%       pow_exp    : -alpha beta^(alpha-1) d^alpha
%
%   The Matern 5/2 expression is algebraically identical to the one used by
%   the R package RobustGaSP, but written in a form that is stable for large u.
%
%   See also RGASP_CORR_1D.

switch lower(kernel_type)
    case 'matern_5_2'
        u   = sqrt(5) * beta * d;
        dlc = -sqrt(5) * d .* (u .* (1 + u)) ./ (3 + 3*u + u.^2);
    case 'matern_3_2'
        u   = sqrt(3) * beta * d;
        dlc = -sqrt(3) * d .* u ./ (1 + u);
    case 'pow_exp'
        if beta <= 0
            error('rgasp:kernel', 'beta must be positive for pow_exp.');
        end
        dlc = -alpha * beta^(alpha - 1) * d.^alpha;
    otherwise
        error('rgasp:kernel', 'Unknown kernel_type ''%s''.', kernel_type);
end
end
