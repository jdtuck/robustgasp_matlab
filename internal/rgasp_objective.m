function [f, g, aux] = rgasp_objective(param, m)
%RGASP_OBJECTIVE  Negative log marginal posterior / likelihood and its gradient.
%
%   [f, g, aux] = RGASP_OBJECTIVE(param, m)
%
%   param   Optimization variable.  param = log(beta) when the nugget is not
%           estimated, and param = [log(beta); log(eta)] when it is.  The log
%           parameterization xi_l = log(beta_l) = -log(gamma_l) is the
%           "robust" parameterization recommended in Gu, Wang and Berger
%           (2018): the marginal posterior mode in this parameterization is
%           bounded away from both 0 and infinity under mild conditions.
%
%   m       Model struct with fields
%             input       n-by-p design (always present)
%             R0          1-by-p cell of distance matrices, or {} when the
%                         MEX path is in use
%             isotropic   logical
%             X           n-by-q trend matrix ([] when zero_mean)
%             Y           n-by-k output matrix (k = 1 for rgasp)
%             kernel_type 1-by-p cellstr
%             alpha       1-by-p roughness parameters
%             nugget      fixed nugget (used when nugget_est is false)
%             nugget_est  logical
%             method      'post_mode' | 'mmle' | 'mle'
%             CL, a, b    jointly robust prior quantities (post_mode only)
%
%   Returns the NEGATIVE of
%
%     post_mode / mmle  (marginal likelihood, mean and variance integrated out
%                        with pi(theta,sigma^2) ~ 1/sigma^2; eq. (2.5)-(2.7)):
%
%       l = -k*sum(log(diag(L))) - k*sum(log(diag(LX)))
%           - ((n-q)/2) * sum_i log(S2_i)                    [+ log pi_JR]
%
%     mle               (profile likelihood, eq. (2.9)):
%
%       l = -k*sum(log(diag(L))) - (n/2) * sum_i log(S2_i)
%
%   where R~ = R(beta) + eta*I = L L', S2_i = y_i' Q y_i and
%   Q = R~^{-1} - R~^{-1} X (X' R~^{-1} X)^{-1} X' R~^{-1}.
%
%   For k > 1 this is exactly the PP GaSP (parallel partial GaSP) marginal
%   likelihood: the range parameters are shared across the k output
%   coordinates while the mean vector and variance of each coordinate are
%   integrated out separately.
%
%   The gradient is taken with respect to param (i.e. the chain rule factor
%   beta_l, and eta, is already applied).
%
%   Analytic gradient (with W_l = Rdot_l * Q, Rdot_l = dR/dbeta_l):
%
%       dl/dbeta_l = -(k/2) tr(W_l)
%                    + ((n-q)/2) * sum_i (Qy_i' Rdot_l Qy_i) / S2_i
%
%   which is the derivative of eq. (2.7); the trace term for 'mle' uses
%   tr(R~^{-1} Rdot_l) instead of tr(W_l).
%
%   IMPLEMENTATION NOTE.  The trend correction enters only through the low
%   rank factor B with B*B' = R~^{-1}X (X'R~^{-1}X)^{-1} X'R~^{-1}, so neither
%   that product nor Q is ever formed as an n-by-n matrix:
%
%       tr(Rdot_l Q) = sum_ij Ri_ij Rdot_ij - sum_c B(:,c)' Rdot_l B(:,c)
%       QY           = Ri*Y - B*(B'*Y)
%
%   Together with recovering the nugget-free R from R~ on the diagonal, this
%   keeps the working set at three n-by-n arrays (R~, its Cholesky factor and
%   its inverse) regardless of p.

n  = size(m.Y, 1);
k  = size(m.Y, 2);
p  = numel(m.kernel_type);
nugget_est = m.nugget_est;
use_mex = rgasp_have_mex() && isfield(m, 'input') && ~isempty(m.input);
if isfield(m, 'force_matlab') && m.force_matlab, use_mex = false; end
if ~use_mex && (~isfield(m, 'R0') || isempty(m.R0))
    m.R0 = rgasp_R0(m.input, m.input, m.isotropic);
end

if nugget_est
    beta = exp(param(1:p));
    nu   = exp(param(p+1));
else
    beta = exp(param(:));
    nu   = m.nugget;
end
beta = beta(:);

% ---- R~ = R(beta) + eta I, built in place ----------------------------
if use_mex
    Rt = rgasp_corr_mex(m.input, m.input, beta, ...
                        rgasp_kernel_code(m.kernel_type), m.alpha(:), m.isotropic);
else
    Rt = rgasp_corr(m.R0, beta, m.kernel_type, m.alpha);
end
idx = 1:(n+1):n*n;
Rt(idx) = Rt(idx) + nu;                      % now it really is R~

[U, chol_fail] = chol(Rt);                   % R~ = U'U
if chol_fail ~= 0
    f = 1e100; g = zeros(numel(param), 1); aux = struct(); return
end

% Guard against numerically singular correlation matrices: the marginal
% likelihood is flat as beta -> 0, where R degenerates to a matrix of ones,
% and the optimizer must be pushed back rather than allowed to invert it.
dU = diag(U);
if min(dU) <= 0 || max(dU)/min(dU) > 1e9
    f = 1e100; g = zeros(numel(param), 1); aux = struct(); return
end

is_mle = strcmpi(m.method, 'mle');
logdetL = sum(log(dU));

Ri = U \ (U' \ eye(n));                      % R~^{-1}

if isempty(m.X)
    q = 0;
    B = zeros(n, 0);
    logdetLX = 0;
    QY = Ri * m.Y;
else
    q = size(m.X, 2);
    Rinv_X = Ri * m.X;
    W = m.X' * Rinv_X;
    W = (W + W')/2;
    [UW, cf2] = chol(W);
    if cf2 ~= 0
        f = 1e100; g = zeros(numel(param),1); aux = struct(); return
    end
    logdetLX = sum(log(diag(UW)));
    B  = Rinv_X / UW;                        % B*B' = Rinv_X W^{-1} Rinv_X'
    QY = Ri * m.Y - B * (B' * m.Y);
end

S2 = sum(m.Y .* QY, 1)';
if any(S2 <= 0) || any(~isfinite(S2))
    f = 1e100; g = zeros(numel(param),1); aux = struct(); return
end

if is_mle
    coef = n;
    ll   = -k*logdetL - (coef/2)*sum(log(S2));
else
    coef = n - q;
    ll   = -k*logdetL - k*logdetLX - (coef/2)*sum(log(S2));
end

% ---- jointly robust prior (Gu 2019); 'ref_approx' in RobustGaSP -------
add_prior = strcmpi(m.method, 'post_mode');
t = NaN;
if add_prior
    t  = sum(m.CL(:) .* beta) + nu;
    ll = ll + m.a*log(t) - m.b*t;
end

f = -ll;

aux = struct();
if nargout < 2
    if nargout >= 3, aux = pack_aux(); end
    return
end

% ---- gradient ---------------------------------------------------------
% For 'mle' the trace term is tr(R~^{-1} Rdot) and the trend correction is
% dropped; otherwise it is tr(Q Rdot) = tr(R~^{-1} Rdot) - sum_c B_c' Rdot B_c.
use_trend_term = ~is_mle && q > 0;

g = zeros(numel(param), 1);

if use_mex
    if use_trend_term, Buse = B; else, Buse = zeros(n,0); end
    [tr_l, quad_l] = rgasp_gradterms_mex(m.input, Rt, nu, Ri, Buse, QY, S2, ...
        beta, rgasp_kernel_code(m.kernel_type), m.alpha(:), m.isotropic, ...
        use_trend_term);
    g(1:p) = -( -(k/2)*tr_l + (coef/2)*quad_l ) .* beta;
else
    if use_trend_term
        Tmat = Ri - B*B';
    else
        Tmat = Ri;
    end
    Rnf = Rt; Rnf(idx) = Rnf(idx) - nu;        % nugget-free correlation
    for l = 1:p
        dlc  = rgasp_dlogcorr_1d(m.R0{l}, beta(l), m.kernel_type{l}, m.alpha(l));
        dR   = Rnf .* dlc;                     % dR/dbeta_l
        tr_t = sum(sum(Tmat .* dR));
        quad = sum(sum(QY .* (dR * QY), 1)' ./ S2);
        g(l) = -( -(k/2)*tr_t + (coef/2)*quad ) * beta(l);
    end
end

if nugget_est
    if use_trend_term
        tr_t = trace(Ri) - sum(sum(B.^2));     % tr(Q), Rdot = I
    else
        tr_t = trace(Ri);
    end
    quad = sum(sum(QY .* QY, 1)' ./ S2);
    g(p+1) = -( -(k/2)*tr_t + (coef/2)*quad ) * nu;
end

if add_prior
    dprior_dbeta = m.a*m.CL(:)/t - m.b*m.CL(:);
    g(1:p) = g(1:p) - dprior_dbeta .* beta;
    if nugget_est
        g(p+1) = g(p+1) - (m.a/t - m.b) * nu;
    end
end

if nargout >= 3
    aux = pack_aux();
end

    function a_ = pack_aux()
        a_.U    = U;
        a_.L    = U';
        a_.S2   = S2;
        a_.q    = q;
        a_.beta = beta;
        a_.nu   = nu;
        if isempty(m.X)
            a_.LX = [];
            a_.theta_hat = zeros(0, k);
        else
            a_.LX = UW';
            a_.theta_hat = UW \ (UW' \ (m.X' * (Ri * m.Y)));
        end
        if is_mle
            a_.sigma2_hat = S2 / n;
        else
            a_.sigma2_hat = S2 / (n - q);
        end
    end
end
