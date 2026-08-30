function pred = gasp_predict(model, testing_input, varargin)
%GASP_PREDICT  Predictive distribution of a fitted (pp)GaSP emulator.
%
%   Shared implementation behind RGASP_PREDICT and PPGASP_PREDICT.
%
%   With the mean and variance integrated out, the predictive distribution at
%   a new input x* is a Student t (eq. (2.10)-(2.13) of Gu, Wang and Berger,
%   2018):
%
%     y(x*) | y^D, gamma_hat  ~  T( yhat(x*), sigma2_hat * c**(x*), n - q )
%
%     yhat(x*) = h(x*) theta_hat + r(x*)' R~^{-1} (y^D - H theta_hat)
%     c**(x*)  = c(x*,x*) + eta - r' R~^{-1} r
%                + (h(x*) - H' R~^{-1} r)' (H' R~^{-1} H)^{-1} (h(x*) - H' R~^{-1} r)
%
%   The eta term is included when 'intervalData' is true (default), i.e. when
%   the interval should cover a *noisy observation* rather than the noise-free
%   mean surface.
%
%   For method 'mle' the plug-in normal predictive distribution is used and
%   the trend-uncertainty term is dropped, matching RobustGaSP.

o = struct('trend', [], 'intervalData', true, 'level', 0.95);
if numel(varargin) == 1 && isstruct(varargin{1})
    f = fieldnames(varargin{1});
    for i = 1:numel(f), o.(f{i}) = varargin{1}.(f{i}); end
else
    for i = 1:2:numel(varargin)
        fn = fieldnames(o);
        j = find(strcmpi(fn, varargin{i}), 1);
        if isempty(j), error('rgasp:options','Unknown option ''%s''.', varargin{i}); end
        o.(fn{j}) = varargin{i+1};
    end
end

if isvector(testing_input) && model.p == 1
    testing_input = testing_input(:);
end
if size(testing_input,2) ~= model.p
    error('rgasp:dim', 'testing_input must have %d columns.', model.p);
end

n  = model.num_obs;
q  = model.q;
k  = model.k;
nt = size(testing_input, 1);

if model.zero_mean
    X_testing = zeros(nt, 0);
else
    X_testing = o.trend;
    if isempty(X_testing)
        if q ~= 1
            error('rgasp:trend', ...
                ['A testing trend matrix must be supplied via ''trend'' when ', ...
                 'the fitted model uses a non-constant mean.']);
        end
        X_testing = ones(nt, 1);
    end
    if size(X_testing,1) ~= nt || size(X_testing,2) ~= q
        error('rgasp:dim', 'trend for prediction must be %d-by-%d.', nt, q);
    end
end

L = model.L;
r = rgasp_corr_inputs(testing_input, model.input, model.beta_hat, ...
                      model.kernel_type, model.alpha, model.isotropic);  % nt-by-n

% r' * R~^{-1}
rt_R_inv = (L' \ (L \ r'))';                     % nt-by-n
rtRinvr  = sum(rt_R_inv .* r, 2);                % nt-by-1

is_mle = strcmpi(model.method, 'mle');
use_trend_var = (q > 0) && ~is_mle;

if use_trend_var
    Rinv_X = L' \ (L \ model.X);                 % n-by-q
    Hdiff  = X_testing - r * Rinv_X;             % nt-by-q
    LX = model.LX;
    tmp = (LX' \ (LX \ Hdiff'));                 % q-by-nt
    diff2 = sum(Hdiff' .* tmp, 1)';              % nt-by-1
else
    diff2 = zeros(nt,1);
end

nug = 0;
if o.intervalData, nug = model.nugget; end

c_star_star = abs(1 + nug - rtRinvr + diff2);

if q > 0
    resid = model.output - model.X * model.theta_hat;
    mu = X_testing * model.theta_hat + rt_R_inv * resid;
else
    mu = rt_R_inv * model.output;
end

sigma2 = model.sigma2_hat(:)';                   % 1-by-k
var_raw = c_star_star * sigma2;                  % nt-by-k

df = n - q;
lo = (1 - o.level)/2;
hi = 1 - lo;
if is_mle
    z_lo = rgasp_norminv(lo);
    z_hi = rgasp_norminv(hi);
    sd = sqrt(var_raw);
else
    z_lo = rgasp_tinv(lo, df);
    z_hi = rgasp_tinv(hi, df);
    if df > 2
        sd = sqrt(var_raw * df/(df-2));
    else
        sd = sqrt(var_raw);
    end
end

pred = struct();
pred.mean    = mu;
pred.sd      = sd;
pred.lower95 = mu + sqrt(var_raw) * z_lo;
pred.upper95 = mu + sqrt(var_raw) * z_hi;
pred.df      = df;
pred.level   = o.level;
pred.c_star_star = c_star_star;
end
