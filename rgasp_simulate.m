function S = rgasp_simulate(model, testing_input, num_sample, varargin)
%RGASP_SIMULATE  Draw sample paths from the predictive distribution.
%
%   S = RGASP_SIMULATE(model, testing_input, num_sample)
%   S = RGASP_SIMULATE(model, testing_input, num_sample, 'Name', Value, ...)
%
%   Returns an nt-by-num_sample matrix whose columns are joint draws from the
%   posterior predictive process
%
%       y(x*_1:nt) | y^D, gamma_hat ~ MVT( yhat, sigma2_hat * Sigma*, n - q )
%
%   with
%       Sigma* = C(x*,x*) + eta I - r R~^{-1} r'
%                + (H* - r R~^{-1} H) (H' R~^{-1} H)^{-1} (H* - r R~^{-1} H)'
%
%   OPTIONS: 'trend', 'intervalData' (default true), 'seed'.
%
%   Only scalar-output (rgasp) models are supported.
%
%   See also RGASP_PREDICT.

if ~strcmp(model.type, 'rgasp')
    error('rgasp:type', 'rgasp_simulate supports scalar-output models only.');
end
if nargin < 3 || isempty(num_sample), num_sample = 1; end

o = struct('trend', [], 'intervalData', true, 'seed', []);
for i = 1:2:numel(varargin)
    fn = fieldnames(o);
    j = find(strcmpi(fn, varargin{i}), 1);
    if isempty(j), error('rgasp:options','Unknown option ''%s''.', varargin{i}); end
    o.(fn{j}) = varargin{i+1};
end
if ~isempty(o.seed)
    try
        rng(o.seed);                                   % MATLAB
    catch
        rand('state', o.seed); randn('state', o.seed); %#ok<RAND>  % Octave
    end
end

if isvector(testing_input) && model.p == 1, testing_input = testing_input(:); end
nt = size(testing_input,1);
n  = model.num_obs;
q  = model.q;
L  = model.L;

if model.zero_mean
    X_testing = zeros(nt,0);
else
    X_testing = o.trend;
    if isempty(X_testing)
        if q ~= 1
            error('rgasp:trend', 'Supply ''trend'' for a non-constant mean model.');
        end
        X_testing = ones(nt,1);
    end
end

r = rgasp_corr_inputs(testing_input, model.input, model.beta_hat, ...
                      model.kernel_type, model.alpha, model.isotropic);   % nt-by-n
C = rgasp_corr_inputs(testing_input, testing_input, model.beta_hat, ...
                      model.kernel_type, model.alpha, model.isotropic);   % nt-by-nt

rt_R_inv = (L' \ (L \ r'))';
Sigma = C - rt_R_inv * r';
if o.intervalData
    Sigma = Sigma + model.nugget * eye(nt);
end

if q > 0 && ~strcmpi(model.method,'mle')
    Rinv_X = L' \ (L \ model.X);
    Hdiff  = X_testing - r * Rinv_X;
    LX = model.LX;
    Sigma = Sigma + Hdiff * (LX' \ (LX \ Hdiff'));
end

if q > 0
    mu = X_testing * model.theta_hat + rt_R_inv * (model.output - model.X*model.theta_hat);
else
    mu = rt_R_inv * model.output;
end

Sigma = model.sigma2_hat * (Sigma + Sigma')/2;
% jitter until positive definite
jit = 0;
[U, flag] = chol(Sigma);
while flag ~= 0
    if jit == 0, jit = 1e-12 * max(diag(Sigma)); else, jit = jit * 10; end
    Sigma = Sigma + jit*eye(nt);
    [U, flag] = chol(Sigma);
    if jit > max(diag(Sigma))
        error('rgasp:simulate', 'Predictive covariance is not positive definite.');
    end
end

% Gaussian part: MVNRND (Statistics Toolbox) when available.
if rgasp_toolbox('mvnrnd')
    G = mvnrnd(zeros(1, nt), Sigma, num_sample)';   % nt-by-num_sample
else
    G = U' * randn(nt, num_sample);
end

if strcmpi(model.method, 'mle')
    S = repmat(mu, 1, num_sample) + G;
else
    % Multivariate t: scale the Gaussian draw by sqrt(df/W), W ~ chi2(df).
    df = n - q;
    if rgasp_toolbox('chi2rnd')
        W = chi2rnd(df, 1, num_sample);
    else
        W = sum(randn(df, num_sample).^2, 1);
    end
    S = repmat(mu, 1, num_sample) + G .* repmat(sqrt(df ./ W), nt, 1);
end
end
