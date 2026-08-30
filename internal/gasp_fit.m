function model = gasp_fit(design, response, varargin)
%GASP_FIT  Shared fitting engine for RGASP (k = 1) and PPGASP (k > 1).
%
%   Not intended to be called directly -- use RGASP or PPGASP.
%
%   Implements marginal posterior mode estimation of the range parameters
%   (and, optionally, the nugget-variance ratio) for the Gaussian stochastic
%   process emulator of
%
%     Gu, M., Wang, X. and Berger, J.O. (2018). Robust Gaussian stochastic
%     process emulation. The Annals of Statistics 46(6A), 3038-3066.
%
%   The estimation is performed in the log inverse-range parameterization
%   xi = log(beta) = -log(gamma), which is the parameterization the paper
%   shows to give robust estimates, with the jointly robust prior
%   (Gu, 2019) on (beta, eta) supplying the tail penalty that removes the
%   degenerate modes of the likelihood.

% ---------------------------------------------------------------- inputs
design = double(design);
if isvector(design), design = design(:); end
Y = double(response);
if isvector(Y) && size(Y,1) ~= size(design,1)
    Y = Y(:);
end
if size(Y,1) ~= size(design,1)
    error('rgasp:dim', 'design and response must have the same number of rows.');
end

n = size(design, 1);
p = size(design, 2);
k = size(Y, 2);

opt = rgasp_options(n, p, varargin{:});

% trend
if opt.zeroMean
    X = zeros(n, 0);
    q = 0;
else
    X = opt.trend;
    if isempty(X), X = ones(n,1); end
    if size(X,1) ~= n
        error('rgasp:dim', 'trend must have %d rows.', n);
    end
    q = size(X, 2);
end
if n - q <= 0
    error('rgasp:dim', 'Need n > q (number of runs greater than trend rank).');
end

% kernels
kernel_type = opt.kernelType;
if ischar(kernel_type), kernel_type = repmat({kernel_type}, 1, p); end
if numel(kernel_type) == 1 && p > 1, kernel_type = repmat(kernel_type, 1, p); end
alpha = opt.alpha(:)';
if numel(alpha) == 1, alpha = repmat(alpha, 1, p); end

isotropic = opt.isotropic;
p_eff = 1;
if ~isotropic, p_eff = p; end
if isotropic
    kernel_type = kernel_type(1);
    alpha = alpha(1);
end

% The p per-dimension distance matrices cost p*n^2 doubles.  The MEX path
% computes the distances on the fly, so they are only built when it is absent
% (and always for the lower-bound search, which needs max|x_il - x_jl|).
use_mex = rgasp_have_mex() && opt.useMex;
if use_mex
    R0 = {};
else
    R0 = rgasp_R0(design, design, isotropic);
end

% ------------------------------------------- jointly robust prior scaling
CL = zeros(p_eff, 1);
if isotropic
    CL(1) = max_dist(design, true, 1) / n;
else
    for l = 1:p_eff
        CL(l) = (max(design(:,l)) - min(design(:,l))) / n^(1/p_eff);
    end
end
CL(CL <= 0) = eps;

a = opt.a;
if isempty(opt.b) || isnan(opt.b)
    b = (a + p_eff) / n^(1/p_eff);
else
    b = opt.b;
end

% ------------------------------------------------------- model container
m = struct();
m.R0          = R0;
m.force_matlab = ~use_mex;
m.input       = design;
m.isotropic   = isotropic;
m.X           = X;
m.Y           = Y;
m.kernel_type = kernel_type;
m.alpha       = alpha;
m.nugget      = opt.nugget;
m.nugget_est  = opt.nuggetEst;
m.method      = opt.method;
m.CL          = CL;
m.a           = a;
m.b           = b;

% --------------------------------------------------------- lower bounds
maxd = zeros(p_eff,1);
for l = 1:p_eff
    maxd(l) = max_dist(design, isotropic, l);
end
LB_cond = rgasp_search_lb(design, maxd, isotropic, kernel_type, alpha, ...
                          opt.nugget, opt.condNumUB);
if opt.lowerBound
    LB = LB_cond;
else
    LB = -inf(p_eff, 1);
end
if opt.nuggetEst
    LB = [LB; -Inf];
end
UB = inf(numel(LB), 1);

% --------------------------------------------- fixed range parameters?
if ~isempty(opt.rangePar) && all(isfinite(opt.rangePar))
    gamma_fixed = opt.rangePar(:);
    if numel(gamma_fixed) ~= p_eff
        error('rgasp:dim', 'rangePar must have %d elements.', p_eff);
    end
    param = log(1 ./ gamma_fixed);
    if opt.nuggetEst
        param = [param; log(max(opt.nugget, 1e-8))];
    end
    [f, ~, aux] = rgasp_objective(param, m);
    model = finalize(param, f, aux);
    model.optim_info = struct('message', 'range parameters supplied by user', ...
                              'iterations', 0, 'funcCount', 1, 'converged', true);
    return
end

% ------------------------------------------------------ starting values
ini = opt.initialValues;
if isempty(ini)
    nini = max(1, opt.numInitialValues);
    beta_ini = zeros(nini, p_eff);
    eta_ini  = zeros(nini, 1);

    base = LB_cond;
    bad  = ~isfinite(base);
    base(bad) = log(1 ./ CL(bad));

    beta_ini(1,:) = 50 * exp(base)';
    eta_ini(1)    = 1e-4;
    if nini > 1
        beta_ini(2,:) = ((a + p_eff) ./ (p_eff * CL * b) / 2)';
        eta_ini(2)    = 2e-4;
    end
    if nini > 2
        % Deterministic pseudo-random restarts (reproducible in MATLAB and
        % Octave); the R package uses set.seed(i)/runif here.
        u = rgasp_lcg((nini-2)*(p_eff+1), 17);
        pos = 0;
        for i = 3:nini
            beta_ini(i,:) = 1e3 * u(pos + (1:p_eff))' ./ CL';
            eta_ini(i)    = 1e-3 * u(pos + p_eff + 1);
            pos = pos + p_eff + 1;
        end
    end
    if opt.nuggetEst
        ini = [log(beta_ini), log(eta_ini)];
    else
        ini = log(beta_ini);
    end
end

% ------------------------------------------------------------ optimize
objfun = @(prm) rgasp_objective(prm, m);
best_f = Inf; best_param = []; best_info = [];

% Near-singular correlation matrices are an expected, and handled, part of
% the search (the likelihood is flat as beta -> 0).  Silence the linear
% algebra warnings for the duration of the optimization only.
ws = rgasp_mute_warnings();

for i = 1:size(ini,1)
    x0 = ini(i,:)';
    x0 = min(max(x0, LB), UB);
    x0(~isfinite(x0)) = 0;
    switch lower(opt.optimizer)
        case 'fmincon'
            so = rgasp_optimoptions('fmincon', opt);
            [xi, fi, ef, out] = fmincon(objfun, x0, [], [], [], [], ...
                                        LB, UB, [], so);
            infoi = solver_info('fmincon', ef, out);
        case 'fminunc'
            so = rgasp_optimoptions('fminunc', opt);
            [xi, fi, ef, out] = fminunc(objfun, x0, so);
            xi = min(max(xi, LB), UB);
            infoi = solver_info('fminunc', ef, out);
        case 'lbfgs'
            [xi, fi, infoi] = rgasp_lbfgs(objfun, x0, LB, UB, opt.maxEval, ...
                                          struct('verbose', opt.verbose));
        case 'neldermead'
            [xi, fi, infoi] = rgasp_neldermead(objfun, x0, LB, UB, opt.maxEval);
        otherwise
            error('rgasp:optimizer', 'Unknown optimizer ''%s''.', opt.optimizer);
    end
    if fi < best_f
        best_f = fi; best_param = xi; best_info = infoi;
    end
end

warning(ws);

if isempty(best_param)
    error('rgasp:optim', 'Optimization failed for every starting value.');
end

[f, ~, aux] = rgasp_objective(best_param, m);
model = finalize(best_param, f, aux);
model.optim_info = best_info;

% ----------------------------------------------------------------------
    function mo = finalize(param, fval, aux)
        mo = struct();
        mo.input       = design;
        mo.output      = Y;
        mo.X           = X;
        mo.zero_mean   = opt.zeroMean;
        mo.num_obs     = n;
        mo.p           = p;
        mo.p_eff       = p_eff;
        mo.q           = q;
        mo.k           = k;
        mo.method      = opt.method;
        mo.kernel_type = kernel_type;
        mo.alpha       = alpha;
        mo.isotropic   = isotropic;
        mo.CL          = CL;
        mo.a           = a;
        mo.b           = b;
        mo.LB          = LB;
        mo.nugget_est  = opt.nuggetEst;
        mo.beta_hat    = aux.beta;
        mo.range_hat   = 1 ./ aux.beta;
        mo.nugget      = aux.nu;
        mo.theta_hat   = aux.theta_hat;
        mo.sigma2_hat  = aux.sigma2_hat;
        mo.L           = aux.L;
        mo.LX          = aux.LX;
        mo.log_post    = -fval;
        mo.param       = param;
        if k == 1
            mo.type = 'rgasp';
        else
            mo.type = 'ppgasp';
        end
    end
end

% ----------------------------------------------------------------------
function d = max_dist(design, isotropic, l)
%   max |x_il - x_jl| (or the largest Euclidean distance when isotropic),
%   computed without materializing an n-by-n matrix.
if isotropic
    n = size(design,1);
    d = 0;
    for i = 1:n
        df = bsxfun(@minus, design, design(i,:));
        d = max(d, sqrt(max(sum(df.^2, 2))));
    end
else
    d = max(design(:,l)) - min(design(:,l));
end
end

% ----------------------------------------------------------------------
function info = solver_info(name, exitflag, out)
%   Normalize the Optimization Toolbox output struct.  MATLAB reports
%   iterations/funcCount, Octave's optim package reports niter/nobjf.
info = struct('message', name, 'iterations', NaN, 'funcCount', NaN, ...
              'converged', exitflag > 0, 'exitflag', exitflag);
if isstruct(out)
    if isfield(out, 'iterations'), info.iterations = out.iterations;
    elseif isfield(out, 'niter'),  info.iterations = out.niter; end
    if isfield(out, 'funcCount'),  info.funcCount = out.funcCount;
    elseif isfield(out, 'nobjf'),  info.funcCount = out.nobjf; end
    if isfield(out, 'message') && ~isempty(out.message)
        info.message = out.message;
    end
end
end
