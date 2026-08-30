function LB = rgasp_search_lb(design, maxd, isotropic, kernel_type, alpha, ...
                              nugget, cond_num_ub)
%RGASP_SEARCH_LB  Lower bound for log(beta) that keeps R well conditioned.
%
%   LB = RGASP_SEARCH_LB(design, maxd, isotropic, kernel_type, alpha, ...
%                        nugget, cond_num_ub)
%
%   The marginal posterior of the range parameters is flat as gamma -> Inf
%   (beta -> 0), where the correlation matrix degenerates towards a matrix of
%   ones.  RobustGaSP therefore constrains the optimizer with a data-driven
%   lower bound on log(beta): the smallest common "correlation floor" q such
%   that setting
%
%       beta_l = -log(q) / max_ij |x_il - x_jl|
%
%   gives a correlation matrix whose condition number equals cond_num_ub
%   (10^16 by default).  The floor q is parameterized as q = expit(z) and z is
%   found by a golden-section search on [-5, 12], matching search_LB_prob() in
%   the R package.
%
%   maxd is the p-vector of maximum per-dimension distances (or the maximum
%   Euclidean distance when isotropic), supplied by the caller so that no
%   n-by-n distance matrix has to exist.
%
%   Returns a p-vector of lower bounds on log(beta).

if nargin < 7 || isempty(cond_num_ub), cond_num_ub = 1e16; end

p = numel(maxd);
maxd = maxd(:);

% On a large design the conditioning search would dominate the fit, and the
% bound only needs to be roughly right; subsample the rows for the search.
n = size(design, 1);
nsub = min(n, 400);
if nsub < n
    step = max(1, floor(n/nsub));
    sub = design(1:step:n, :);
else
    sub = design;
end

obj = @(z) lb_objective(z, sub, maxd, isotropic, kernel_type, alpha, ...
                        nugget, cond_num_ub, p);
zopt = golden_section(obj, -5, 12, 1e-6, 200);

prob = 1/(1 + exp(-zopt));
LB = zeros(p, 1);
for l = 1:p
    if maxd(l) <= 0
        LB(l) = -Inf;
    else
        LB(l) = log(-log(prob) / maxd(l));
    end
end
end

% ----------------------------------------------------------------------
function v = lb_objective(z, sub, maxd, isotropic, kernel_type, alpha, ...
                          nugget, cond_num_ub, p)
prob = 1/(1 + exp(-z));
b = zeros(p,1);
for l = 1:p
    if maxd(l) <= 0
        b(l) = 1;
    else
        b(l) = -log(prob) / maxd(l);
    end
end
ns = size(sub,1);
R = rgasp_corr_inputs(sub, sub, b, kernel_type, alpha, isotropic) + nugget*eye(ns);
c = cond(R);
if ~isfinite(c), c = 1e300; end
v = (log(c) - log(cond_num_ub))^2;   % log scale: far better conditioned search
end

function xopt = golden_section(f, a, b, tol, maxit)
gr = (sqrt(5) - 1)/2;
c = b - gr*(b - a);
d = a + gr*(b - a);
fc = f(c); fd = f(d);
for it = 1:maxit
    if abs(b - a) < tol, break; end
    if fc < fd
        b = d; d = c; fd = fc;
        c = b - gr*(b - a); fc = f(c);
    else
        a = c; c = d; fc = fd;
        d = a + gr*(b - a); fd = f(d);
    end
end
if fc < fd, xopt = c; else, xopt = d; end
end
