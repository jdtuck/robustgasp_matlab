function [x, fval, info] = rgasp_lbfgs(fun, x0, lb, ub, maxeval, opts)
%RGASP_LBFGS  Self-contained bound-constrained limited-memory BFGS.
%
%   [x, fval, info] = RGASP_LBFGS(fun, x0, lb, ub, maxeval, opts)
%
%   Minimizes fun (which returns [f, g]) subject to lb <= x <= ub using a
%   projected L-BFGS iteration with an Armijo backtracking line search along
%   the projected arc.  This mirrors the low-storage quasi-Newton method
%   (Nocedal, 1980) used by RobustGaSP, and keeps the toolbox dependency-free
%   so it also runs under GNU Octave.
%
%   opts fields (all optional):
%       m        memory size                       (default 10)
%       tol_g    projected-gradient tolerance      (default 1e-6)
%       tol_f    relative function tolerance       (default 1e-12)
%       verbose  print iterations                  (default false)
%
%   info has fields: iterations, funcCount, message, converged.

if nargin < 6 || isempty(opts), opts = struct(); end
mem     = getopt(opts, 'm', 10);
tol_g   = getopt(opts, 'tol_g', 1e-6);
tol_f   = getopt(opts, 'tol_f', 1e-12);
verbose = getopt(opts, 'verbose', false);

x0 = x0(:);
nx = numel(x0);
if isempty(lb), lb = -inf(nx,1); end
if isempty(ub), ub =  inf(nx,1); end
lb = lb(:); ub = ub(:);

x = min(max(x0, lb), ub);
[f, g] = fun(x);
nfev = 1;

S = zeros(nx, 0);
Yv = zeros(nx, 0);
rho = zeros(1, 0);

info.converged = false;
info.message   = 'maximum number of evaluations reached';
iter = 0;

while nfev < maxeval
    iter = iter + 1;

    % projected gradient (for the stopping rule)
    pg = g;
    pg(x <= lb + 0 & g > 0) = 0;
    pg(x >= ub - 0 & g < 0) = 0;
    if norm(pg, inf) < tol_g
        info.converged = true;
        info.message = 'projected gradient below tolerance';
        break
    end

    % ---- two-loop recursion ------------------------------------------
    d = -lbfgs_direction(g, S, Yv, rho);

    % zero out directions that push further into an active bound
    atlb = (x <= lb) & (d < 0);
    atub = (x >= ub) & (d > 0);
    d(atlb | atub) = 0;

    if all(d == 0) || (g' * d) >= 0
        d = -pg;                       % fall back to steepest descent
        if all(d == 0)
            info.converged = true;
            info.message = 'no feasible descent direction';
            break
        end
    end

    % ---- projected Armijo backtracking -------------------------------
    t = 1;
    c1 = 1e-4;
    accepted = false;
    for ls = 1:30
        xn = min(max(x + t*d, lb), ub);
        s  = xn - x;
        if norm(s, inf) < 1e-16
            break
        end
        [fn, gn] = fun(xn);
        nfev = nfev + 1;
        if fn <= f + c1 * (g' * s) || fn < f - abs(f)*1e-14
            accepted = true;
            break
        end
        t = t / 2;
        if nfev >= maxeval, break; end
    end

    if ~accepted
        info.message = 'line search failed';
        break
    end

    df = f - fn;
    yv = gn - g;
    x = xn; f = fn; g = gn;

    % ---- curvature-safe memory update --------------------------------
    sty = s' * yv;
    if sty > 1e-12 * norm(s) * norm(yv) && norm(s) > 0
        S   = [S, s];      Yv = [Yv, yv];      rho = [rho, 1/sty]; %#ok<AGROW>
        if size(S,2) > mem
            S(:,1) = []; Yv(:,1) = []; rho(1) = [];
        end
    end

    if verbose
        fprintf('  iter %3d  f = %.10g  |pg| = %.3g  step = %.3g\n', ...
                iter, f, norm(pg,inf), t);
    end

    if df >= 0 && df <= tol_f * max(1, abs(f))
        info.converged = true;
        info.message = 'function value converged';
        break
    end
end

fval = f;
info.iterations = iter;
info.funcCount  = nfev;
end

% ----------------------------------------------------------------------
function d = lbfgs_direction(g, S, Y, rho)
kmem = size(S, 2);
q = g;
a = zeros(kmem, 1);
for i = kmem:-1:1
    a(i) = rho(i) * (S(:,i)' * q);
    q = q - a(i) * Y(:,i);
end
if kmem > 0
    gamma = (S(:,end)' * Y(:,end)) / (Y(:,end)' * Y(:,end));
else
    gamma = 1;
end
if ~isfinite(gamma) || gamma <= 0, gamma = 1; end
r = gamma * q;
for i = 1:kmem
    b = rho(i) * (Y(:,i)' * r);
    r = r + S(:,i) * (a(i) - b);
end
d = r;
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
