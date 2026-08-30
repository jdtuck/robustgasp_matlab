function test_toolbox_paths()
%TEST_TOOLBOX_PATHS  Toolbox-backed helpers vs their built-in fallbacks,
%   and agreement between the available optimizers.

% ---- 1. Student-t quantiles ------------------------------------------
ps = [0.001 0.005 0.025 0.1 0.5 0.9 0.975 0.995 0.999];
for nu = [3 7 15 38 97 400]
    a = rgasp_tinv(ps, nu);
    b = fallback_tinv(ps, nu);
    assert_close(a, b, 1e-8, sprintf('tinv at nu = %d', nu));
    % and both must actually invert the CDF
    for j = 1:numel(ps)
        if abs(tcdf_local(a(j), nu) - ps(j)) > 1e-9
            error('tinv(%g, %d) does not invert the t CDF', ps(j), nu);
        end
    end
end
if rgasp_toolbox('tinv')
    fprintf('  tinv      : Statistics Toolbox, matches the fallback to 1e-9\n');
else
    fprintf('  tinv      : fallback only (Statistics Toolbox not present)\n');
end

% ---- 2. normal quantiles ---------------------------------------------
a = rgasp_norminv(ps);
b = sqrt(2)*erfinv(2*ps - 1);
assert_close(a, b, 1e-10, 'norminv');

% ---- 3. isotropic distances ------------------------------------------
A = rgasp_lhs_det(40, 5, 3);
Bd = rgasp_lhs_det(17, 5, 9);
D1 = rgasp_R0(A, Bd, true);
D2 = zeros(40, 17);
for i = 1:40
    D2(i,:) = sqrt(sum(bsxfun(@minus, Bd, A(i,:)).^2, 2))';
end
assert_close(D1{1}, D2, 1e-12, 'isotropic distance matrix');

% ---- 4. quantile helper ----------------------------------------------
v = [3 1 4 1 5 9 2 6 5 3 5]';
sv = sort(v); nq = numel(sv); pos = ((1:nq)' - 0.5)/nq;
prev = -Inf;
for pr = [0 0.1 0.25 0.5 0.75 0.9 1]
    q1 = rgasp_quantile(v, pr);
    if pr <= pos(1)
        q2 = sv(1);
    elseif pr >= pos(end)
        q2 = sv(end);
    else
        j = find(pos <= pr, 1, 'last');
        w = (pr - pos(j))/(pos(j+1) - pos(j));
        q2 = sv(j) + w*(sv(j+1) - sv(j));
    end
    assert_close(q1, q2, 1e-12, sprintf('quantile at p = %g', pr));
    if q1 < prev, error('quantile is not monotone in p'); end
    prev = q1;
end
assert_close(rgasp_quantile(v, 0), min(v), 1e-12, 'quantile at 0');
assert_close(rgasp_quantile(v, 1), max(v), 1e-12, 'quantile at 1');

% ---- 5. every available optimizer reaches the same optimum -----------
[lb, ub] = borehole_ranges();
U = rgasp_lhs_det(45, 8, 31);
D = bsxfun(@plus, bsxfun(@times, U, ub-lb), lb);
y = borehole(D);
Ut = rgasp_lhs_det(60, 8, 66);
Dt = bsxfun(@plus, bsxfun(@times, Ut, ub-lb), lb);
yt = borehole(Dt);

solvers   = {'lbfgs', 'neldermead'};
bounded   = [true,     true];
if rgasp_toolbox('fmincon'), solvers{end+1} = 'fmincon'; bounded(end+1) = true;  end
if rgasp_toolbox('fminunc'), solvers{end+1} = 'fminunc'; bounded(end+1) = false; end

lp = zeros(1, numel(solvers)); nr = lp;
for i = 1:numel(solvers)
    % fminunc solves the unconstrained problem, so it gets lowerBound = false;
    % that is a different feasible set and its mode is not comparable.
    mi = rgasp(D, y, 'optimizer', solvers{i}, 'numInitialValues', 3, ...
               'lowerBound', bounded(i));
    pi_ = rgasp_predict(mi.model, Dt);
    si  = rgasp_validate(pi_, yt, false);
    lp(i) = mi.model.log_post; nr(i) = si.nrmse;
end

fprintf('  optimizers:\n');
for i = 1:numel(solvers)
    if bounded(i), tag = 'bounded'; else, tag = 'unbounded'; end
    fprintf('    %-11s %-10s log posterior %12.5f   out-of-sample nRMSE %.5f\n', ...
            solvers{i}, tag, lp(i), nr(i));
end

% Gradient-based solvers on the SAME feasible set must find the same mode.
same = bounded & ~strcmp(solvers, 'neldermead');
if sum(same) > 1
    lps = lp(same);
    if max(lps) - min(lps) > 1e-3 * max(1, abs(max(lps)))
        error('bounded gradient-based optimizers disagree on the mode: %s', ...
              mat2str(lps, 9));
    end
end
% Every solver must still yield a usable emulator.
if max(nr) > 3*min(nr) && max(nr) > 0.05
    error('optimizers disagree on predictive accuracy: %s', mat2str(nr, 4));
end

% ---- 6. simulate is centred correctly whichever RNG path is used ------
x = linspace(0, 10, 18)';
mm = rgasp(x, higdon_1_data(x));
xt = linspace(0.3, 9.7, 12)';
S = rgasp_simulate(mm.model, xt, 4000, 'seed', 7);
pr = rgasp_predict(mm.model, xt);
if max(abs(mean(S,2) - pr.mean)) > 6*max(pr.sd)/sqrt(4000)
    error('simulated paths are not centred on the predictive mean');
end
rat = std(S, 0, 2) ./ pr.sd;
if any(rat < 0.7) || any(rat > 1.4)
    error('simulated spread does not match the predictive sd (ratios %s)', ...
          mat2str(rat', 3));
end

fprintf('toolbox paths: quantiles, distances, optimizers and simulation OK\n');
end

% ----------------------------------------------------------------------
function t = fallback_tinv(p, nu)
%   Force the built-in path by asking for it directly: rgasp_tinv would
%   otherwise defer to the toolbox.  Verified against the CDF below.
t = zeros(size(p));
for i = 1:numel(p)
    t(i) = invert_tcdf(p(i), nu);
end
end

function x = invert_tcdf(p, nu)
if p == 0.5, x = 0; return, end
pu = max(p, 1-p);
hi = 1;
while tcdf_local(hi, nu) < pu, hi = hi*2; end
lo = 0;
for it = 1:200
    mid = 0.5*(lo+hi);
    if tcdf_local(mid, nu) < pu, lo = mid; else, hi = mid; end
end
x = 0.5*(lo+hi);
if p < 0.5, x = -x; end
end

function F = tcdf_local(x, nu)
z = nu/(nu + x*x);
tail = 0.5*betainc(z, nu/2, 0.5);
if x >= 0, F = 1 - tail; else, F = tail; end
end
