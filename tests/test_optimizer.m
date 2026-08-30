function test_optimizer()
%TEST_OPTIMIZER  Bound-constrained L-BFGS on problems with known solutions.

% --- 1. quadratic, unconstrained ---------------------------------------
A = [3 1; 1 2]; bvec = [1; -2];
f1 = @(x) deal(0.5*x'*A*x - bvec'*x, A*x - bvec);
[x, fv] = rgasp_lbfgs(f1, [5; -5], [], [], 500, struct());
assert_close(x, A\bvec, 1e-6, 'quadratic minimizer');

% --- 2. Rosenbrock, unconstrained --------------------------------------
[x, fv] = rgasp_lbfgs(@rosen, [-1.2; 1], [], [], 5000, struct());
assert_close(x, [1;1], 1e-4, 'Rosenbrock minimizer');
if fv > 1e-8, error('Rosenbrock objective not minimized (f = %g)', fv); end

% --- 3. active lower bound ---------------------------------------------
lb = [2; -Inf];
[x, ~] = rgasp_lbfgs(f1, [5; 5], lb, [], 500, struct());
if abs(x(1) - 2) > 1e-8
    error('lower bound not active: x1 = %g', x(1));
end
% with x1 fixed at 2 the optimum in x2 solves 2*x2 + 1*2 = -2
assert_close(x(2), -2, 1e-6, 'constrained second coordinate');

% --- 4. Nelder-Mead fallback -------------------------------------------
[x, ~] = rgasp_neldermead(@rosen, [-1.2; 1], [], [], 20000);
assert_close(x, [1;1], 1e-2, 'Nelder-Mead on Rosenbrock');

fprintf('optimizer: quadratic, Rosenbrock, active bound and simplex fallback OK\n');
end

function [f, g] = rosen(x)
f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
g = [-400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1)); 200*(x(2)-x(1)^2)];
end
