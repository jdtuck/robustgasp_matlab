function test_rgasp_fit_predict()
%TEST_RGASP_FIT_PREDICT  Fitting and prediction for scalar-output models.

% ---- 1. interpolation: the emulator reproduces the training data -------
x = linspace(0, 10, 15)';
y = higdon_1_data(x);
m = rgasp(x, y);
pr = rgasp_predict(m, x);
assert_close(pr.mean, y, 1e-7, 'emulator must interpolate the design points');
if max(pr.sd) > 1e-5
    error('predictive sd at design points should vanish (max = %.3g)', max(pr.sd));
end

% ---- 2. predictive formulas reproduced by direct linear algebra --------
xt = linspace(0, 10, 37)';
pr = rgasp_predict(m, xt);
[mu_ref, var_ref] = brute_force_predict(m, xt);
assert_close(pr.mean, mu_ref, 1e-10, 'predictive mean');
% (the variance is compared at a looser tolerance: the reference uses an
%  explicit inverse of an ill-conditioned R, the implementation uses Cholesky)
assert_close(pr.lower95, mu_ref + sqrt(var_ref)*rgasp_tinv(0.025, m.num_obs-m.q), ...
             1e-5, 'lower 95% bound');
assert_close(pr.upper95, mu_ref + sqrt(var_ref)*rgasp_tinv(0.975, m.num_obs-m.q), ...
             1e-5, 'upper 95% bound');
df = m.num_obs - m.q;
assert_close(pr.sd, sqrt(var_ref*df/(df-2)), 1e-5, 'predictive sd');

% ---- 3. accuracy on a dense grid --------------------------------------
xt = linspace(0, 10, 500)';
pr = rgasp_predict(m, xt);
st = rgasp_validate(pr, higdon_1_data(xt), false);
if st.nrmse > 0.05
    error('1-d emulator too inaccurate: normalized RMSE = %.4g', st.nrmse);
end
if st.coverage < 0.90
    error('1-d emulator under-covers: %.3f', st.coverage);
end

% ---- 4. eight-dimensional borehole function ---------------------------
[lb, ub] = borehole_ranges();
n = 80; nt = 300;
Xd = rgasp_lhs_det(n, 8, 11);
Xt = rgasp_lhs_det(nt, 8, 999);
Dd = bsxfun(@plus, bsxfun(@times, Xd, ub-lb), lb);
Dt = bsxfun(@plus, bsxfun(@times, Xt, ub-lb), lb);
yd = borehole(Dd); ytr = borehole(Dt);
mb = rgasp(Dd, yd);
pb = rgasp_predict(mb, Dt);
sb = rgasp_validate(pb, ytr, false);
if sb.nrmse > 0.02
    error('borehole emulator too inaccurate: normalized RMSE = %.4g', sb.nrmse);
end
if sb.coverage < 0.90
    error('borehole emulator under-covers: %.3f', sb.coverage);
end

% ---- 5. nugget estimation on noisy data -------------------------------
rng_u = rgasp_lcg(60, 3);
xr = 10*rgasp_lcg(60, 21);
noise = 0.1 * sqrt(-2*log(rng_u)) .* cos(2*pi*rgasp_lcg(60, 5));   % Box-Muller
yn = higdon_1_data(xr) + noise;
mn = rgasp(xr, yn, 'nuggetEst', true, 'numInitialValues', 3);
if mn.nugget <= 0 || ~isfinite(mn.nugget)
    error('nugget estimation failed: eta = %g', mn.nugget);
end
noise_var_hat = mn.nugget * mn.sigma2_hat;
if noise_var_hat < 0.002 || noise_var_hat > 0.08
    error('estimated noise variance %.4g is far from the true 0.01', noise_var_hat);
end

% ---- 6. mmle and mle also run, and interpolate -------------------------
for meth = {'mmle','mle'}
    mm = rgasp(x, y, 'method', meth{1});
    pm = rgasp_predict(mm, x);
    assert_close(pm.mean, y, 1e-6, sprintf('%s interpolation', meth{1}));
end

% ---- 7. zero mean, non-constant trend, isotropic, other kernels --------
m0 = rgasp(x, y, 'zeroMean', true);
p0 = rgasp_predict(m0, x);
assert_close(p0.mean, y, 1e-6, 'zero-mean interpolation');

H  = [ones(15,1), x];
mt = rgasp(x, y, 'trend', H);
pt = rgasp_predict(mt, x, 'trend', H);
assert_close(pt.mean, y, 1e-6, 'linear-trend interpolation');

X2 = rgasp_lhs_det(30, 2, 4);
y2 = limetal_2_data(X2);
mi = rgasp(X2, y2, 'isotropic', true);
if numel(mi.beta_hat) ~= 1
    error('isotropic model should have a single range parameter');
end
for kt = {'matern_3_2','pow_exp'}
    mk = rgasp(X2, y2, 'kernelType', kt{1});
    pk = rgasp_predict(mk, X2);
    assert_close(pk.mean, y2, 1e-5, sprintf('%s interpolation', kt{1}));
end

% ---- 8. simulate ------------------------------------------------------
S = rgasp_simulate(m, xt(1:20), 200, 'seed', 1);
if any(size(S) ~= [20 200])
    error('rgasp_simulate returned the wrong shape');
end
pr20 = rgasp_predict(m, xt(1:20));
if max(abs(mean(S,2) - pr20.mean)) > 5*max(pr20.sd)/sqrt(200) + 1e-8
    error('simulated paths are not centred on the predictive mean');
end

fprintf('rgasp: interpolation, predictive formulas, borehole accuracy,\n');
fprintf('       nugget estimation, alternative kernels and simulation OK\n');
end

% ----------------------------------------------------------------------
function [mu, v] = brute_force_predict(m, xt)
%   Independent re-derivation of eq. (2.10)-(2.13) with plain inverses.
R0 = rgasp_R0(m.input);
R  = rgasp_corr(R0, m.beta_hat, m.kernel_type, m.alpha) + m.nugget*eye(m.num_obs);
Ri = inv(R);
X  = m.X;
theta = (X'*Ri*X) \ (X'*Ri*m.output);
sigma2 = (m.output - X*theta)' * Ri * (m.output - X*theta) / (m.num_obs - m.q);
r0 = rgasp_R0(xt, m.input);
r  = rgasp_corr(r0, m.beta_hat, m.kernel_type, m.alpha);
Xt = ones(size(xt,1), 1);
mu = Xt*theta + r*Ri*(m.output - X*theta);
v  = zeros(size(xt,1),1);
for i = 1:size(xt,1)
    ri = r(i,:)';
    hd = Xt(i,:)' - X'*Ri*ri;
    cst = 1 + m.nugget - ri'*Ri*ri + hd'*((X'*Ri*X)\hd);
    v(i) = sigma2 * abs(cst);
end
end
