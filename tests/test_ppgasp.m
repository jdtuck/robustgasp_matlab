function test_ppgasp()
%TEST_PPGASP  Parallel partial GaSP for vector-valued output.

% ---- 1. k = 1 must reproduce rgasp exactly ----------------------------
x = linspace(0, 10, 15)';
y = higdon_1_data(x);
m1 = rgasp(x, y);
m2 = ppgasp(x, y);
assert_close(m2.beta_hat, m1.model.beta_hat, 1e-12, 'ppgasp with k=1 range parameters');
assert_close(m2.sigma2_hat, m1.model.sigma2_hat, 1e-12, 'ppgasp with k=1 variance');
p1 = rgasp_predict(m1.model, x); p2 = ppgasp_predict(m2, x);
assert_close(p2.mean, p1.mean, 1e-12, 'ppgasp with k=1 predictive mean');
assert_close(p2.sd,   p1.sd,   1e-12, 'ppgasp with k=1 predictive sd');

% ---- 2. environmental spill model: 200 outputs per run ----------------
lb = [7,   0.02, 0.01, 30.01];
ub = [13,  0.12, 3.00, 30.295];
nd = 40; nt = 25;
Ud = rgasp_lhs_det(nd, 4, 2);
Ut = rgasp_lhs_det(nt, 4, 77);
Dd = bsxfun(@plus, bsxfun(@times, Ud, ub-lb), lb);
Dt = bsxfun(@plus, bsxfun(@times, Ut, ub-lb), lb);
Yd = environ_4_data(Dd);
Yt = environ_4_data(Dt);

mp = ppgasp(Dd, Yd);
if mp.k ~= size(Yd,2)
    error('ppgasp did not record the number of output coordinates');
end
if numel(mp.sigma2_hat) ~= mp.k
    error('ppgasp should estimate one variance per output coordinate');
end
pp = ppgasp_predict(mp, Dt);
if any(size(pp.mean) ~= [nt, mp.k])
    error('ppgasp_predict returned the wrong shape');
end
sp = rgasp_validate(pp, Yt, false);
if sp.nrmse > 0.05
    error('ppgasp emulator too inaccurate: normalized RMSE = %.4g', sp.nrmse);
end
if sp.coverage < 0.85
    error('ppgasp under-covers: %.3f', sp.coverage);
end

% ---- 3. interpolation at the design points ---------------------------
pd = ppgasp_predict(mp, Dd);
assert_close(pd.mean, Yd, 1e-5, 'ppgasp interpolation');

% ---- 4. shared range parameters, separate means/variances -------------
%   Fixing gamma at the PP GaSP estimate and fitting each output column on
%   its own must reproduce the PP GaSP column-wise mean and variance.
cols = [1, 37, 120, mp.k];
for c = cols
    mc = rgasp(Dd, Yd(:,c), 'rangePar', mp.range_hat, 'nugget', mp.nugget);
    assert_close(mc.model.sigma2_hat, mp.sigma2_hat(c), 1e-8, ...
                 sprintf('column %d variance', c));
    pc = rgasp_predict(mc.model, Dt);
    assert_close(pc.mean, pp.mean(:,c), 1e-8, sprintf('column %d mean', c));
    assert_close(pc.sd,   pp.sd(:,c),   1e-8, sprintf('column %d sd', c));
end

% ---- 5. the joint likelihood is not the sum of independent fits --------
%   (sanity check that the range parameters really are shared)
ma = rgasp(Dd, Yd(:,1));
if abs(ma.model.beta_hat(1) - mp.beta_hat(1)) < 1e-12
    error('PP GaSP range estimate coincides exactly with a single-column fit');
end

fprintf('ppgasp: k=1 equivalence, field emulation accuracy, interpolation\n');
fprintf('        and shared-range/separate-variance structure OK\n');
end
