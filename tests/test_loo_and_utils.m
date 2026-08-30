function test_loo_and_utils()
%TEST_LOO_AND_UTILS  Leave-one-out CV, inert-input detection, validation.

% ---- 1. leave-one-out on a smooth 1-d function ------------------------
x = linspace(0, 10, 20)';
y = higdon_1_data(x);
m = rgasp(x, y);
loo = rgasp_loo(m.model);

if numel(loo.mean) ~= 20
    error('rgasp_loo returned the wrong number of predictions');
end
% LOO predictions must NOT equal the data (that would mean the point was
% still in the training set) but must still be accurate.
if max(abs(loo.resid)) < 1e-10
    error('leave-one-out prediction reproduced the held-out point exactly');
end
if loo.rmse > 0.15*std(y,1)
    error('leave-one-out RMSE too large: %.4g', loo.rmse);
end
if loo.coverage < 0.8
    error('leave-one-out coverage too low: %.3f', loo.coverage);
end

% ---- 2. a directly recomputed LOO point must match ---------------------
i = 7;
idx = [1:i-1, i+1:20];
mi = rgasp(x(idx), y(idx), 'rangePar', m.model.range_hat, 'nugget', m.model.nugget);
pi_ = rgasp_predict(mi.model, x(i));
assert_close(pi_.mean, loo.mean(i), 1e-8, 'LOO mean vs refitted model');

% ---- 3. inert inputs: two pure noise-free dummy coordinates ------------
n = 50;
U = rgasp_lhs_det(n, 4, 33);
yy = limetal_2_data(U(:,1:2));       % coordinates 3 and 4 do nothing
mi2 = rgasp(U, yy, 'numInitialValues', 3);
P = rgasp_inert_inputs(mi2.model, 0.1, false);
if numel(P) ~= 4
    error('rgasp_inert_inputs returned the wrong length');
end
assert_close(sum(P), 4, 1e-10, 'normalized inverse ranges must average to 1');
if ~(P(3) < 0.1 && P(4) < 0.1)
    error('dummy inputs were not detected as inert: P = %s', mat2str(P', 4));
end
if ~(P(1) > 0.1 && P(2) > 0.1)
    error('active inputs were wrongly flagged as inert: P = %s', mat2str(P', 4));
end

% ---- 4. borehole: r and Tl are known to be nearly inert ---------------
[lb, ub] = borehole_ranges();
Ub = rgasp_lhs_det(60, 8, 5);
Db = bsxfun(@plus, bsxfun(@times, Ub, ub-lb), lb);
mb = rgasp(Db, borehole(Db), 'numInitialValues', 3);
Pb = rgasp_inert_inputs(mb.model, 0.1, false);
if Pb(1) < max(Pb([2 3 5]))
    error('borehole: rw should dominate the normalized inverse ranges');
end

% ---- 5. validation metrics ------------------------------------------
xt = linspace(0, 10, 200)';
pr = rgasp_predict(m.model, xt);
st = rgasp_validate(pr, higdon_1_data(xt), false);
inside = (higdon_1_data(xt) >= pr.lower95) & (higdon_1_data(xt) <= pr.upper95);
assert_close(st.coverage, mean(inside), 1e-12, 'coverage definition');
assert_close(st.length, mean(pr.upper95 - pr.lower95), 1e-12, 'length definition');
assert_close(st.rmse, sqrt(mean((higdon_1_data(xt)-pr.mean).^2)), 1e-12, 'rmse');

fprintf('utilities: leave-one-out CV, inert-input detection and validation OK\n');
end
