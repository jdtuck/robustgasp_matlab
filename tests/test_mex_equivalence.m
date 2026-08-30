function test_mex_equivalence()
%TEST_MEX_EQUIVALENCE  The C++ path must agree with the pure MATLAB path.
%
%   Skipped with a message when the MEX files have not been compiled.

if ~rgasp_have_mex()
    fprintf('mex: not compiled -- skipped (run rgasp_compile_mex)\n');
    return
end

% ---- 1. correlation matrices ------------------------------------------
n1 = 60; n2 = 25; p = 4;
A = rgasp_lhs_det(n1, p, 12);
Bd = rgasp_lhs_det(n2, p, 34);
beta = [0.4; 1.7; 0.05; 3.2];
alpha = [1.9 1.5 2.0 1.2];

for kt = {'matern_5_2','matern_3_2','pow_exp'}
    K = repmat(kt, 1, p);
    Rm = rgasp_corr(rgasp_R0(A, Bd), beta, K, alpha);
    Rc = rgasp_corr_mex(A, Bd, beta, rgasp_kernel_code(K), alpha(:), false);
    assert_close(Rm, Rc, 1e-13, sprintf('%s cross-correlation', kt{1}));

    Rm = rgasp_corr(rgasp_R0(A), beta, K, alpha);
    Rc = rgasp_corr_mex(A, A, beta, rgasp_kernel_code(K), alpha(:), false);
    assert_close(Rm, Rc, 1e-13, sprintf('%s self-correlation', kt{1}));
end

% mixed kernels across dimensions
K = {'matern_5_2','pow_exp','matern_3_2','matern_5_2'};
Rm = rgasp_corr(rgasp_R0(A, Bd), beta, K, alpha);
Rc = rgasp_corr_mex(A, Bd, beta, rgasp_kernel_code(K), alpha(:), false);
assert_close(Rm, Rc, 1e-13, 'mixed-kernel correlation');

% isotropic
Rm = rgasp_corr(rgasp_R0(A, Bd, true), beta(1), {'matern_5_2'}, alpha(1));
Rc = rgasp_corr_mex(A, Bd, beta(1), 3, alpha(1), true);
assert_close(Rm, Rc, 1e-13, 'isotropic correlation');

% ---- 2. objective and gradient ----------------------------------------
n = 40; p = 3;
X0 = rgasp_lhs_det(n, p, 7);
y  = dettepepel_3_data(X0);
Y3 = [y, 2*y + sin(5*X0(:,1)), cos(3*X0(:,2))];

maxerr_f = 0; maxerr_g = 0;
for kt = {'matern_5_2','matern_3_2','pow_exp'}
  for meth = {'post_mode','mmle','mle'}
    for nug = [false true]
      for zm = [false true]
        for iso = [false true]
          for k = [1 3]
            if k == 1, Yc = y; else, Yc = Y3; end
            core = build_core(X0, Yc, kt{1}, meth{1}, nug, zm, iso);
            pe = numel(core.kernel_type);
            for xi0 = [-1.2 0.4 1.8]
                param = xi0*ones(pe,1);
                if nug, param = [param; log(0.02)]; end %#ok<AGROW>
                core.force_matlab = true;
                [f1, g1] = rgasp_objective(param, core);
                core.force_matlab = false;
                core.R0 = {};
                [f2, g2] = rgasp_objective(param, core);
                core.R0 = rgasp_R0(X0, X0, iso);
                ef = abs(f1-f2)/max(1,abs(f1));
                eg = max(abs(g1-g2))/max(1,max(abs(g1)));
                maxerr_f = max(maxerr_f, ef);
                maxerr_g = max(maxerr_g, eg);
                if ~(ef < 1e-11 && eg < 1e-9)
                    error(['MEX/MATLAB mismatch: %s %s nugget=%d zeroMean=%d ' ...
                           'iso=%d k=%d xi=%g -> df=%.3g dg=%.3g'], ...
                          kt{1}, meth{1}, nug, zm, iso, k, xi0, ef, eg);
                end
            end
          end
        end
      end
    end
  end
end

% ---- 3. end-to-end fits must agree ------------------------------------
[lb, ub] = borehole_ranges();
U = rgasp_lhs_det(50, 8, 21);
D = bsxfun(@plus, bsxfun(@times, U, ub-lb), lb);
yb = borehole(D);
Ut = rgasp_lhs_det(40, 8, 88);
Dt = bsxfun(@plus, bsxfun(@times, Ut, ub-lb), lb);

cases = { {'numInitialValues',3}, ...
          {'nuggetEst',true,'numInitialValues',3}, ...
          {'method','mmle'}, ...
          {'zeroMean',true}, ...
          {'isotropic',true}, ...
          {'kernelType','pow_exp'} };
% The objective and the gradient agree to machine precision (checked above),
% but a full fit is an optimizer trajectory: on the nearly flat ridge left by
% inert inputs a 1e-14 difference can move where the solver stops.  The fitted
% MODEL is what has to agree, so the tolerances here are on the attained log
% posterior and on the predictions, not on the last digit of beta.
sy = std(yb, 1);
for c = 1:numel(cases)
    m1 = rgasp(D, yb, cases{c}{:}, 'useMex', false);
    m2 = rgasp(D, yb, cases{c}{:}, 'useMex', true);
    assert_close(m1.log_post, m2.log_post, 1e-7, sprintf('case %d log posterior', c));
    assert_close(m1.beta_hat, m2.beta_hat, 1e-3, sprintf('case %d beta', c));
    assert_close(m1.sigma2_hat, m2.sigma2_hat, 1e-3, sprintf('case %d sigma2', c));
    p1 = rgasp_predict(m1, Dt); p2 = rgasp_predict(m2, Dt);
    if max(abs(p1.mean - p2.mean))/sy > 1e-5
        error('case %d: predictive means differ by %.2g sd(y)', ...
              c, max(abs(p1.mean - p2.mean))/sy);
    end
    assert_close(p1.sd, p2.sd, 1e-3, sprintf('case %d predictive sd', c));
end

% ppgasp end to end
Ue = rgasp_lhs_det(30, 4, 2);
lbe = [7, 0.02, 0.01, 30.010]; ube = [13, 0.12, 3.00, 30.295];
De = bsxfun(@plus, bsxfun(@times, Ue, ube-lbe), lbe);
Ye = environ_4_data(De);
q1 = ppgasp(De, Ye, 'useMex', false);
q2 = ppgasp(De, Ye, 'useMex', true);
assert_close(q1.log_post, q2.log_post, 1e-7, 'ppgasp log posterior');
assert_close(q1.beta_hat, q2.beta_hat, 1e-3, 'ppgasp beta');
assert_close(q1.sigma2_hat, q2.sigma2_hat, 1e-3, 'ppgasp sigma2');

fprintf('mex: correlation, objective, gradient and end-to-end fits agree\n');
fprintf('     with the pure MATLAB paths (max |df| %.2g, max |dg| %.2g)\n', ...
        maxerr_f, maxerr_g);
end

% ----------------------------------------------------------------------
function m = build_core(X0, Y, kernel, method, nugget_est, zero_mean, iso)
n = size(X0,1); p = size(X0,2);
pe = p; if iso, pe = 1; end
m = struct();
m.input = X0;
m.isotropic = iso;
m.R0 = rgasp_R0(X0, X0, iso);
if zero_mean
    m.X = zeros(n,0);
else
    m.X = [ones(n,1), X0(:,1)];
end
m.Y = Y;
m.kernel_type = repmat({kernel}, 1, pe);
m.alpha = 1.9*ones(1,pe);
m.nugget = 0.01;
m.nugget_est = nugget_est;
m.method = method;
m.CL = ones(pe,1)*0.3;
m.a = 0.2;
m.b = (0.2 + pe)/n^(1/pe);
m.force_matlab = false;
end
