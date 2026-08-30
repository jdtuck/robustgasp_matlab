function test_kernels()
%TEST_KERNELS  Correlation functions and their analytic derivatives.

d = [0 0.3 1.7; 0.3 0 2.0; 1.7 2.0 0];

% ---- closed forms ----------------------------------------------------
beta = 0.7;
u = sqrt(5)*beta*d;
assert_close(rgasp_corr_1d(d, beta, 'matern_5_2', 1.9), ...
             (1+u+u.^2/3).*exp(-u), 1e-14, 'matern 5/2 value');
u3 = sqrt(3)*beta*d;
assert_close(rgasp_corr_1d(d, beta, 'matern_3_2', 1.9), ...
             (1+u3).*exp(-u3), 1e-14, 'matern 3/2 value');
assert_close(rgasp_corr_1d(d, beta, 'pow_exp', 1.9), ...
             exp(-(beta*d).^1.9), 1e-14, 'pow exp value');

% correlation matrices must have unit diagonal and be symmetric
for kt = {'matern_5_2','matern_3_2','pow_exp'}
    C = rgasp_corr_1d(d, beta, kt{1}, 1.9);
    assert_close(diag(C), ones(3,1), 1e-14, 'unit diagonal');
    assert_close(C, C', 1e-14, 'symmetry');
    if max(abs(C(:))) > 1 + 1e-12
        error('correlation exceeds 1 for %s', kt{1});
    end
end

% ---- d log c / d beta against central differences ---------------------
h = 1e-6;
for kt = {'matern_5_2','matern_3_2','pow_exp'}
    for beta = [0.05 0.5 3 20]
        analytic = rgasp_dlogcorr_1d(d, beta, kt{1}, 1.9);
        cp = rgasp_corr_1d(d, beta+h, kt{1}, 1.9);
        cm = rgasp_corr_1d(d, beta-h, kt{1}, 1.9);
        numeric = (log(cp) - log(cm)) / (2*h);
        mask = isfinite(numeric) & (cp > 1e-250);
        assert_close(analytic(mask), numeric(mask), 1e-5, ...
            sprintf('dlogc/dbeta for %s at beta=%g', kt{1}, beta));
    end
end

% ---- separable product ----------------------------------------------
X = [0 0; 1 0.5; 0.3 0.9];
R0 = rgasp_R0(X);
b  = [1.3; 0.4];
R  = rgasp_corr(R0, b, {'matern_5_2','matern_3_2'}, [1.9 1.9]);
Rref = rgasp_corr_1d(R0{1}, b(1), 'matern_5_2', 1.9) .* ...
       rgasp_corr_1d(R0{2}, b(2), 'matern_3_2', 1.9);
assert_close(R, Rref, 1e-14, 'separable product');

% ---- isotropic distances --------------------------------------------
R0iso = rgasp_R0(X, X, true);
assert_close(R0iso{1}(1,2), sqrt(1^2 + 0.5^2), 1e-14, 'isotropic distance');

fprintf('kernels: values, symmetry, derivatives and product structure OK\n');
end
