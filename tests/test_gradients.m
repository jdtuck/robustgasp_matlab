function test_gradients()
%TEST_GRADIENTS  Analytic gradient of the negative log marginal posterior
%   checked against central finite differences, over every combination of
%   estimation method, kernel, nugget setting, trend and output dimension.

n = 12;
p = 3;
X0 = rgasp_lhs_det(n, p, 7);
y  = dettepepel_3_data(X0);
Y2 = [y, 2*y + sin(5*X0(:,1)), cos(3*X0(:,2))];   % vector output (k = 3)

kernels  = {'matern_5_2','matern_3_2','pow_exp'};
methods  = {'post_mode','mmle','mle'};
maxerr   = 0;

for ik = 1:numel(kernels)
  for im = 1:numel(methods)
    for nug = [false true]
      for zm = [false true]
        for k = [1 3]
          if k == 1, Yc = y; else, Yc = Y2; end
          m = build_core(X0, Yc, kernels{ik}, methods{im}, nug, zm);
          np = p + double(nug);
          % a few parameter vectors spanning small / medium / large beta
          P = [ -2*ones(1,p),  log(1e-3);
                 0.3*[1 -0.4 0.9], log(2e-2);
                 1.5*ones(1,p), log(1e-1) ];
          for r = 1:size(P,1)
            param = P(r, 1:p)';
            if nug, param = [param; P(r, p+1)]; end %#ok<AGROW>
            [f0, g] = rgasp_objective(param, m);
            gn = zeros(np,1);
            for j = 1:np
              h = 1e-6 * max(1, abs(param(j)));
              e = zeros(np,1); e(j) = h;
              fp = rgasp_objective(param+e, m);
              fm = rgasp_objective(param-e, m);
              gn(j) = (fp - fm)/(2*h);
            end
            den = max(1, max(abs(gn)));
            err = max(abs(g - gn))/den;
            maxerr = max(maxerr, err);
            if ~(err < 2e-5)
              error(['gradient mismatch: kernel=%s method=%s nugget=%d ' ...
                     'zeroMean=%d k=%d row=%d  err=%.3g\n analytic=%s\n numeric =%s'], ...
                 kernels{ik}, methods{im}, nug, zm, k, r, err, ...
                 mat2str(g', 6), mat2str(gn', 6));
            end
          end
        end
      end
    end
  end
end

fprintf('gradients: 108 configurations checked, max relative error %.3g\n', maxerr);
end

% ----------------------------------------------------------------------
function m = build_core(X0, Y, kernel, method, nugget_est, zero_mean)
n = size(X0,1); p = size(X0,2);
m = struct();
m.R0 = rgasp_R0(X0);
if zero_mean
    m.X = zeros(n,0);
else
    m.X = [ones(n,1), X0(:,1)];
end
m.Y = Y;
m.kernel_type = repmat({kernel}, 1, p);
m.alpha = 1.9*ones(1,p);
m.nugget = 0.01;
m.nugget_est = nugget_est;
m.method = method;
m.CL = ((max(X0,[],1) - min(X0,[],1))'/n^(1/p));
m.a = 0.2;
m.b = (0.2 + p)/n^(1/p);
end
