function T = bench_rgasp(ns, p, k)
%BENCH_RGASP  End-to-end fit timing, pure MATLAB vs the C++ acceleration.
%
%   T = bench_rgasp([250 500 1000 2000], 8)
%   T = bench_rgasp([250 500], 6, 50)      % PP GaSP with k = 50 outputs
%
%   Reports the wall-clock time of a complete rgasp/ppgasp fit (identical
%   settings, identical results) with 'useMex' false and true, and the peak
%   number of n-by-n double arrays each path has to hold.
%
%   Also reports the cost of one objective+gradient evaluation, which is what
%   the optimizer actually repeats.

if nargin < 1 || isempty(ns), ns = [250 500 1000]; end
if nargin < 2 || isempty(p),  p  = 8; end
if nargin < 3 || isempty(k),  k  = 1; end

have = rgasp_have_mex();
if ~have
    fprintf(['MEX not compiled -- run rgasp_compile_mex() to see the ' ...
             'comparison.\n']);
end

fprintf('\n p = %d, k = %d, Matern 5/2, constant mean, 2 restarts\n', p, k);
fprintf('%6s | %10s %10s %8s | %10s %10s %8s | %6s %6s\n', ...
        'n', 'fit .m', 'fit mex', 'x', 'obj+g .m', 'obj+g mex', 'x', ...
        'RAM .m', 'RAMmex');
fprintf('%6s | %10s %10s %8s | %10s %10s %8s | %6s %6s\n', ...
        '', '(s)', '(s)', '', '(s)', '(s)', '', '(MB)', '(MB)');

T = zeros(numel(ns), 5);
for ii = 1:numel(ns)
    n = ns(ii);
    X0 = rgasp_lhs_det(n, p, 5);
    y0 = friedman_5_data(X0(:,1:min(5,p)));
    if k == 1
        Y = y0;
    else
        Y = bsxfun(@times, y0, 1:k) + repmat(sin((1:k)/3), n, 1);
    end

    t = tic;
    if k == 1, m1 = rgasp(X0, Y, 'useMex', false);
    else,      m1 = ppgasp(X0, Y, 'useMex', false); end
    t1 = toc(t);

    t2 = NaN; m2 = m1;
    if have
        t = tic;
        if k == 1, m2 = rgasp(X0, Y, 'useMex', true);
        else,      m2 = ppgasp(X0, Y, 'useMex', true); end
        t2 = toc(t);
        % The objective and gradient agree to machine precision, but on a
        % flat ridge the optimizer can stop at slightly different points, so
        % compare the attained log posterior rather than the raw parameters.
        rel = abs(m1.log_post - m2.log_post)/max(1, abs(m1.log_post));
        if rel > 1e-4
            warning('bench:mismatch', ...
                'log posteriors differ by %.2g at n = %d', rel, n);
        end
    end

    % single objective + gradient evaluation
    core = make_core(X0, Y, p, n);
    param = zeros(p,1);
    core.force_matlab = true;  core.R0 = rgasp_R0(X0);
    reps = max(1, round(3e8/n^3));
    t = tic; for r = 1:reps, [~,~] = rgasp_objective(param, core); end
    t3 = toc(t)/reps;
    t4 = NaN;
    if have
        core.force_matlab = false; core.R0 = {};
        t = tic; for r = 1:reps, [~,~] = rgasp_objective(param, core); end
        t4 = toc(t)/reps;
    end

    mb = n*n*8/1024^2;
    ram_m   = (3 + p) * mb;      % R~, chol, inverse, plus p distance matrices
    ram_mex = 3 * mb;

    fprintf('%6d | %10.3f %10.3f %7.1fx | %10.4f %10.4f %7.1fx | %6.0f %6.0f\n', ...
            n, t1, t2, t1/t2, t3, t4, t3/t4, ram_m, ram_mex);
    T(ii,:) = [n t1 t2 t3 t4];
end
fprintf('\n');
end

function core = make_core(X0, Y, p, n)
core = struct();
core.input = X0;
core.isotropic = false;
core.R0 = {};
core.X = ones(n,1);
core.Y = Y;
core.kernel_type = repmat({'matern_5_2'}, 1, p);
core.alpha = 1.9*ones(1,p);
core.nugget = 1e-6;
core.nugget_est = false;
core.method = 'post_mode';
core.CL = ((max(X0)-min(X0))'/n^(1/p));
core.a = 0.2;
core.b = (0.2+p)/n^(1/p);
core.force_matlab = false;
end
