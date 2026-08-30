function rgasp_toolbox_report()
%RGASP_TOOLBOX_REPORT  Which toolbox functions RobustGaSP-MATLAB will use.
%
%   RGASP_TOOLBOX_REPORT() lists, for every optional dependency, whether the
%   toolbox function is available and what the toolbox falls back to when it
%   is not.

rows = {
 'fmincon',    'Optimization',  'default optimizer for the posterior mode', 'rgasp_lbfgs (built-in projected L-BFGS)'
 'fminunc',    'Optimization',  'unconstrained fits (lowerBound = false)',  'rgasp_lbfgs'
 'optimoptions','Optimization', 'solver option objects',                    'optimset'
 'tinv',       'Statistics',    'Student-t predictive quantiles',           'betaincinv-based quantile'
 'norminv',    'Statistics',    'normal quantiles (method = mle)',          'erfinv-based quantile'
 'lhsdesign',  'Statistics',    'Latin hypercube designs in the demos',     'rgasp_lhs_det (deterministic maximin LHS)'
 'mvnrnd',     'Statistics',    'joint draws in rgasp_simulate',            'Cholesky + randn'
 'chi2rnd',    'Statistics',    'the t scaling in rgasp_simulate',          'sum of squared randn'
 'pdist2',     'Statistics',    'isotropic distance matrices',              'explicit bsxfun loop'
 'quantile',   'Statistics',    'box plot summaries in the demos',          'linear interpolation'
 'boxplot',    'Statistics',    'demo_robustness figure',                   'hand-drawn box plot'
};

fprintf('\nRobustGaSP-MATLAB optional dependencies\n');
fprintf('%s\n', repmat('-', 1, 92));
fprintf('%-14s %-13s %-42s %s\n', 'function', 'toolbox', 'used for', 'status');
fprintf('%s\n', repmat('-', 1, 92));
nmiss = 0;
for i = 1:size(rows,1)
    ok = rgasp_toolbox(rows{i,1}, true);
    if ok
        st = 'available';
    else
        st = ['MISSING -> ' rows{i,4}];
        nmiss = nmiss + 1;
    end
    fprintf('%-14s %-13s %-42s %s\n', rows{i,1}, rows{i,2}, rows{i,3}, st);
end
fprintf('%s\n', repmat('-', 1, 92));
if nmiss == 0
    fprintf('All optional toolbox functions are available and will be used.\n');
else
    fprintf('%d function(s) unavailable; the listed fallbacks will be used instead.\n', nmiss);
end
if rgasp_have_mex()
    fprintf('C++ acceleration: compiled and active.\n\n');
else
    fprintf('C++ acceleration: not compiled (run rgasp_compile_mex).\n\n');
end
end
