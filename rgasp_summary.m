function rgasp_summary(model)
%RGASP_SUMMARY  Print a summary of a fitted (pp)GaSP emulator.
%
%   RGASP_SUMMARY(model)

fprintf('\n');
if strcmp(model.type,'rgasp')
    fprintf('Robust GaSP emulator (scalar output)\n');
else
    fprintf('Parallel partial GaSP emulator (%d output coordinates)\n', model.k);
end
fprintf('------------------------------------------------------------\n');
fprintf('Number of design points        : %d\n', model.num_obs);
fprintf('Input dimension                : %d\n', model.p);
kt = model.kernel_type;
fprintf('Correlation function           : %s%s\n', kt{1}, ...
        repmat_note(kt));
if any(strcmpi(kt,'pow_exp'))
    fprintf('Roughness alpha                : %s\n', num2str(model.alpha, ' %.3g'));
end
fprintf('Estimation method              : %s\n', model.method);
if strcmpi(model.method,'post_mode')
    fprintf('Prior                          : jointly robust (a=%.3g, b=%.4g)\n', model.a, model.b);
end
fprintf('Mean function                  : %s\n', ...
        tern(model.zero_mean, 'zero mean', sprintf('%d basis function(s)', model.q)));
fprintf('log marginal posterior at mode : %.6f\n', model.log_post);
fprintf('\nEstimated range parameters (gamma = 1/beta):\n');
fprintf('  %12.6g', model.range_hat); fprintf('\n');
fprintf('Estimated inverse ranges (beta):\n');
fprintf('  %12.6g', model.beta_hat); fprintf('\n');
fprintf('Nugget-variance ratio (eta)    : %.6g%s\n', model.nugget, ...
        tern(model.nugget_est, ' (estimated)', ' (fixed)'));
if numel(model.sigma2_hat) == 1
    fprintf('Estimated variance sigma^2     : %.6g\n', model.sigma2_hat);
else
    fprintf('Estimated variances sigma_i^2  : min %.4g / median %.4g / max %.4g\n', ...
        min(model.sigma2_hat), median(model.sigma2_hat), max(model.sigma2_hat));
end
if model.q > 0
    fprintf('Trend coefficients (theta_hat) :\n');
    if size(model.theta_hat,2) == 1
        fprintf('  %12.6g', model.theta_hat); fprintf('\n');
    else
        fprintf('  [%d-by-%d matrix]\n', size(model.theta_hat,1), size(model.theta_hat,2));
    end
end
fprintf('------------------------------------------------------------\n\n');
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function s = repmat_note(kt)
if numel(unique(kt)) > 1
    s = ' (mixed across dimensions)';
else
    s = '';
end
end
