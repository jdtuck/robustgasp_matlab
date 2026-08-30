function stats = rgasp_validate(pred, y_true, verbose)
%RGASP_VALIDATE  Out-of-sample criteria used in Section 5 of the paper.
%
%   stats = RGASP_VALIDATE(pred, y_true)
%
%   pred    output of RGASP_PREDICT / PPGASP_PREDICT
%   y_true  held-out outputs, same size as pred.mean
%
%   Returns
%     .rmse       root mean squared error
%     .nrmse      RMSE normalized by the standard deviation of y_true
%     .coverage   proportion of held-out values inside the 95% interval
%     .length     average length of the 95% interval
%     .p_ci       normalized average interval length,
%                 mean(length) / sd(y_true)
%
%   Gu, Wang and Berger (2018) compare emulators using exactly these three
%   quantities: predictions should have small RMSE, coverage close to the
%   nominal 95%, and short intervals.

if nargin < 3 || isempty(verbose), verbose = true; end

m = pred.mean;
y_true = reshape(y_true, size(m));

err = y_true - m;
stats.rmse = sqrt(mean(err(:).^2));
sy = std(y_true(:), 1);
stats.nrmse = stats.rmse / sy;

inside = (y_true >= pred.lower95) & (y_true <= pred.upper95);
stats.coverage = mean(inside(:));
len = pred.upper95 - pred.lower95;
stats.length = mean(len(:));
stats.p_ci = stats.length / sy;

if verbose
    fprintf('RMSE            : %.6g\n', stats.rmse);
    fprintf('Normalized RMSE : %.6g\n', stats.nrmse);
    fprintf('95%% coverage    : %.4f\n', stats.coverage);
    fprintf('Avg 95%% length  : %.6g\n', stats.length);
end
end
