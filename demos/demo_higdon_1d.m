function demo_higdon_1d(savefig_dir)
%DEMO_HIGDON_1D  One-dimensional emulation of the Higdon (2002) function.
%
%   demo_higdon_1d()            show the figure
%   demo_higdon_1d(outdir)      also save it as PNG in outdir
%
%   Fits an RGASP emulator to 15 equally spaced runs, plots the predictive
%   mean and the 95% predictive band, and reports leave-one-out diagnostics.

if nargin < 1, savefig_dir = ''; end

x = linspace(0, 10, 15)';
y = higdon_1_data(x);

model = rgasp(x, y);
rgasp_summary(model);

xt = linspace(0, 10, 400)';
yt = higdon_1_data(xt);
pred = rgasp_predict(model, xt);

fprintf('Out-of-sample criteria on a dense grid:\n');
rgasp_validate(pred, yt);

loo = rgasp_loo(model);
fprintf('\nLeave-one-out RMSE %.5g, 95%% coverage %.3f\n', loo.rmse, loo.coverage);

fig = figure('Visible', 'off');
subplot(2,1,1);
fill([xt; flipud(xt)], [pred.lower95; flipud(pred.upper95)], [0.85 0.90 0.97], ...
     'EdgeColor', 'none'); hold on
plot(xt, yt, 'k-', 'LineWidth', 1.2);
plot(xt, pred.mean, 'b--', 'LineWidth', 1.4);
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
legend('95% predictive band', 'truth', 'emulator mean', 'design points', ...
       'Location', 'southwest');
xlabel('s'); ylabel('y(s)');
title(sprintf('RGASP emulator, Higdon (2002) function (n = %d, gamma = %.3f)', ...
      numel(x), model.range_hat));
box on

subplot(2,1,2);
errorbar(y, loo.mean, 2*loo.sd, 'o'); hold on
lims = [min(y) max(y)];
plot(lims, lims, 'k-');
xlabel('held-out output'); ylabel('leave-one-out prediction');
title('Leave-one-out cross validation');
box on

if ~isempty(savefig_dir)
    f = fullfile(savefig_dir, 'demo_higdon_1d.png');
    print(fig, f, '-dpng', '-r120');
    fprintf('Figure written to %s\n', f);
else
    set(fig, 'Visible', 'on');
end
end
