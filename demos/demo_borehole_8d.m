function demo_borehole_8d(savefig_dir)
%DEMO_BOREHOLE_8D  Emulating the 8-dimensional borehole function.
%
%   Fits an RGASP emulator on a maximin Latin hypercube of 80 runs, validates
%   it on 500 held-out runs with the criteria used in Section 5 of the paper
%   (RMSE, 95% coverage, average interval length), and runs the inert-input
%   diagnostic.  The borehole function is dominated by rw; r and Tl are
%   nearly inactive.

if nargin < 1, savefig_dir = ''; end

[lb, ub] = borehole_ranges();
n = 80; nt = 500;

Ud = rgasp_lhs(n, 8, 11);
Ut = rgasp_lhs(nt, 8, 999);
Dd = bsxfun(@plus, bsxfun(@times, Ud, ub-lb), lb);
Dt = bsxfun(@plus, bsxfun(@times, Ut, ub-lb), lb);
yd = borehole(Dd);
yt = borehole(Dt);

model = rgasp(Dd, yd, 'numInitialValues', 3);
rgasp_summary(model);

pred = rgasp_predict(model, Dt);
fprintf('Held-out validation (%d points):\n', nt);
rgasp_validate(pred, yt);

fprintf('\nInert-input diagnostic:\n');
P = rgasp_inert_inputs(model);

names = {'rw','r','Tu','Hu','Tl','Hl','L','Kw'};
fig = figure('Visible','off');
subplot(1,2,1);
plot(yt, pred.mean, '.', 'MarkerSize', 8); hold on
lims = [min(yt) max(yt)];
plot(lims, lims, 'k-');
xlabel('borehole output'); ylabel('emulator prediction');
title(sprintf('Held-out predictions (RMSE = %.3g)', ...
      sqrt(mean((yt-pred.mean).^2))));
box on

subplot(1,2,2);
bar(P); hold on
plot([0 9], [0.1 0.1], 'r--');
set(gca, 'XTick', 1:8, 'XTickLabel', names);
ylabel('normalized inverse range P_l');
title('Inert-input diagnostic (threshold 0.1)');
box on

if ~isempty(savefig_dir)
    f = fullfile(savefig_dir, 'demo_borehole_8d.png');
    print(fig, f, '-dpng', '-r120');
    fprintf('Figure written to %s\n', f);
else
    set(fig, 'Visible', 'on');
end
end
