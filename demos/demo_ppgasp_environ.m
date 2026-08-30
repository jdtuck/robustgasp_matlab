function demo_ppgasp_environ(savefig_dir)
%DEMO_PPGASP_ENVIRON  PP GaSP emulation of a vector-valued simulator.
%
%   The environmental spill model of Bliznyuk et al. (2008) returns a
%   concentration field C(s,t) on a 5 x 40 space-time grid (200 outputs) for
%   every choice of the four inputs (M, D, L, tau).  PPGASP shares one set of
%   range parameters across all 200 output coordinates while giving each its
%   own mean and variance, so a single 40-run design suffices.

if nargin < 1, savefig_dir = ''; end

lb = [7,   0.02, 0.01, 30.010];
ub = [13,  0.12, 3.00, 30.295];
nd = 40; nt = 50;

Ud = rgasp_lhs(nd, 4, 2);
Ut = rgasp_lhs(nt, 4, 77);
Dd = bsxfun(@plus, bsxfun(@times, Ud, ub-lb), lb);
Dt = bsxfun(@plus, bsxfun(@times, Ut, ub-lb), lb);
Yd = environ_4_data(Dd);
Yt = environ_4_data(Dt);

fprintf('Training runs: %d,  outputs per run: %d\n', size(Yd,1), size(Yd,2));
model = ppgasp(Dd, Yd, 'numInitialValues', 3);
rgasp_summary(model);

pred = ppgasp_predict(model, Dt);
fprintf('Held-out validation over all %d x %d values:\n', nt, size(Yt,2));
rgasp_validate(pred, Yt);

s = [0.5 1 1.5 2 2.5];
t = 0.3:0.3:12;
idx = 3;                          % show one held-out run
Ytrue = reshape(Yt(idx,:), numel(t), numel(s));
Ypred = reshape(pred.mean(idx,:), numel(t), numel(s));
Ylo   = reshape(pred.lower95(idx,:), numel(t), numel(s));
Yhi   = reshape(pred.upper95(idx,:), numel(t), numel(s));

fig = figure('Visible','off');
for j = 1:numel(s)
    subplot(2,3,j);
    fill([t'; flipud(t')], [Ylo(:,j); flipud(Yhi(:,j))], [0.87 0.91 0.97], ...
         'EdgeColor','none'); hold on
    plot(t, Ytrue(:,j), 'k-', 'LineWidth', 1.1);
    plot(t, Ypred(:,j), 'b--', 'LineWidth', 1.2);
    xlabel('t'); ylabel('C(s,t)');
    title(sprintf('s = %.1f', s(j)));
    box on
end
subplot(2,3,6);
plot(Yt(:), pred.mean(:), '.', 'MarkerSize', 4); hold on
lims = [min(Yt(:)) max(Yt(:))];
plot(lims, lims, 'k-');
xlabel('simulator'); ylabel('PP GaSP');
title('All held-out outputs'); box on

if ~isempty(savefig_dir)
    f = fullfile(savefig_dir, 'demo_ppgasp_environ.png');
    print(fig, f, '-dpng', '-r120');
    fprintf('Figure written to %s\n', f);
else
    set(fig, 'Visible', 'on');
end
end
