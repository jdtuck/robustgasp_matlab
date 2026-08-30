function demo_noisy_nugget(savefig_dir)
%DEMO_NOISY_NUGGET  GaSP with a noise term (Section 4 of the paper).
%
%   Observations y = f(x) + epsilon are emulated with the covariance
%   sigma^2 (c(.,.) + eta I).  The nugget-variance ratio eta is estimated
%   jointly with the range parameters from the marginal posterior, using the
%   jointly robust prior on (beta, eta).  The demo compares the fit with and
%   without the noise term.

if nargin < 1, savefig_dir = ''; end

n = 60;
x = 10*rgasp_lcg(n, 21);
u1 = rgasp_lcg(n, 3); u2 = rgasp_lcg(n, 5);
noise = 0.1 * sqrt(-2*log(u1)) .* cos(2*pi*u2);      % Box-Muller, sd = 0.1
y = higdon_1_data(x) + noise;

fprintf('=== interpolator (eta fixed at 0) ===\n');
m0 = rgasp(x, y);
rgasp_summary(m0);

fprintf('=== noise term estimated ===\n');
m1 = rgasp(x, y, 'nuggetEst', true, 'numInitialValues', 3);
rgasp_summary(m1);
fprintf('Estimated noise variance  eta*sigma^2 = %.5g   (true 0.01)\n', ...
        m1.nugget * m1.sigma2_hat);

xt = linspace(0, 10, 400)';
yt = higdon_1_data(xt);
p0 = rgasp_predict(m0, xt);
p1 = rgasp_predict(m1, xt);

fprintf('\nAccuracy against the noise-free truth:\n');
fprintf('  no nugget    : '); s0 = rgasp_validate(p0, yt, false);
fprintf('RMSE %.5g, coverage %.3f\n', s0.rmse, s0.coverage);
fprintf('  nugget est.  : '); s1 = rgasp_validate(p1, yt, false);
fprintf('RMSE %.5g, coverage %.3f\n', s1.rmse, s1.coverage);

fig = figure('Visible','off');
ttl = {'eta fixed at 0 (interpolation)', 'eta estimated (smoothing)'};
P = {p0, p1};
for j = 1:2
    subplot(2,1,j);
    fill([xt; flipud(xt)], [P{j}.lower95; flipud(P{j}.upper95)], ...
         [0.87 0.91 0.97], 'EdgeColor','none'); hold on
    plot(xt, yt, 'k-', 'LineWidth', 1.1);
    plot(xt, P{j}.mean, 'b--', 'LineWidth', 1.3);
    plot(x, y, 'r.', 'MarkerSize', 9);
    xlabel('s'); ylabel('y'); title(ttl{j}); box on
    ylim([min(y)-0.4, max(y)+0.4]);
end

if ~isempty(savefig_dir)
    f = fullfile(savefig_dir, 'demo_noisy_nugget.png');
    print(fig, f, '-dpng', '-r120');
    fprintf('Figure written to %s\n', f);
else
    set(fig, 'Visible', 'on');
end
end
