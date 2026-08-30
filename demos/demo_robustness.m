function demo_robustness(savefig_dir)
%DEMO_ROBUSTNESS  Reproduce the central comparison of arXiv:1708.04738.
%
%   Panel 1 -- the shape of the objective.  Along the ray xi = log(beta) the
%   profile likelihood and the marginal likelihood become perfectly FLAT as
%   the range parameters shrink to zero, so their maximizer can sit anywhere
%   on that plateau; there the correlation matrix degenerates to the identity
%   and the emulator collapses to its fitted mean.  The jointly robust prior
%   adds a term -b*sum(C_l beta_l) that grows without bound, so the marginal
%   posterior mode is always interior.
%
%   Panel 2 -- consequences.  Over repeated Latin hypercube designs for the
%   Lim et al. (2002) function embedded in four inputs (two of which are
%   inert), the profile likelihood occasionally produces a useless emulator
%   while the marginal posterior mode never does.

if nargin < 1, savefig_dir = ''; end

% ------------------------------------------------ objective along a ray
U = rgasp_lhs(16, 4, 301);
y = limetal_2_data(U(:,1:2));
n = 16; p = 4;
core = struct('R0', {rgasp_R0(U)}, 'X', ones(n,1), 'Y', y, ...
    'kernel_type', {repmat({'matern_5_2'},1,p)}, 'alpha', 1.9*ones(1,p), ...
    'nugget', 0, 'nugget_est', false, 'method', 'mle', ...
    'CL', ((max(U)-min(U))'/n^(1/p)), 'a', 0.2, 'b', (0.2+p)/n^(1/p));

xis = linspace(-6, 10, 90);
F = nan(numel(xis), 3);
meths = {'mle','mmle','post_mode'};
for i = 1:numel(xis)
    for j = 1:3
        core.method = meths{j};
        v = rgasp_objective(xis(i)*ones(p,1), core);
        if v < 1e99, F(i,j) = v; end
    end
end

fprintf('Objective at xi = 5, 8, 11, 14 (log inverse range):\n');
for j = 1:3
    core.method = meths{j};
    vals = arrayfun(@(z) rgasp_objective(z*ones(p,1), core), [5 8 11 14]);
    fprintf('  %-10s %s\n', meths{j}, sprintf('%14.6g', vals));
end
fprintf(['  -> the two likelihood criteria are constant to machine precision;\n' ...
         '     the marginal posterior grows without bound.\n\n']);

% ------------------------------------------------ repeated designs
%   Whether a particular solver actually walks onto the plateau depends on the
%   solver, so both are shown: the built-in projected L-BFGS (identical in
%   every installation) and, when the Optimization Toolbox is present, its
%   solver.  The plateau itself -- panel 1 -- is a property of the objective
%   and does not depend on either.
Ut = rgasp_lhs(400, 4, 999);
yt = limetal_2_data(Ut(:,1:2));
nrep = 20; nd = 16;

solvers = {'lbfgs'};
labels  = {'built-in projected L-BFGS'};
if rgasp_toolbox('fmincon') || rgasp_toolbox('fminunc')
    solvers{end+1} = 'auto';
    labels{end+1}  = 'Optimization Toolbox solver';
end

res_all = cell(1, numel(solvers));
for s = 1:numel(solvers)
    res = zeros(nrep, 3);
    for r = 1:nrep
        Ur = rgasp_lhs(nd, 4, 300+r);
        yr = limetal_2_data(Ur(:,1:2));
        for j = 1:3
            mj = rgasp(Ur, yr, 'method', meths{j}, 'lowerBound', false, ...
                       'numInitialValues', 3, 'optimizer', solvers{s});
            pj = rgasp_predict(mj, Ut);
            sj = rgasp_validate(pj, yt, false);
            res(r,j) = sj.nrmse;
        end
    end
    res_all{s} = res;
    fprintf('Normalized out-of-sample RMSE over %d designs (%s):\n', nrep, labels{s});
    fprintf('              %10s %10s %10s\n', meths{:});
    fprintf('  median      %10.4f %10.4f %10.4f\n', median(res));
    fprintf('  mean        %10.4f %10.4f %10.4f\n', mean(res));
    fprintf('  worst case  %10.4f %10.4f %10.4f\n', max(res));
    fprintf('  # > 0.3     %10d %10d %10d\n\n', sum(res > 0.3, 1));
end
res = res_all{1};

fig = figure('Visible','off');
subplot(1,2,1);
plot(xis, F(:,1), 'r-', 'LineWidth', 1.3); hold on
plot(xis, F(:,2), 'g-.', 'LineWidth', 1.3);
plot(xis, F(:,3), 'b--', 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
xlabel('\xi = log(1/\gamma), common to all four inputs');
ylabel('negative log objective');
legend({'profile likelihood', 'marginal likelihood', ...
        'marginal posterior (JR prior)'}, 'Location', 'northwest');
title('Flat likelihood plateau as \gamma \rightarrow 0'); box on

subplot(1,2,2);
boxplot_fallback(res, {'MLE','MMLE','post. mode'});
ylabel('normalized out-of-sample RMSE');
title(sprintf('%d repeated designs (%s)', nrep, 'built-in L-BFGS')); box on

if ~isempty(savefig_dir)
    f = fullfile(savefig_dir, 'demo_robustness.png');
    print(fig, f, '-dpng', '-r120');
    fprintf('Figure written to %s\n', f);
else
    set(fig, 'Visible', 'on');
end
end

% ----------------------------------------------------------------------
function boxplot_fallback(res, labels)
%   BOXPLOT from the Statistics Toolbox when available, otherwise a minimal
%   hand-drawn equivalent so the demo still runs.
if rgasp_toolbox('boxplot')
    boxplot(res, 'Labels', labels);
    return
end
hold on
for j = 1:size(res,2)
    v = sort(res(:,j));
    q1 = rgasp_quantile(v, 0.25);
    q2 = rgasp_quantile(v, 0.50);
    q3 = rgasp_quantile(v, 0.75);
    w = 0.28;
    plot([j-w j+w j+w j-w j-w], [q1 q1 q3 q3 q1], 'k-', 'LineWidth', 1.2);
    plot([j-w j+w], [q2 q2], 'b-', 'LineWidth', 1.6);
    plot([j j], [min(v) q1], 'k-'); plot([j j], [q3 max(v)], 'k-');
    plot(j + 0.05*randn(size(v)), v, 'r.', 'MarkerSize', 7);
end
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels);
xlim([0.4, numel(labels)+0.6]);
end
