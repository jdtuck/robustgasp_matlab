function test_paper_behaviour()
%TEST_PAPER_BEHAVIOUR  Reproduce the two central claims of arXiv:1708.04738.
%
%   (i)  The profile likelihood and the marginal likelihood are FLAT as the
%        range parameters go to zero (xi = log(1/gamma) -> +Inf), so their
%        maximizers can sit on that degenerate plateau, where the emulator
%        reverts to the fitted mean.  Adding the jointly robust prior makes
%        the objective increase without bound in that direction, so the
%        marginal posterior mode is always interior.
%
%   (ii) Over repeated Latin hypercube designs the marginal posterior mode
%        never produces the catastrophic predictions the profile likelihood
%        occasionally produces (Section 5 of the paper).

% ---------------------------------------------------------------- (i)
U = rgasp_lhs_det(16, 4, 301);
y = limetal_2_data(U(:,1:2));    % coordinates 3 and 4 are inert
n = 16; p = 4;
core = struct('R0', {rgasp_R0(U)}, 'X', ones(n,1), 'Y', y, ...
    'kernel_type', {repmat({'matern_5_2'},1,p)}, 'alpha', 1.9*ones(1,p), ...
    'nugget', 0, 'nugget_est', false, 'method', 'mle', ...
    'CL', ((max(U)-min(U))'/n^(1/p)), 'a', 0.2, 'b', (0.2+p)/n^(1/p));

xis = [5 8 11 14];
f_mle = zeros(size(xis)); f_mmle = f_mle; f_post = f_mle;
for i = 1:numel(xis)
    xi = xis(i)*ones(p,1);
    core.method = 'mle';       f_mle(i)  = rgasp_objective(xi, core);
    core.method = 'mmle';      f_mmle(i) = rgasp_objective(xi, core);
    core.method = 'post_mode'; f_post(i) = rgasp_objective(xi, core);
end

if max(abs(diff(f_mle))) > 1e-6
    error(['the profile likelihood should be flat as gamma -> 0 ' ...
           '(values %s)'], mat2str(f_mle, 8));
end
if max(abs(diff(f_mmle))) > 1e-6
    error('the marginal likelihood should be flat as gamma -> 0 (%s)', ...
          mat2str(f_mmle, 8));
end
if ~all(diff(f_post) > 0) || f_post(end) < 100*f_post(1)
    error('the jointly robust prior must penalize gamma -> 0 (%s)', ...
          mat2str(f_post, 6));
end

% the mode really is interior: the objective increases in both directions
m = rgasp(U, y, 'numInitialValues', 3);
xi_hat = log(m.beta_hat);
core.method = 'post_mode';
f_hat = rgasp_objective(xi_hat, core);
for delta = [-2 -1 1 2]
    fd = rgasp_objective(xi_hat + delta, core);
    if fd <= f_hat
        error('objective not minimized at the reported mode (delta = %g)', delta);
    end
end
[~, ghat] = rgasp_objective(xi_hat, core);
active = (xi_hat <= m.LB(1:p) + 1e-10);
if max(abs(ghat(~active))) > 1e-3 * max(1, abs(f_hat))
    error('gradient at the mode is not small: %s', mat2str(ghat', 4));
end

% --------------------------------------------------------------- (ii)
%   The plateau makes the likelihood maximizer ill-posed, but whether a given
%   solver actually walks onto it is a property of the solver.  The assertion
%   below is therefore made with the toolbox-independent built-in optimizer,
%   which is identical everywhere; the Optimization Toolbox solver is reported
%   alongside it when present, and only has to keep post_mode out of trouble.
Ut = rgasp_lhs_det(400, 4, 999);
yt = limetal_2_data(Ut(:,1:2));
nrep = 20; nd = 16;
names = {'mle','mmle','post_mode'};

solvers = {'lbfgs'};
if rgasp_toolbox('fminunc') || rgasp_toolbox('fmincon')
    solvers{end+1} = 'auto';
end

res_all = cell(1, numel(solvers));
for s = 1:numel(solvers)
    res = zeros(nrep, 3);
    for r = 1:nrep
        Ur = rgasp_lhs_det(nd, 4, 300+r);
        yr = limetal_2_data(Ur(:,1:2));
        for j = 1:3
            mj = rgasp(Ur, yr, 'method', names{j}, 'lowerBound', false, ...
                       'numInitialValues', 3, 'optimizer', solvers{s});
            pj = rgasp_predict(mj, Ut);
            sj = rgasp_validate(pj, yt, false);
            res(r,j) = sj.nrmse;
        end
    end
    res_all{s} = res;

    if strcmp(solvers{s}, 'auto')
        label = 'Optimization Toolbox solver';
    else
        label = 'built-in projected L-BFGS';
    end
    fprintf('  normalized RMSE over %d designs, %s:\n', nrep, label);
    fprintf('             %10s %10s %10s\n', names{:});
    fprintf('    median   %10.4f %10.4f %10.4f\n', median(res));
    fprintf('    mean     %10.4f %10.4f %10.4f\n', mean(res));
    fprintf('    max      %10.4f %10.4f %10.4f\n', max(res));
    fprintf('    failures %10d %10d %10d\n', sum(res > 0.3, 1));
end

res = res_all{1};
fail = sum(res > 0.3, 1);
if fail(3) > 0
    error('the marginal posterior mode produced %d degenerate fits', fail(3));
end
if fail(1) == 0
    error(['the profile likelihood was expected to degenerate on at least ' ...
           'one design; the comparison is not exercising the failure mode']);
end
if mean(res(:,3)) > 0.5*mean(res(:,1))
    error('posterior mode is not clearly better than the profile likelihood');
end

for s = 2:numel(res_all)
    r2 = res_all{s};
    if sum(r2(:,3) > 0.3) > 0
        error('posterior mode degenerated under the toolbox solver');
    end
    if mean(r2(:,3)) > mean(r2(:,1))
        error('posterior mode is worse than the profile likelihood under the toolbox solver');
    end
end

fprintf('paper behaviour: flat likelihood plateau, interior posterior mode\n');
fprintf('                 and robustness over repeated designs reproduced\n');
end
