function loo = rgasp_loo(model)
%RGASP_LOO  Leave-one-out cross validation for a fitted RGASP emulator.
%
%   loo = RGASP_LOO(model) returns a struct with
%       .mean   n-by-1 leave-one-out predictive means
%       .sd     n-by-1 leave-one-out predictive standard deviations
%       .resid  observed minus predicted
%       .rmse   root mean squared LOO error
%       .coverage / .length  empirical coverage and mean length of the
%                            nominal 95% LOO intervals
%
%   The range parameters are held fixed at their full-data estimates (as in
%   RobustGaSP's plot/leave_one_out_rgasp); only the mean, the variance and
%   the interpolation weights are refitted on each subset.
%
%   See also RGASP, RGASP_VALIDATE.

if ~strcmp(model.type, 'rgasp')
    error('rgasp:type', 'rgasp_loo supports scalar-output models only.');
end

n = model.num_obs;
q = model.q;
y = model.output;
Rt = model.L * model.L';           % R(beta_hat) + eta I

mu = zeros(n,1);
s2 = zeros(n,1);

for i = 1:n
    idx = [1:i-1, i+1:n];
    r_sub  = Rt(idx, i);
    R_sub  = Rt(idx, idx);
    U      = chol(R_sub);
    Rinv_r = U \ (U' \ r_sub);                    % R_sub^{-1} r
    y_sub  = y(idx);

    if q == 0
        mu(i) = Rinv_r' * y_sub;
        sigma2_hat = (y_sub' * (U \ (U' \ y_sub))) / (n - 1);
        c_star = Rt(i,i) - Rinv_r' * r_sub;
        s2(i)  = sigma2_hat * c_star;
    else
        X_sub  = model.X(idx, :);
        Rinv_X = U \ (U' \ X_sub);
        A      = X_sub' * Rinv_X;
        UX     = chol((A+A')/2);
        theta  = UX \ (UX' \ (Rinv_X' * y_sub));
        tilde  = y_sub - X_sub * theta;
        mu(i)  = model.X(i,:) * theta + Rinv_r' * tilde;

        c_star = Rt(i,i) - Rinv_r' * r_sub;
        if strcmpi(model.method, 'mle')
            sigma2_hat = (tilde' * (U \ (U' \ tilde))) / (n - 1);
            s2(i) = sigma2_hat * c_star;
        else
            sigma2_hat = (tilde' * (U \ (U' \ tilde))) / (n - 1 - q);
            h_hat = model.X(i,:)' - X_sub' * Rinv_r;
            c_ss  = c_star + h_hat' * (UX \ (UX' \ h_hat));
            s2(i) = sigma2_hat * c_ss;
        end
    end
end

s2 = abs(s2);
loo = struct();
loo.mean  = mu;
loo.sd    = sqrt(s2);
loo.resid = y - mu;
loo.rmse  = sqrt(mean(loo.resid.^2));

df = n - 1 - q;
if strcmpi(model.method,'mle'), df = n - 1; end
if df > 0
    tq = rgasp_tinv(0.975, df);
else
    tq = rgasp_norminv(0.975);
end
loo.lower95 = mu - tq*loo.sd;
loo.upper95 = mu + tq*loo.sd;
loo.coverage = mean(y >= loo.lower95 & y <= loo.upper95);
loo.length   = mean(loo.upper95 - loo.lower95);
loo.standardized_resid = loo.resid ./ loo.sd;
end
