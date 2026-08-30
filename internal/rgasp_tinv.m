function t = rgasp_tinv(p, nu)
%RGASP_TINV  Quantile function of the Student t distribution.
%
%   Uses TINV from the Statistics and Machine Learning Toolbox when it is
%   available.
%
%   The fallback inverts the t CDF directly.  With
%       F(t) = 1 - 0.5*I_{nu/(nu+t^2)}(nu/2, 1/2)   for t >= 0,
%   the quantile is found by bracketing and bisection followed by a few
%   Newton steps on the density.  (Going through BETAINCINV in one shot is
%   tempting but is not reliable in the far tails -- for nu = 38 and
%   p = 0.001 it can return a root whose actual tail probability is 0.015.)

if any(p(:) <= 0 | p(:) >= 1)
    error('rgasp:tinv', 'p must lie strictly between 0 and 1.');
end
p = double(p);

if rgasp_toolbox('tinv')
    t = tinv(p, nu);
    return
end

t = zeros(size(p));
for i = 1:numel(p)
    pi_ = p(i);
    if pi_ == 0.5
        t(i) = 0;
        continue
    end
    pu = max(pi_, 1 - pi_);          % upper-tail probability, >= 0.5
    val = upper_quantile(pu, nu);
    if pi_ < 0.5, t(i) = -val; else, t(i) = val; end
end
end

% ----------------------------------------------------------------------
function x = upper_quantile(pu, nu)
%   Solve F(x) = pu for x >= 0.

% bracket
hi = 1;
while tcdf_local(hi, nu) < pu
    hi = hi * 2;
    if hi > 1e12
        x = hi; return
    end
end
lo = 0;

% bisection to a safe bracket
for it = 1:200
    mid = 0.5*(lo + hi);
    if tcdf_local(mid, nu) < pu, lo = mid; else, hi = mid; end
    if hi - lo < 1e-13*max(1, hi), break, end
end
x = 0.5*(lo + hi);

% Newton polish using the exact density
logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(nu*pi);
for it = 1:5
    fx = tcdf_local(x, nu) - pu;
    dens = exp(logc - ((nu+1)/2)*log1p(x*x/nu));
    if dens <= 0 || ~isfinite(dens), break, end
    step = fx/dens;
    xn = x - step;
    if ~isfinite(xn) || xn < lo || xn > hi, break, end
    x = xn;
    if abs(step) < 1e-14*max(1, abs(x)), break, end
end
end

function F = tcdf_local(x, nu)
%   CDF of the t distribution via the regularized incomplete beta function.
z = nu / (nu + x*x);
tail = 0.5 * betainc(z, nu/2, 0.5);      % P(T > |x|)
if x >= 0
    F = 1 - tail;
else
    F = tail;
end
end
