function z = rgasp_norminv(p, mu, sigma)
%RGASP_NORMINV  Standard normal quantile function.
%
%   Uses NORMINV from the Statistics and Machine Learning Toolbox when it is
%   available; the fallback is the equivalent erfinv expression.

if nargin < 2 || isempty(mu), mu = 0; end
if nargin < 3 || isempty(sigma), sigma = 1; end

if rgasp_toolbox('norminv')
    z = norminv(p, mu, sigma);
else
    z = mu + sigma * sqrt(2) * erfinv(2*p - 1);
end
end
