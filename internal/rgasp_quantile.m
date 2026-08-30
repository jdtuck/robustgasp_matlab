function q = rgasp_quantile(v, pr)
%RGASP_QUANTILE  Sample quantile.
%
%   Uses QUANTILE from the Statistics and Machine Learning Toolbox when it is
%   available.  The fallback reproduces the same definition MATLAB uses
%   (Hyndman and Fan type 5): the sorted values are placed at cumulative
%   probabilities (i - 0.5)/n, linearly interpolated in between, and clamped
%   to the extremes outside that range.

if rgasp_toolbox('quantile')
    q = quantile(v(:), pr);
    return
end

v = sort(v(:));
n = numel(v);
q = zeros(size(pr));
pos = ((1:n)' - 0.5) / n;
for i = 1:numel(pr)
    p = pr(i);
    if p <= pos(1)
        q(i) = v(1);
    elseif p >= pos(end)
        q(i) = v(end);
    else
        j = find(pos <= p, 1, 'last');
        w = (p - pos(j)) / (pos(j+1) - pos(j));
        q(i) = v(j) + w*(v(j+1) - v(j));
    end
end
end
