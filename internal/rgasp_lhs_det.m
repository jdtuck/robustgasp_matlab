function X = rgasp_lhs_det(n, p, seed, num_try)
%RGASP_LHS_DET  Deterministic maximin Latin hypercube design on [0,1]^p.
%
%   X = RGASP_LHS_DET(n, p)
%   X = RGASP_LHS_DET(n, p, seed, num_try)
%
%   Used by the test suite, where a design that is byte-identical in every
%   MATLAB and Octave version is what makes the recorded results reproducible.
%   RGASP_LHS is the user-facing generator and prefers LHSDESIGN.
%
%   Generates num_try random LHDs (points at cell centres) and keeps the one
%   with the largest minimum inter-point distance.  Latin hypercube designs
%   are the design class the robustness results of Gu, Wang and Berger (2018)
%   are stated for.  A small self-contained linear congruential generator is
%   used so that demos and tests reproduce identically in MATLAB and Octave.

if nargin < 3 || isempty(seed), seed = 0; end
if nargin < 4 || isempty(num_try), num_try = 20; end

U = rgasp_lcg(n*p*num_try, seed);
best = []; best_d = -Inf;
pos = 0;
for it = 1:num_try
    Xc = zeros(n, p);
    for l = 1:p
        u = U(pos + (1:n)); pos = pos + n;
        [~, perm] = sort(u);
        Xc(perm, l) = ((1:n)' - 0.5) / n;
    end
    d = min_pairwise_dist(Xc);
    if d > best_d
        best_d = d; best = Xc;
    end
end
X = best;
end

function d = min_pairwise_dist(X)
n = size(X,1);
d = Inf;
for i = 1:n-1
    df = bsxfun(@minus, X(i+1:end,:), X(i,:));
    d = min(d, min(sum(df.^2, 2)));
end
d = sqrt(d);
end
