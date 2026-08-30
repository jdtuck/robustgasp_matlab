function u = rgasp_lcg(n, seed)
%RGASP_LCG  Reproducible uniform(0,1) stream (Numerical Recipes LCG).
%
%   u = RGASP_LCG(n, seed) returns an n-by-1 vector.  Used by the demos and
%   tests so that results are identical across MATLAB and Octave versions;
%   it is not intended for serious Monte Carlo work.

if nargin < 2 || isempty(seed), seed = 0; end
m = 2^32;
s = mod(floor(seed)*2654435761 + 1013904223, m);
u = zeros(n, 1);
for i = 1:n
    s = mod(1664525*s + 1013904223, m);
    u(i) = s / m;
end
end
