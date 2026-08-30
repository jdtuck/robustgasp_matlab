function y = environ_4_data(x, s, t)
%ENVIRON_4_DATA  Environmental spill model (Bliznyuk et al., 2008).
%
%   y = ENVIRON_4_DATA(x)          uses the default space/time grid
%   y = ENVIRON_4_DATA(x, s, t)    uses the supplied grids
%
%   x is n-by-4 with columns (M, D, L, tau):
%       M    mass of pollutant       [7, 13]
%       D    diffusion rate          [0.02, 0.12]
%       L    location of 2nd spill   [0.01, 3]
%       tau  time of 2nd spill       [30.01, 30.295]
%
%   Returns an n-by-(ns*nt) matrix of concentrations C(s,t), i.e. a vector
%   ("field") output per run -- the natural input for PPGASP.
if nargin < 2 || isempty(s), s = [0.5 1 1.5 2 2.5]; end
if nargin < 3 || isempty(t), t = 0.3:0.3:12; end
x = reshape(x, [], 4);
n = size(x,1);
[S, T] = meshgrid(s, t);          % nt-by-ns
S = S(:)'; T = T(:)';             % 1-by-(ns*nt)
y = zeros(n, numel(S));
for i = 1:n
    M = x(i,1); D = x(i,2); L = x(i,3); tau = x(i,4);
    term1 = M ./ sqrt(4*pi*D*T) .* exp(-S.^2 ./ (4*D*T));
    term2 = zeros(size(T));
    idx = T > tau;
    term2(idx) = M ./ sqrt(4*pi*D*(T(idx)-tau)) .* ...
                 exp(-(S(idx)-L).^2 ./ (4*D*(T(idx)-tau)));
    y(i,:) = sqrt(4*pi) * (term1 + term2);
end
end
