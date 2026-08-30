function X = rgasp_lhs(n, p, seed, num_try)
%RGASP_LHS  Maximin Latin hypercube design on [0,1]^p.
%
%   X = RGASP_LHS(n, p)
%   X = RGASP_LHS(n, p, seed)
%   X = RGASP_LHS(n, p, seed, num_try)
%
%   Uses LHSDESIGN from the Statistics and Machine Learning Toolbox with the
%   'maximin' criterion.  Latin hypercube designs are the design class the
%   robustness results of Gu, Wang and Berger (2018) are stated for.
%
%   With a seed supplied the random stream is set first, so a given seed
%   reproduces a given design within one MATLAB installation.  Where the
%   toolbox is absent, RGASP_LHS_DET provides a deterministic maximin LHS
%   instead; that generator is also what the test suite calls directly, so
%   that recorded test results are identical across MATLAB and Octave.
%
%   See also RGASP_LHS_DET, LHSDESIGN.

if nargin < 3, seed = []; end
if nargin < 4 || isempty(num_try), num_try = 20; end

if ~rgasp_toolbox('lhsdesign')
    if isempty(seed), seed = 0; end
    X = rgasp_lhs_det(n, p, seed, num_try);
    return
end

if ~isempty(seed)
    try
        rng(seed);                                  % MATLAB
    catch
        rand('state', seed); randn('state', seed);  %#ok<RAND>
    end
end
X = lhsdesign(n, p, 'Criterion', 'maximin', 'Iterations', num_try);
end
