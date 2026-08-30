function model = ppgasp(design, response, varargin)
%PPGASP  Parallel partial Gaussian stochastic process emulator (vector output).
%
%   model = PPGASP(design, response)
%   model = PPGASP(design, response, 'Name', Value, ...)
%
%   For simulators that return a vector (a spatial field, a time series, ...)
%   at every run.  response is n-by-k, one row per run, one column per output
%   coordinate.  All k coordinates share the same correlation function over
%   the input space; the mean vector and the variance of each coordinate are
%   different and are marginalized out separately with the objective prior
%   pi(theta_i, sigma_i^2) ~ 1/sigma_i^2.  The resulting marginal likelihood is
%
%     L ~ |R|^{-k/2} |X' R^{-1} X|^{-k/2} prod_i (S_i^2)^{-(n-q)/2},
%
%   which is maximized (times the jointly robust prior) over log(beta).
%
%   All name/value options of RGASP apply.
%
%   See also PPGASP_PREDICT, RGASP.

if size(response,1) ~= size(design,1)
    error('rgasp:dim', 'response must have one row per design point.');
end
model = gasp_fit(design, response, varargin{:});
end
