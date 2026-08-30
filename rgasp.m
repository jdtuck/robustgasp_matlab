function model = rgasp(design, response, varargin)
%RGASP  Robust Gaussian stochastic process emulator (scalar output).
%
%   model = RGASP(design, response)
%   model = RGASP(design, response, 'Name', Value, ...)
%
%   Fits the Gaussian stochastic process (GaSP) emulator
%
%       y(x) = h(x) theta + z(x),   z ~ GaSP(0, sigma^2 c(.,.))
%
%   with a separable (product) correlation function, estimating the range
%   parameters gamma_1..gamma_p by the mode of the marginal posterior in the
%   log inverse-range parameterization xi_l = log(1/gamma_l), using the
%   jointly robust prior.  This is the "robust" estimator studied in
%
%     Gu, M., Wang, X. and Berger, J.O. (2018).  Robust Gaussian stochastic
%     process emulation.  The Annals of Statistics 46(6A), 3038-3066.
%     (arXiv:1708.04738)
%
%   and implemented in the R package RobustGaSP.
%
%   INPUTS
%     design    n-by-p matrix of design points (one computer-model run per row)
%     response  n-by-1 vector of scalar outputs
%
%   NAME/VALUE OPTIONS (defaults match RobustGaSP)
%     'trend'            n-by-q basis matrix h(x^D).  Default ones(n,1).
%     'zeroMean'         true to force mu(x) = 0.  Default false.
%     'nugget'           fixed nugget-variance ratio eta.  Default 0.
%     'nuggetEst'        true to estimate eta.  Default false.
%     'rangePar'         fix gamma at these values (skips optimization).
%     'method'           'post_mode' (default, marginal posterior mode with the
%                        jointly robust prior), 'mmle' (marginal likelihood,
%                        no prior) or 'mle' (profile likelihood).
%     'a','b'            jointly robust prior hyper-parameters.  Defaults
%                        a = 0.2 and b = (a+p)/n^(1/p).
%     'kernelType'       'matern_5_2' (default), 'matern_3_2' or 'pow_exp'.
%                        May be a 1-by-p cell array for mixed kernels.
%     'alpha'            roughness for 'pow_exp'.  Default 1.9.
%     'isotropic'        true for a single range parameter.  Default false.
%     'lowerBound'       constrain log(beta) from below.  Default true.
%     'numInitialValues' number of optimization restarts.  Default 2.
%     'maxEval'          optimizer budget.  Default max(30, 20+5p).
%     'optimizer'        'lbfgs' (default), 'neldermead' or 'fmincon'.
%
%   OUTPUT
%     model  struct with (among others)
%       .beta_hat    estimated inverse range parameters
%       .range_hat   estimated range parameters gamma = 1/beta
%       .nugget      nugget-variance ratio used/estimated
%       .theta_hat   generalized least squares trend coefficients
%       .sigma2_hat  estimated variance
%       .log_post    value of the maximized log marginal posterior
%       .L, .LX      Cholesky factors reused by RGASP_PREDICT
%
%   EXAMPLE
%     x  = linspace(0,10,15)';
%     y  = higdon_1_data(x);
%     m  = rgasp(x, y);
%     xt = linspace(0,10,200)';
%     pr = rgasp_predict(m, xt);
%
%   See also RGASP_PREDICT, RGASP_SIMULATE, RGASP_LOO, RGASP_INERT_INPUTS,
%            PPGASP.

if size(response, 2) ~= 1 && numel(response) ~= size(design,1)
    error('rgasp:dim', 'rgasp expects a scalar output; use ppgasp for vector output.');
end
response = response(:);
model = gasp_fit(design, response, varargin{:});
end
