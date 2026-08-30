function pred = ppgasp_predict(model, testing_input, varargin)
%PPGASP_PREDICT  Predictive distribution for a PP GaSP emulator.
%
%   pred = PPGASP_PREDICT(model, testing_input, ...)
%
%   Fields of pred are nt-by-k matrices: each column corresponds to one
%   output coordinate and carries its own sigma_i^2, while the correlation
%   term c**(x*) is shared across coordinates.  The predictive distribution
%   of each coordinate is Student t with n - q degrees of freedom.
%
%   See also PPGASP, RGASP_PREDICT.

pred = gasp_predict(model, testing_input, varargin{:});
end
