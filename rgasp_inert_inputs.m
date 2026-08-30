function P_hat = rgasp_inert_inputs(model, threshold, verbose)
%RGASP_INERT_INPUTS  Detect inert (inactive) inputs from the fitted ranges.
%
%   P_hat = RGASP_INERT_INPUTS(model)
%   P_hat = RGASP_INERT_INPUTS(model, threshold, verbose)
%
%   Computes the normalized inverse range parameters
%
%       P_l = p * C_l beta_l / sum_i (C_i beta_i),      l = 1..p
%
%   which average to 1.  Following the RobustGaSP package, an input with
%   P_l < threshold (default 0.1) is flagged as likely inert.  The jointly
%   robust prior is what makes this diagnostic work: it shrinks the inverse
%   range of a genuinely inactive input towards zero instead of leaving it
%   at an arbitrary likelihood-flat value.
%
%   See also RGASP.

if nargin < 2 || isempty(threshold), threshold = 0.1; end
if nargin < 3 || isempty(verbose), verbose = true; end

w = model.CL(:) .* model.beta_hat(:);
P_hat = model.p_eff * w / sum(w);

if verbose
    fprintf('Estimated normalized inverse range parameters:\n');
    fprintf('  %8.4f', P_hat); fprintf('\n');
    idx = find(P_hat < threshold);
    if isempty(idx)
        fprintf('No input is suspected to be inert.\n');
    else
        fprintf('Inputs suspected to be inert:');
        fprintf(' %d', idx); fprintf('\n');
    end
end
end
