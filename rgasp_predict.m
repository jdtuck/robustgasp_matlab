function pred = rgasp_predict(model, testing_input, varargin)
%RGASP_PREDICT  Predictive mean, sd and 95% interval for an RGASP emulator.
%
%   pred = RGASP_PREDICT(model, testing_input)
%   pred = RGASP_PREDICT(model, testing_input, 'Name', Value, ...)
%
%   OPTIONS
%     'trend'         nt-by-q basis matrix at the testing inputs.  Required
%                     when the model was fitted with a non-constant trend.
%     'intervalData'  true (default) gives intervals for a noisy observation
%                     (adds the nugget); false gives intervals for the mean.
%     'level'         credible level.  Default 0.95.
%
%   OUTPUT struct with fields mean, sd, lower95, upper95, df.
%
%   See also RGASP, RGASP_SIMULATE, PPGASP_PREDICT.

if ~strcmp(model.type, 'rgasp')
    error('rgasp:type', 'Use ppgasp_predict for a ppgasp model.');
end
pred = gasp_predict(model, testing_input, varargin{:});
end
