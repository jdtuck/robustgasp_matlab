function test_rgasp_model_object()
%TEST_RGASP_MODEL_OBJECT  The rgasp_model classdef wrapper and its predict method.
%
%   rgasp() returns an rgasp_model object whose .model field carries the
%   fitted struct used throughout the rest of the suite.  This test covers
%   the object API itself: the constructor, the stored struct, and the
%   Monte-Carlo predict method that draws posterior sample paths.

% ---- 1. rgasp returns an rgasp_model wrapping the fitted struct --------
x = linspace(0, 10, 15)';
y = higdon_1_data(x);
obj = rgasp(x, y);

if ~isa(obj, 'rgasp_model')
    error('rgasp must return an rgasp_model object, got %s', class(obj));
end
if ~isstruct(obj.model)
    error('rgasp_model.model must be the fitted struct');
end
if ~strcmp(obj.model.type, 'rgasp')
    error('wrapped model has the wrong type: %s', obj.model.type);
end

% the constructor stores the struct verbatim and leaves samples empty
m2 = rgasp_model(obj.model);
if ~isequaln(m2.model, obj.model)
    error('rgasp_model constructor did not store the model unchanged');
end
if ~isempty(m2.samples)
    error('a freshly constructed rgasp_model must have empty samples');
end

% ---- 2. predict draws sample paths of the requested shape --------------
xt = linspace(0, 10, 12)';
B  = 400;
P  = obj.predict(xt, 'B', B, 'idxSamples', 1:B);
if any(size(P) ~= [B numel(xt)])
    error('predict returned shape %s, expected [%d %d]', ...
          mat2str(size(P)), B, numel(xt));
end

% ---- 3. the sample paths are centred on the predictive mean -----------
pr = rgasp_predict(obj.model, xt);
mc_mean = mean(P, 1)';                       % Monte-Carlo mean per test point
tol = 5*max(pr.sd)/sqrt(B) + 1e-8;
if max(abs(mc_mean - pr.mean)) > tol
    error(['predict sample paths are not centred on the predictive mean ' ...
           '(max deviation %.3g > %.3g)'], max(abs(mc_mean - pr.mean)), tol);
end

% at the design points the process interpolates, so the draws must collapse
Pd = obj.predict(x, 'B', B);
if max(std(Pd, 0, 1)) > 1e-4
    error('sample paths do not collapse at the (noise-free) design points');
end
assert_close(mean(Pd, 1)', y, 1e-4, 'predict mean at design points');

% ---- 4. idxSamples selects a subset of the drawn paths -----------------
sel = [1 5 9];
Ps  = obj.predict(xt, 'B', B, 'idxSamples', sel);
if any(size(Ps) ~= [numel(sel) numel(xt)])
    error('idxSamples subset returned the wrong shape: %s', mat2str(size(Ps)));
end

% ---- 5. defaults: no name/value args at all ---------------------------
Pdef = obj.predict(xt);                      % B defaults to 500
if any(size(Pdef) ~= [500 numel(xt)])
    error('predict default B should be 500, got %d rows', size(Pdef,1));
end

fprintf('rgasp_model: object wrapper, constructor and Monte-Carlo predict OK\n');
end
