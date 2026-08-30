function ws = rgasp_mute_warnings()
%RGASP_MUTE_WARNINGS  Temporarily silence linear-algebra conditioning warnings.
%
%   ws = RGASP_MUTE_WARNINGS();  ... ;  warning(ws);
%
%   During optimization the correlation matrix is deliberately evaluated at
%   parameter values where it is nearly singular (the marginal likelihood is
%   flat as beta -> 0); those states are detected and penalized inside
%   RGASP_OBJECTIVE, so the accompanying solver warnings are noise.

ws = warning();
ids = {'MATLAB:nearlySingularMatrix', 'MATLAB:singularMatrix', ...
       'MATLAB:illConditionedMatrix', 'MATLAB:rankDeficientMatrix', ...
       'Octave:singular-matrix', 'Octave:nearly-singular-matrix', ...
       'Octave:singular-matrix-div'};
for i = 1:numel(ids)
    try %#ok<TRYNC>
        warning('off', ids{i});
    end
end
end
