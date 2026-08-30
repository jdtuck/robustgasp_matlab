function ok = rgasp_compile_mex(verbose)
%RGASP_COMPILE_MEX  Build the optional C++ acceleration for RobustGaSP-MATLAB.
%
%   rgasp_compile_mex()          build (quiet on success)
%   ok = rgasp_compile_mex(true) build with compiler output
%
%   The toolbox is fully functional without this step -- every MEX function
%   has a pure-MATLAB fallback that produces identical results.  Compiling is
%   worthwhile once n reaches a few hundred and matters a lot in the
%   thousands, because it removes the two dominant costs of a fit:
%
%     * building the separable correlation matrix (p exponentials per entry),
%     * the per-dimension gradient terms,
%
%   and, just as importantly, it stops the p per-dimension distance matrices
%   from ever being allocated (p*n^2 doubles).
%
%   You need a C++ compiler configured for MEX ("mex -setup C++" in MATLAB).
%   OpenMP is used when the compiler accepts it, and silently skipped
%   otherwise.  Run this once; the built files live in internal/mex.
%
%   See also RGASP_HAVE_MEX, RGASP.

if nargin < 1 || isempty(verbose), verbose = false; end

here = fileparts(mfilename('fullpath'));
src  = fullfile(here, 'internal', 'mex');
files = {'rgasp_corr_mex.cpp', 'rgasp_gradterms_mex.cpp'};

isOctave = exist('OCTAVE_VERSION', 'builtin') == 5;

% Flag sets tried in order: OpenMP first, then a plain build.
if isOctave
    flagsets = {{'-fopenmp'}, {}};
else
    flagsets = {{'openmp'}, {}};
end

addpath(src);
old = pwd;
cd(src);
restore = onCleanup(@() cd(old));

ok = false;
for fs = 1:numel(flagsets)
    try
        for i = 1:numel(files)
            build_one(files{i}, flagsets{fs}, isOctave, verbose);
        end
        ok = true;
        if fs == 1
            fprintf('RobustGaSP-MATLAB: MEX acceleration built (with OpenMP).\n');
        else
            fprintf(['RobustGaSP-MATLAB: MEX acceleration built ' ...
                     '(no OpenMP; single threaded).\n']);
        end
        break
    catch err
        if fs == numel(flagsets)
            fprintf(2, 'RobustGaSP-MATLAB: MEX build failed (%s).\n', err.message);
            fprintf(2, ['The toolbox still works -- it will use the pure ' ...
                        'MATLAB code paths.\n']);
        end
    end
end

rgasp_have_mex(true);   % refresh the cached availability flag
end

% ----------------------------------------------------------------------
function build_one(file, flags, isOctave, verbose)
if isOctave
    oldC = getenv('CXXFLAGS'); oldL = getenv('LDFLAGS');
    c = '-O3'; l = '';
    if ~isempty(flags), c = [c ' ' flags{1}]; l = flags{1}; end
    setenv('CXXFLAGS', c);
    setenv('LDFLAGS', l);
    cleanup = onCleanup(@() restore_env(oldC, oldL));
    if verbose
        mex(file);
    else
        evalc('mex(file)');
    end
else
    args = {'-O', file};
    if ~isempty(flags) && strcmp(flags{1}, 'openmp')
        if ispc
            args = [{'-O', 'COMPFLAGS=$COMPFLAGS /openmp'}, {file}];
        else
            args = [{'-O', 'CXXFLAGS=$CXXFLAGS -fopenmp', ...
                     'LDFLAGS=$LDFLAGS -fopenmp'}, {file}];
        end
    end
    if verbose
        mex(args{:});
    else
        evalc('mex(args{:})');
    end
end
end

function restore_env(c, l)
setenv('CXXFLAGS', c);
setenv('LDFLAGS', l);
end
