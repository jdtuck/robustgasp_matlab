function ok = run_all_tests(mode)
%RUN_ALL_TESTS  Run the RobustGaSP-MATLAB test suite.
%
%   ok = RUN_ALL_TESTS()             run as installed (toolboxes if present)
%   ok = RUN_ALL_TESTS('fallback')   force every built-in fallback, so the
%                                    suite also covers a machine with no
%                                    Statistics or Optimization Toolbox
%
%   Returns true if every test passes.

if nargin < 1, mode = 'default'; end

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_robustgasp.m'));

if strcmpi(mode, 'fallback')
    rgasp_use_toolboxes(false);
    fprintf('Running with the toolboxes DISABLED (built-in fallbacks only).\n');
else
    rgasp_use_toolboxes(true);
end
restore = onCleanup(@() rgasp_use_toolboxes(true));

tests = {@test_kernels, @test_gradients, @test_optimizer, ...
         @test_toolbox_paths, @test_mex_equivalence, ...
         @test_rgasp_fit_predict, @test_ppgasp, ...
         @test_loo_and_utils, @test_paper_behaviour};

npass = 0; nfail = 0; failures = {};
for i = 1:numel(tests)
    name = func2str(tests{i});
    fprintf('\n=== %s ===\n', name);
    try
        tests{i}();
        fprintf('--- %s PASSED\n', name);
        npass = npass + 1;
    catch err
        fprintf(2, '--- %s FAILED: %s\n', name, err.message);
        if ~isempty(err.stack)
            fprintf(2, '    at %s line %d\n', err.stack(1).name, err.stack(1).line);
        end
        nfail = nfail + 1;
        failures{end+1} = name; %#ok<AGROW>
    end
end

fprintf('\n============================================\n');
fprintf('%d passed, %d failed\n', npass, nfail);
for i = 1:numel(failures)
    fprintf('  FAILED: %s\n', failures{i});
end
fprintf('============================================\n');
ok = (nfail == 0);
end
