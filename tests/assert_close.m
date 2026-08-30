function assert_close(a, b, tol, msg)
%ASSERT_CLOSE  Assert two arrays agree to within a relative/absolute tolerance.
if nargin < 3 || isempty(tol), tol = 1e-8; end
if nargin < 4, msg = ''; end
a = a(:); b = b(:);
if numel(a) ~= numel(b)
    error('assert:size', 'Size mismatch (%d vs %d). %s', numel(a), numel(b), msg);
end
den = max(1, max(abs(a), abs(b)));
err = max(abs(a-b) ./ den);
if ~(err <= tol)
    error('assert:close', 'Max relative error %.3g exceeds tol %.3g. %s', err, tol, msg);
end
end
