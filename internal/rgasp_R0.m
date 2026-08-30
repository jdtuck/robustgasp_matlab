function R0 = rgasp_R0(input1, input2, isotropic)
%RGASP_R0  Absolute-distance matrices, one per input dimension.
%
%   R0 = RGASP_R0(input1)                 distances within a single design
%   R0 = RGASP_R0(input1, input2)         cross distances (rows = input1)
%   R0 = RGASP_R0(input1, input2, true)   single Euclidean-distance matrix
%
%   input1 is n1-by-p, input2 is n2-by-p.  For the anisotropic case R0 is a
%   1-by-p cell array with R0{l}(i,j) = |input1(i,l) - input2(j,l)|.  For the
%   isotropic case R0 is a 1-by-1 cell array holding the Euclidean distance.

if nargin < 2 || isempty(input2)
    input2 = input1;
end
if nargin < 3 || isempty(isotropic)
    isotropic = false;
end

p = size(input1, 2);
if isotropic
    if rgasp_toolbox('pdist2')
        R0 = {pdist2(input1, input2)};            % Statistics Toolbox
    else
        d2 = zeros(size(input1,1), size(input2,1));
        for l = 1:p
            d2 = d2 + bsxfun(@minus, input1(:,l), input2(:,l).').^2;
        end
        R0 = {sqrt(d2)};
    end
else
    R0 = cell(1, p);
    for l = 1:p
        R0{l} = abs(bsxfun(@minus, input1(:,l), input2(:,l).'));
    end
end
end
