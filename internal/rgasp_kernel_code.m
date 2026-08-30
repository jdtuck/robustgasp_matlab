function code = rgasp_kernel_code(kernel_type)
%RGASP_KERNEL_CODE  Numeric kernel codes used by the MEX functions.
%
%   1 = pow_exp, 2 = matern_3_2, 3 = matern_5_2

if ischar(kernel_type), kernel_type = {kernel_type}; end
code = zeros(numel(kernel_type), 1);
for i = 1:numel(kernel_type)
    switch lower(kernel_type{i})
        case 'pow_exp',    code(i) = 1;
        case 'matern_3_2', code(i) = 2;
        case 'matern_5_2', code(i) = 3;
        otherwise
            error('rgasp:kernel', 'Unknown kernel_type ''%s''.', kernel_type{i});
    end
end
end
