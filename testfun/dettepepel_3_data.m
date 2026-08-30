function y = dettepepel_3_data(x)
%DETTEPEPEL_3_DATA  3-d curved test function of Dette and Pepelyshev (2010).
%   x in [0,1]^3, n-by-3 input, n-by-1 output.
x = reshape(x, [], 3);
x1 = x(:,1); x2 = x(:,2); x3 = x(:,3);
y = 4*(x1 - 2 + 8*x2 - 8*x2.^2).^2 + (3 - 4*x2).^2 ...
    + 16*sqrt(x3 + 1).*(2*x3 - 1).^2;
end
