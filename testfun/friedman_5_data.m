function y = friedman_5_data(x)
%FRIEDMAN_5_DATA  5-d test function of Friedman (1991), x in [0,1]^5.
%   Inputs 1-5 are active but x3 enters quadratically and x5 weakly, which
%   makes this a convenient example for the inert-input diagnostic.
x = reshape(x, [], 5);
y = 10*sin(pi*x(:,1).*x(:,2)) + 20*(x(:,3) - 0.5).^2 + 10*x(:,4) + 5*x(:,5);
end
