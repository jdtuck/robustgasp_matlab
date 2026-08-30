function y = limetal_2_data(x)
%LIMETAL_2_DATA  2-d test function of Lim et al. (2002), x in [0,1]^2.
%   Accepts an n-by-2 matrix and returns an n-by-1 vector.
x = reshape(x, [], 2);
res1 = 30 + 5*x(:,1).*sin(5*x(:,1));
res2 = 4 + exp(-5*x(:,2));
y = (res1.*res2 - 100)/6;
end
