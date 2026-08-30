function [x, fval, info] = rgasp_neldermead(fun, x0, lb, ub, maxeval)
%RGASP_NELDERMEAD  Derivative-free simplex search with box projection.
%
%   Provided as a fallback optimizer (RobustGaSP offers 'nelder-mead' too).
%   The objective is evaluated at projected points so the bounds are honoured.

x0 = x0(:); nx = numel(x0);
if isempty(lb), lb = -inf(nx,1); end
if isempty(ub), ub =  inf(nx,1); end
proj = @(z) min(max(z, lb(:)), ub(:));

f = @(z) fun_only(fun, proj(z));

% initial simplex
V = repmat(x0, 1, nx+1);
for i = 1:nx
    if x0(i) ~= 0
        V(i, i+1) = 1.05 * x0(i);
    else
        V(i, i+1) = 0.00025;
    end
end
F = zeros(1, nx+1);
for i = 1:nx+1, F(i) = f(V(:,i)); end
nfev = nx + 1;

rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
while nfev < maxeval
    [F, idx] = sort(F); V = V(:, idx);
    if max(abs(F(2:end) - F(1))) < 1e-10 * max(1, abs(F(1))) && ...
       max(max(abs(V(:,2:end) - repmat(V(:,1),1,nx)))) < 1e-8
        break
    end
    xbar = mean(V(:,1:nx), 2);
    xr = xbar + rho*(xbar - V(:,end)); fr = f(xr); nfev = nfev + 1;
    if fr < F(1)
        xe = xbar + chi*(xr - xbar); fe = f(xe); nfev = nfev + 1;
        if fe < fr, V(:,end) = xe; F(end) = fe; else, V(:,end) = xr; F(end) = fr; end
    elseif fr < F(nx)
        V(:,end) = xr; F(end) = fr;
    else
        if fr < F(end)
            xc = xbar + psi*(xr - xbar); fc = f(xc); nfev = nfev + 1;
            if fc <= fr, V(:,end) = xc; F(end) = fc; else, [V,F,nfev] = shrink(V,F,f,sigma,nfev,nx); end
        else
            xcc = xbar - psi*(xbar - V(:,end)); fcc = f(xcc); nfev = nfev + 1;
            if fcc < F(end), V(:,end) = xcc; F(end) = fcc; else, [V,F,nfev] = shrink(V,F,f,sigma,nfev,nx); end
        end
    end
end

[fval, i] = min(F);
x = proj(V(:,i));
info = struct('iterations', nfev, 'funcCount', nfev, ...
              'message', 'nelder-mead terminated', 'converged', true);
end

function [V,F,nfev] = shrink(V,F,f,sigma,nfev,nx)
for i = 2:nx+1
    V(:,i) = V(:,1) + sigma*(V(:,i) - V(:,1));
    F(i) = f(V(:,i)); nfev = nfev + 1;
end
end

function v = fun_only(fun, z)
v = fun(z);
end
