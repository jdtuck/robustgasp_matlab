function y = borehole(x)
%BOREHOLE  8-d borehole flow-rate model (Worley, 1987).
%
%   y = BOREHOLE(x) with x an n-by-8 matrix of columns
%       rw  radius of borehole (m)        [0.05, 0.15]
%       r   radius of influence (m)       [100, 50000]
%       Tu  transmissivity, upper aquifer [63070, 115600]
%       Hu  potentiometric head, upper    [990, 1110]
%       Tl  transmissivity, lower aquifer [63.1, 116]
%       Hl  potentiometric head, lower    [700, 820]
%       L   length of borehole (m)        [1120, 1680]
%       Kw  hydraulic conductivity        [9855, 12045]
%
%   Returns the water flow rate in m^3/yr.  Inputs r and Tl are nearly inert.
x  = reshape(x, [], 8);
rw = x(:,1); r = x(:,2); Tu = x(:,3); Hu = x(:,4);
Tl = x(:,5); Hl = x(:,6); L = x(:,7); Kw = x(:,8);
res1 = 2*pi*Tu.*(Hu - Hl);
res2 = log(r./rw);
res3 = 1 + (2*L.*Tu)./(res2.*rw.^2.*Kw) + Tu./Tl;
y = res1./res2./res3;
end
