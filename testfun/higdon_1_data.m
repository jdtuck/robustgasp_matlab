function y = higdon_1_data(s)
%HIGDON_1_DATA  1-d test function of Higdon (2002), s in [0, 10].
%   y = sin(2*pi*s/10) + 0.2*sin(2*pi*s/2.5)
s = s(:);
y = sin(2*pi*s/10) + 0.2*sin(2*pi*s/2.5);
end
