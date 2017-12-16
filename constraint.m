function [c,ceq] = constraint(x)
global sigma_w rho_w rho_core mu_r a_w V
l = x(1);
r = x(2);
n = x(3);

R = n*sigma_w*(2*pi*r)/a_w; %Total resistance
c = [rho_core*pi*r^2*l + a_w*(2*pi*r*n)*rho_w - 0.1 , V^2/R - 0.2];
ceq = [];
end
