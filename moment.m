function [m] = moment(x)
l = x(1);
r = x(2);
n = x(3);
global sigma_w rho_w rho_core mu_r a_w V

N_d = 4*(log(l/r)-1)/((l/r)^2-4*log(l/r));    %Demagnetization factor
K = 1+(mu_r-1)/(1+N_d*(mu_r-1));
R = n*sigma_w*(2*pi*r)/a_w; %Total resistance

m = -K*n*V*pi*r^2/R;

end