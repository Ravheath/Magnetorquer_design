clear;
clc;
d = 0.079e-3;%1.5*0.00031; %diameter of wire
L = 0.09;%0.15; %Length of rod
n_layer = 6;  %number of layers
D = 0.04;  %Diameter of rod
mu_r = 1;% 75e3; %Relative ermeability of rod  
rho = 8900;  %mass density of wire
sigma = 1.68e-8;%1.7e-8; %resistivity of wire
V = 3.3;    %voltage
a_w = 0.00501e-6;%0.25*pi*d^2;  %cross sectional area of wire

N_d = 4*(log(2*L/D)-1)/((2*L/D)^2-4*log(2*L/D));    %Demagnetization factor
K = 1+(mu_r-1)/(1+N_d*(mu_r-1));
m = D*(V*a_w*K)/(4*sigma);  %magnetic moment

Power = V^2*a_w*d/(sigma*pi*n_layer*L*D); %power dissipated
Mass = rho*L*pi*pi*n_layer*D*d/4;
%ezplot('(3.3^3)*pi*d^5
