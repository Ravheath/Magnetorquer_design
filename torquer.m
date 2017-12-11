clear all;
d=0.00031;
L=0.07;
n_layer=2;
D=0.009;
mu_r=4.5;
ro=8900;
sigma=1.7e-8;
V=3.3;
a_w=pi*d*d/4;
N_d=4*(log(2*L/D)-1)/((2*L/D)^2-4*log(2*L/D));
K=1+(mu_r-1)/(1+N_d*(mu_r-1));
m=D*(V*a_w*K)/(4*sigma);


Power=V^2*a_w*d/(sigma*pi*n_layer*L*D);
Mass=ro*L*pi*pi*n_layer*D*d/4;
%ezplot('(3.3^3)*pi*d^5
