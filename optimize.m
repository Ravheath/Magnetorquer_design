clc
global sigma_w rho_w rho_core mu_r a_w V
mu_r = 2000; %relative permeability
V = 5; %Voltage
sigma_w = 1.55e-8;    %resistivity
a_w = 7.97e-9;   %gauge of wire
rho_w = 8.93e3;
rho_core = 8.74e3;

x0 = [5e-2,5e-3,1000];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0, 3e-3, 0];
ub = [0.15,5e-3,10e3];
options = optimoptions('fmincon','TolX',1e-10);
x = fmincon(@moment,x0,A,b,Aeq,beq,lb,ub,@constraint,options);

Mass = rho_core*pi*x(2)^2*x(1) + a_w*(2*pi*x(2)*x(3))*rho_w;
power = V^2/(x(3)*sigma_w*(2*pi*x(2))/a_w);
m = moment(x);