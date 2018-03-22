from numpy import array, zeros, pi
from math import sqrt, log
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import scipy.io as sio

def Inductance(l,n,d):	#Used model from https://www.eeweb.com/tools/rectangle-loop-inductance
	I = 4e-7*l*(2*sqrt(2)-4+2*log(1.65685*l/d))*(n**2)
	return I


#This code generates array x with columns: length width number of turns mass power moment current (all in SI)
#Array is sorted in decreasing order of moment
#Parameters for Pratham torquer from Integration CDR:
# Copper wire diameter = 0.28mm
# Resistivity = 1.7e-8
# Mass density = 8900 kg/m3

t1 = time.time()
#Wire paramters 
V = 3.3	#voltage applied to torquer
mu_r = 1	#relative permiability of core (air in this case)
sigma_w = 1.7e-8	#Resistivity of wire (in SI units)
rho_w = 8.93e3
diameter_w = 0.143e-3#0.28e-3	#cross sectional diameter of wire
a_w = pi*0.25*diameter_w**2	#gauge of wire (cross sectional area in m^2)

l_i = 4e-2	#min length (in m)
l_f = 8.5e-2	#max length	(in m)
n_i = 10	#min turns
n_f = 300	#max turns
dl = 1e-3	#incremental step in length
dn = 1	#Incremental step in turns
size = int((l_f-l_i+dl)*(n_f-n_i+dn)/(dl*dn))	#size of array

#[length;  number of turns;  Resistance;Mass;  power;	moment;	current;	 inductance;   time constant  ]
x = zeros((size,9),float)
k = 0
y = []

for p in range (1,int((l_f-l_i)/dl+2)):
	for q in range (1,int((n_f-n_i)/dn+2)):
		x[k,0] = l_i + dl*(float(p)-1)	#length	(in m)
		x[k,1] = n_i + dn*(float(q)-1)	#number of turns
		x[k,2] = (4*sigma_w/a_w)*x[k,0]*x[k,1]   #Resistance
		x[k,3] = (pi*rho_w*diameter_w**2)*x[k,0]*x[k,1]	#mass = (pi*rho*d^2)*l*n	(in kg)
		x[k,4] = (5/4)*(V**2)/(x[k,2])	#power = 5/4*V^2/R	(in W)
		x[k,5] = V*(x[k,1]/x[k,2])*x[k,0]**2	# V(n/R)l^2(in Am^2)
		x[k,6] = V/x[k,2] #current
		x[k,7] = Inductance(x[k,0],x[k,1],diameter_w) 	#Inductance
		
		x[k,8] = x[k,7]/x[k,2]	#time constant

		if x[k,3] < 30e-3 and  x[k,5] > 0.01 :
			y.append(list(x[k,:]))
			
		k = k+1
		
print x.shape[0] , size, k ,len(y)

x = x[x[:,4].argsort()[::1]]	#sort array in ascending order w.r.t. power 
z = array(y)
z = z[z[:,4].argsort()[::1]]

#print type(z)
#list(z)
#z.tolist()
#print type(z)
#z = [['length','turns','Resistance','Mass','power','moment','current','inductance','time constant']] + z

t2 = time.time()
print t2-t1
print z[1,:]
#sio.savemat('moment.mat', mdict={'x':x})
sio.savemat('sq_l_n_allowed.mat', mdict={'z':z})



