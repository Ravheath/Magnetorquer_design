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
l_f = 8.0e-2	#max length	(in m)
w_i = 1e-3	#min width	(in m)
w_f = 1.0e-2	#max width	(in m)
dl = 1e-3	#incremental step in length
dw = 1e-3	#Incremental step in width

size = int((l_f-l_i+dl)*(w_f-w_i+dw)/(dl*dw))	#size of array
#[length;  width;	number of turns;  Mass;  power;	moment;	current;	 inductance;  resistance; time constant  ]
x = zeros((size,10),float)
k = 0
y = []
for p in range (1,int((l_f-l_i)/dl+2)):
	for q in range (1,int((w_f-w_i)/dw+2)):
		x[k,0] = l_i + dl*(float(p)-1)	#length	(in m)
		x[k,1] = w_i + dw*(float(q)-1)	#width	(in m)
		x[k,2] = int(x[k,1]/diameter_w) #Number of turns = width / diameter of wire 
		
		x[k,8] = (4*sigma_w/a_w)*x[k,0]*x[k,2]   #Resistance
		x[k,3] = (4*rho_w*a_w/diameter_w)*x[k,0]*x[k,1]	#mass = (4*rho_w*a_w/d_w)*l*b	(in kg)
		x[k,4] = V**2/(x[k,8])	#power = (0.25*V^2*a_w*d_w/sigma_w)/(lb)	(in W)
		x[k,5] = x[k,2]*V*x[k,0]**2/x[k,8] #(0.25*a_w*V/sigma_w)*x[k,0]	#mu = (0.25*a_w*V/sigma_w)*l 	(in Am^2)
		x[k,6] = V/x[k,8] #current
		x[k,7] = Inductance(x[k,0],x[k,2],diameter_w) 	#Inductance
		
		x[k,9] = x[k,7]/x[k,8]
		if x[k,3] < 30e-3 and  x[k,5] > 0.01:
			y.append(list(x[k,:]))
			
		k = k+1
		print p,q
print x.shape[0] , size, k

x = x[x[:,4].argsort()[::1]]	#sort array in ascending order w.r.t. power 
z = array(y)
z = z[z[:,4].argsort()[::1]]

t2 = time.time()
print t2-t1

sio.savemat('moment.mat', mdict={'x':x})
sio.savemat('allowed_1U.mat', mdict={'z':z})


