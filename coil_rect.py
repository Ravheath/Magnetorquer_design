from numpy import array, zeros, pi
from math import sqrt, log
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import scipy.io as sio

def Inductance_rect(l,b,n,d):	#Used model from https://www.eeweb.com/tools/rectangle-loop-inductance
	f = sqrt(l**2+b**2)
	temp = -b*log((b+f)/l)-l*log((l+f)/b)+b*log(4*b/d)+l*log(4*l/d)
	I = 4e-7*(-2*(l+b)+2*f+temp)*n**2
	return I


#This code generates array x with columns: length width number of turns mass power moment current (all in SI)
#Array is sorted in decreasing order of moment
#Parameters for Pratham torquer from Integration CDR:
# Copper wire diameter = 0.28mm
# Resistivity = 1.7e-8
# Mass density = 8900 kg/m3

t1 = time.time()
#Wire parameter 
V = 3.3	#voltage applied to torquer
mu_r = 1	#relative permiability of core (air in this case)
sigma_w = 1.7e-8	#Resistivity of wire (in SI units)
rho_w = 8.93e3
diameter_w = 0.143e-3#0.28e-3	#cross sectional diameter of wire
a_w = pi*0.25*diameter_w**2	#gauge of wire (cross sectional area in m^2)

l_i = 4.0e-2	#min length (in m)
l_f = 18.0e-2	#max length	(in m)
b_i = 4.0e-2
b_f = 8.0e-2
w_i = 1.0e-3	#min width	(in m)
w_f = 1.0e-2	#max width	(in m)
dl = 1e-3	#incremental step in length
dw = 1e-3	#Incremental step in width
db = 1e-3	#Incremental step in breadth

size = int((l_f-l_i+dl)*(w_f-w_i+dw)*(b_f-b_i+db)/(dl*dw*db))	#size of array
#[length;  width;	number of turns;  Mass;  power;	moment;	current;	 inductance;  resistance; time constant  ]
x = zeros((size,11),float)
k = 0
y = []

for p in range (1,int((l_f-l_i)/dl+2)):
	for q in range (1,int((b_f-b_i)/db+2)):
		for r in range (1,int((w_f-w_i)/dw+2)):
			x[k,0] = l_i + dl*(float(p)-1)	#length	(in m)
			x[k,1] = b_i + db*(float(q)-1)	#breadth (in m)
			x[k,2] = w_i + dw*(float(r)-1)	#width (in m)
			x[k,3] = int(x[k,2]/diameter_w)	#number of turns
			x[k,4] = 2*rho_w*a_w*x[k,3]*(x[k,0]+x[k,1])	#Mass
			x[k,9] = 2*sigma_w*x[k,3]*(x[k,0]+x[k,1])/a_w + 20 #Resistance
			x[k,5] = V**2/x[k,9]	#Power
			x[k,6] = x[k,3]*V*x[k,0]*x[k,1]/x[k,9]	#moment
			x[k,7] = V/x[k,9]	#current
			x[k,8] = Inductance_rect(x[k,0],x[k,1],x[k,3],diameter_w)	#inductance
			
			x[k,10] = x[k,8]/x[k,9]	#Time constant
			
			if x[k,4] < 30e-3 and  x[k,6] > 0.01:
				y.append(list(x[k,:]))
			k = k+1

print x.shape[0] , size, k

x = x[x[:,5].argsort()[::1]]	#sort array in ascending order w.r.t. power 
z = array(y)
z = z[z[:,5].argsort()[::1]]

t2 = time.time()
print t2-t1

sio.savemat('moment.mat', mdict={'x':x})
sio.savemat('allowed_2U_rect.mat', mdict={'z':z})