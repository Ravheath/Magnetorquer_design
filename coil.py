from numpy import array, zeros, pi
from math import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import scipy.io as sio
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
diameter_w = 0.5*0.28e-3	#cross sectional diameter of wire
a_w = pi*0.25*diameter_w**2	#gauge of wire (corss sectional area in m^2)

l_i = 3e-2	#min length (in m)
l_f = 9e-2	#max length	(in m)
w_i = 0.5e-3	#min width	(in m)
w_f = 3.0e-2	#max width	(in m)
dl = 1e-3	#incremental step in length
dw = 1e-3	#Incremental step in width

size = int((l_f-l_i+dl)*(w_f-w_i+dw)/(dl*dw))	#size of array

x = zeros((size,7),float)
k = 0

for p in range (1,int((l_f-l_i)/dl+2)):
	for q in range (1,int((w_f-w_i)/dw+2)):
		x[k,0] = l_i + dl*(float(p)-1)	#length	(in m)
		x[k,1] = w_i + dw*(float(q)-1)	#width	(in m)
		x[k,2] = int(x[k,1]/diameter_w) #Number of turns = width / diameter of wire 
		x[k,3] = (4*rho_w*a_w/diameter_w)*x[k,0]*x[k,1]	#mass = (4*rho_w*a_w/d_w)*l*b	(in kg)
		x[k,4] = (0.25*(V**2)*a_w*diameter_w/sigma_w)/(x[k,0]*x[k,1])	#power = (0.25*V^2*a_w*d_w/sigma_w)/(lb)	(in W)
		x[k,5] =  (0.25*a_w*V/sigma_w)*x[k,0]	#mu = (0.25*a_w*V/sigma_w)*l 	(in Am^2)
		x[k,6] = (0.25*V*a_w*diameter_w/sigma_w)/(x[k,0]*x[k,1]) #current
		
		k = k+1
		
print x.shape[0] , size, k

x = x[x[:,5].argsort()[::-1]]	#sort array in descending order w.r.t. moment 

t2 = time.time()
print t2-t1
sio.savemat('moment.mat', mdict={'x':x})
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0]*100,x[:,1],x[:,5])
ax.set_xlabel('length (cm)')
ax.set_ylabel('turns ')
ax.set_zlabel('moment (g)')
plt.show()

plt.plot(x[:,1],x[:,4])
plt.show()

y = zeros((j,6),float)
j = 0
for i in range(size):
	if x[i,2] < 0.02 and x[i,3]<20e-3 and x[i,4]<0.2 and x[i,5]>0.008:
		y[j,:] = x[i,:]
		j = j+1
print j

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(y[:,0]*100,y[:,1],y[:,2]*1e2)
ax.set_xlabel('length (cm)')
ax.set_ylabel('turns ')
ax.set_zlabel('power (mW)')
plt.show()


x = array([[4,5,6],[1,20,3],[7,0.8,9],[10,11,12]])

x = x[x[:,2].argsort()]
'''
