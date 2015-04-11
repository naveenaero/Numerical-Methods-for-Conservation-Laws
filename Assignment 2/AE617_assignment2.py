'''
AE617 : Numerical Methods for Conservation Laws
Assignment 2: System of Linear Equations -----> Acoustic Equations 
'''

'''
Import the necessary libraries
'''
import numpy as np
import pylab as plt
from matplotlib import *


index = int(raw_input('Enter the Question Number for which you want to compute: '));
T = float(raw_input('Enter the time at which you want to see the Result: '));

'''
Define important variables
'''
dx = 0.005
dt = 0.004
dtdx = dt/dx
nx = int(1/dx) + 1
nt = 200

x1 = np.linspace(0,1,nx)
u = np.zeros((nx,nt))
p = np.zeros((nx,nt))
rho0 = np.zeros((nx,nt))
K0 = 1
xo = 0.4
xbar = 0.075
pbar = 0.2


'''
Initialise the Conserved variables
'''
for i in range(nx):
	u[i,0] = 0;
	if x1[i] > xo - xbar and x1[i] < xo + xbar:
		p[i,0] = pbar*(1-((x1[i]-xo)/xbar)**2)**0.5

if index == 1:
	for i in range(nx):
		rho0[i] = 1
elif index == 2:
	for i in range(nx):
		if x1[i] < 0.6:
			rho0[i] = 1
		else:
			rho0[i] = 3


U = np.zeros((2,nx))
Un = np.zeros((2,nx))
U[0,:] = u[:,0].transpose()
U[1,:] = p[:,0].transpose() 

# print np.shape(Un)
# print np.shape(U)

A1 = np.array([[0,1],[1,0]])
A2 = np.array([[0,1],[1.0/3.0,0]])


lamda1= (np.linalg.eig(A1))[0]
lamda2 = (np.linalg.eig(A2))[0]


r1 = (np.linalg.eig(A1))[1]
r2 = (np.linalg.eig(A2))[1]


if lamda1[0] > 0:
	Aplus1 = np.dot(r1,np.dot(np.diag((lamda1[0],0)),np.linalg.inv(r1)))
	Aminus1 = np.dot(r1,np.dot(np.diag((0,lamda1[1])),np.linalg.inv(r1)))
else:
	Aplus1 = np.dot(r1,np.dot(np.diag((lamda1[1],0)),np.linalg.inv(r1)))
	Aminus1 = np.dot(r1,np.dot(np.diag((0,lamda1[0])),np.linalg.inv(r1)))

if lamda2[0] > 0:
	Aplus2 = np.dot(r2,np.dot(np.diag((lamda2[0],0)),np.linalg.inv(r2)))
	Aminus2 = np.dot(r2,np.dot(np.diag((0,lamda2[1])),np.linalg.inv(r2)))
else:
	Aplus2 = np.dot(r2,np.dot(np.diag((lamda2[1],0)),np.linalg.inv(r2)))
	Aminus2 = np.dot(r2,np.dot(np.diag((0,lamda2[0])),np.linalg.inv(r2)))


Z = 1

for t in range(1,nt):
	for x in range(nx):
		# if x == 0:
		# 	# Un[:,x] = U[:,x] - dtdx*(np.dot(Aplus1,U[:,x]-U[:,x-1]) + np.dot(Aminus1,U[:,x+1]-U[:,x]))
		# 	Un[:,x] = U[:,x] - dtdx*(np.dot(Aminus1,U[:,x+1]-U[:,x]))
		# if x == nx-1:
		# 	Un[:,x] = U[:,x] - dtdx*(np.dot(Aplus2,U[:,x]-U[:,x-1]))
		# if x > 0 and x < nx-1:
		# 	if x <= 121:
		# 		Un[:,x] = U[:,x] - dtdx*(np.dot(Aplus1,U[:,x]-U[:,x-1]) + np.dot(Aminus1,U[:,x+1]-U[:,x]))
		# 	else:
		# 		Un[:,x] = U[:,x] - dtdx*(np.dot(Aplus2,U[:,x]-U[:,x-1]) + np.dot(Aminus2,U[:,x+1]-U[:,x]))
		
		u[x,t] = 

	Un[0,0] = 0
	u[:,t] = Un[0,:].transpose()
	p[:,t] = Un[1,:].transpose()
	U = Un		





plt.plot(x1,p[:,int(T/dt)])
plt.show()


