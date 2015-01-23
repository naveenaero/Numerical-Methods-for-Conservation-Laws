from numpy import *
from matplotlib import *
from pylab import *

n1 = 40
n2 = 80
x1 = linspace(-1,1,40);
x2 = linspace(-1,1,80);
dx1 = x1[1]-x1[0];			
dx2 = x2[1]-x2[0];

dtdx = 0.8

dt1 = dx1*0.8
dt2 = dx2*0.8

u = zeros((40,int(100/dt1)));

u2 = zeros((80,int(4/dt2)));
for i in range(80):
	if i < 40:
		if x1[i] > -1.0/3.0 and x1[i] < 1.0/3.0:
			u[i,0] = 1
	if x2[i] > -1.0/3.0 and x2[i] < 1.0/3.0:
		u2[i,0] = 1

f = u

def sigma(i,j,flux,cons_variable):
	if cons_variable[i+1,j]-cons_variable[i,j]==0 or flux[i+1,j]-flux[i,j]==0:
		return 0;
	elif (flux[i+1,j]-flux[i,j])/(cons_variable[i+1,j]-cons_variable[i,j]) > 0:
		return 1;
	else:
		return -1;

def apply_boundary_condn(i,j,u):
	if i==0 and j!=0:
		u[i,j] = u[i+39,j-1];


# def propagate(x,t,speed,u):
# 	if speed >


for t in range(int(80/dt1)):
	for x in range(0,n1-1):
		apply_boundary_condn(x,t,u);
		speed = sigma(x,t,f,u);
		if speed >= 0:
			u[x+1,t+1] = u[x+1,t] - dtdx*(u[x+1,t]-u[x,t]);
		if speed < 0:
			u[x,t+1] = u[x,t] - dtdx*(u[x+1,t]-u[x,t])
		

#print u[:,4]
plot(x1,u[:,int(4/dt1)])
# plot(x2,u2)
grid()
show()


		