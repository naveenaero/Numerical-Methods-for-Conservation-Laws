'''
AE617: Numerical Methods for Conservation Laws
ASSIGNMENT 1 : Code for numeically solving Hyperbolic Partial Differential Equations
Equations listed in Assignment 1 @ https://github.com/naveenaero/Numerical-Methods-for-Conservation-Laws
'''

#Import necessary files libraries for plotting and Animation
from numpy import *
from matplotlib import *
from pylab import *

#Take User input for Grid size, dt/dx and Time upto which the simulation will run
n = int(raw_input('Enter the number of points in the interval [-1,1]:'))
dtdx = float(raw_input('Enter the value of dt/dx:'))
index = int(raw_input('Enter the question number which you want to see? :'))
Time = float(raw_input('Enter the time at which you want to capture the plot(0 to 9):'))

#generate grid
x1 = linspace(-1,1,n)
dx = x1[1]-x1[0]			
dt = dx*dtdx
u = zeros((n,int(10/dt)));	#matrix storing values of conserved variable
f = zeros((n,int(10/dt)));	#matrix storing values of flux function 

print shape(x1),shape(u[:,int(4/dt)])

#Initialise Grid

for i in range(n):
	if x1[i] > -1.0/3.0 and x1[i] < 1.0/3.0:
		u[i,0] = 1

count  = 0;
position = 0;

if index == 3:
	for i in range(n):
		if u[i,0] == 1 and u[i-1,0] == 0:
			position = i;
		if u[i,0] == 1:
			count = count + 1;



	u[position,0] = 0.5;
	u[position+count,0] = 0.5;


# Initialise flux function
if index == 1:
	f = u;
if index == 2:
	f[:,0] = u[:,0]**2/2;
	
if index == 3:
	f[:,0] = u[:,0]*(1-u[:,0]);



# Function for returning the sign of wave speed at each grid point
# Return 1 when +ve, -1 when -ve and 0 otherwise
def calc_sigma(i,j,flux,cons_variable):
	if cons_variable[i+1,j]-cons_variable[i,j]==0 or flux[i+1,j]-flux[i,j]==0:
		return 0;
	elif (flux[i+1,j]-flux[i,j])/(cons_variable[i+1,j]-cons_variable[i,j]) > 0:
		return 1;
	else:
		return -1;

# Function for applying the cyclic boundary condition
def apply_boundary_condn(j,u):
	if j!=0:
		u[0,j] = u[n-1,j];

# Function for updating the flux values after each time instant
def update_flux(j,f,u,index):
	if index == 1:
		f[:,j+1] = u[:,j+1];
	if index == 2:
		f[:,j+1] = u[:,j+1]**2/2;
	if index == 3:
		f[:,j+1] = u[:,j+1]*(1-u[:,j+1]);




for t in range(int(9/dt)):				# Time Loop

	for x in range(0,n-1):				# Space Loop 

		speed = calc_sigma(x,t,f,u);

		if speed > 0:					# If speed is positive, then apply backward difference at the (i+1)th point
			u[x+1,t+1] = u[x+1,t] - dtdx*(f[x+1,t]-f[x,t]);
		if speed < 0:					# If speed is negative, the apply forward difference at the ith point
			u[x,t+1] = u[x,t] - dtdx*(f[x+1,t]-f[x,t])
		if speed == 0:					# If the speed is zero, then leave the (i+1)th point unchanges
			u[x+1,t+1] = u[x+1,t];
			# u[x,t+1] = u[x,t];

	apply_boundary_condn(t+1,u);		# Apply cyclic boundary condition
	update_flux(t,f,u,index);			# Update Flux at the end of each space loop

if Time <= 9:		
	plot(x1,u[:,int(Time/dt)])
	xlabel('Space');
	ylabel('Amplitude');
	title('Q. ' + str(index) + ', grid points ' + str(n) + ', at time = ' +str(Time)+ ' seconds');
	# savefig(filename);
	grid()
	show()
else:
	print "Time out of bounds"


		