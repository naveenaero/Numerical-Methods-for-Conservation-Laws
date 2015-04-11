# Import necessary libraries
from numpy import *
from matplotlib import *
from pylab import *

#Take User input for Grid size, dt/dx and Time upto which the simulation will run
n = int(raw_input('Enter the number of points in the interval [-1,1]:'))
dtdx = float(raw_input('Enter the value of dt/dx:'))
index = int(raw_input('Enter the question number which you want to see? :'))
Time = float(raw_input('Enter the time at which you want to capture the plot(0 to 9):'))
filename = raw_input('Enter the filename with which you want to save the plot generated:');

#generate grid
x1 = linspace(-1,1,n)
dx = x1[1]-x1[0]			
dt = dx*dtdx
u = zeros((n,int(Time/dt)));	#matrix storing values of conserved variable
f = zeros((n,int(Time/dt)));	#matrix storing values of flux function 

#Initialisation for Case 1 & 2
for i in range(n):
	if x1[i] > -1.0/3.0 and x1[i] < 1.0/3.0:
		u[i,0] = 1

count  = 0;
position = 0;

#intialisation for case 3: Jump spread over 3 points
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
	for i in range(n):
		f[i,0] = u[i,0];
if index == 2:
	for i in range(n):
		f[i,0] = u[i,0]**2/2;
if index == 3:
	for i in range(n):
		f[i,0] = u[i,0]*(1-u[i,0]);


# Function for returning the sign of wave speed at each grid point
# Return 1 when +ve, -1 when -ve and 0 otherwise
def calc_sigma(i,flux,cons_variable):
	if cons_variable[i+1]-cons_variable[i]==0 or flux[i+1]-flux[i]==0:
		return 0;
	elif (flux[i+1]-flux[i])/(cons_variable[i+1]-cons_variable[i]) > 0:
		return 1;
	else:
		return -1;

# Function for applying the cyclic boundary condition
def apply_boundary_condn(u):
		a[0] = a[n-1];

# Function for updating the flux values after each time instant
def update_flux(fl,a,index):
	if index == 1:
		for i in range(n):
			fl[i] = a[i];
	if index == 2:
		for i in range(n):
			fl[i] = a[i]**2/2;
	if index == 3:
		for i in range(n):
			fl[i] = a[i]*(1-a[i]);

un = u[:,0];	#store intial condition in un - a temporary vector depicting the state at time t
fl = f[:,0]; 	#store initial flux in fl
a = un;			



for t in range(1,int(Time/dt)):			#Time Loop
	
	for x in range(n-1):				#Space Loop

		speed = calc_sigma(x,fl,un);
		flux_diff = fl[x+1] - fl[x];

		if speed > 0:
			a[x+1] = un[x+1] - dtdx*(flux_diff);
		if speed < 0:
			a[x] = un[x] - dtdx*(flux_diff);
		
	apply_boundary_condn(a);
	u[:,t] = a;
	un = a;
	update_flux(fl,a,index);

		
plot(x1,u[:,int(Time/dt)-1]);
xlabel('Space');
ylabel('Amplitude');
title('Q. ' + str(index) + ', grid points ' + str(n) + ', at time = ' +str(Time)+ ' seconds');
savefig(filename);
grid()
show()

