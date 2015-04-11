from numpy import *
from matplotlib import *
from pylab import *
from numpy.linalg import inv


gam = 1.4
gp1 = gam+1
gm1 = gam-1


n1 = input("Enter the number of points: ");

x = linspace(0,1,n1)

dx = x[1] - x[0]
u = zeros(n1)
rho = zeros(n1)
p = zeros(n1)
U = zeros((3,n1))
Un = zeros((3,n1))
f = zeros((3,n1))
vect = zeros((3,n1))
a = 0
T = 0
index = 0

# initialisation of the conserved variables and flux vector
for i in range(n1):
	if x[i]<=0.3:
		u[i] = 0.75
		p[i] = 1
		rho[i] = 1
		index = i
	else:

		u[i] = 0.125
		p[i] = 0.1
		rho[i] = 0.125
	U[0,i] = rho[i]
	U[1,i] = rho[i]*u[i]
	U[2,i] = p[i]/gm1 + 0.5*rho[i]*u[i]**2	# energy (e)
	f[0,i] = rho[i]*u[i]
	f[1,i] = rho[i]*u[i]**2 + p[i]
	f[2,i] = u[i]*(p[i]+U[2,i])				# use U[2,i] instead of energy
	vect[0,i] = u[i];
	vect[1,i] = p[i];
	vect[2,i] = rho[i];

#definfing the temporary variables
rhat = zeros((3,n1))			
Rp1 = zeros((3,3))				
Rm1 = zeros((3,3))
R = zeros((3,3))
lambdahat = zeros((3,n1))	
eigen_valuep1 = zeros((3,1))	
eigen_valuem1 = zeros((3,1))	
eigen_value = zeros((3,1))	
alphap1 = zeros((3,1))
alpham1	= zeros((3,1))
fp1 = zeros((3,1))
fm1 = zeros((3,1))
gamma = zeros(3)
deltaUp = zeros(3)
deltaUm = zeros(3)
a = 0


Hhatp1 = 0	# Hhat[j+1/2]
Hhatm1 = 0	# Hhat[j-1/2]
uhatp1 = 0	# uhat[j+1/2]
uhatm1 = 0	# uhat[j-1/2]
Hp1 = 0		# H[j+1]
Hm1 = 0		# H[j-1]
H = 0		# H[j]


##########################
fhatp1 = zeros((3,1))
fhatm1 = zeros((3,1))

##########################



'''
Returns the least value of time step
by fixing the CFL value at 0.8
'''
def time_step(vect,dx):

	t_step = zeros((3,1))
	for i in range(n1):
		
		a1 = (gam*vect[1,i]/vect[2,i])**0.5
		t_step[0] = 0.8*dx/(vect[1,i]-a1)
		t_step[1] = 0.8*dx/(vect[1,i])
		t_step[2] = 0.8*dx/(vect[1,i]+a1)
		if i==0:
			temp = amin(abs(t_step))
		else:
			if amin(abs(t_step))<temp:
				temp = amin(abs(t_step))
		
	return temp

'''
Returns the eigen-values and eigen-vector given the vector 'vect' 
and index 'l'. vect[:,l] contains u[l], p[l] and rho[l].
'''
def calc_eigen(vect):					
	R1 = zeros((3,3))
	eigen_value1 = zeros(3)
	a1 = (gam*vect[1]/vect[2])**0.5
	H1 = gam*vect[1]/((gam-1)*vect[2]) + vect[0]**2/2
	R1[:,0] = [1, vect[0]-a1, H1-vect[0]*a1]
	R1[:,1] = [1, vect[0], vect[0]**2/2]
	R1[:,2] = [1, vect[0]+a1, H1+vect[0]*a1]
	eigen_value1 = [vect[0]-a1, vect[0], vect[0]+a1]

	return R1,eigen_value1

Un = U
T=0
n=0
'''
Space and Time loop for computing
the Local Lax friedrich scheme
'''
while T<0.25:
	dt = time_step(vect,dx)
	print "time step:",dt
	for j in range(n1):
		a = (gam*vect[1,j]/vect[2,j])**0.5
		
		if j > 0 and j<n1-1:
			(R,eigen_value) = calc_eigen(vect[:,j])
			(Rp1,eigen_valuep1) = calc_eigen(vect[:,j+1])
			(Rm1,eigen_valuem1) = calc_eigen(vect[:,j-1])
			alphap1 = [ max(abs(eigen_value[0]),abs(eigen_valuep1[0])),
						max(abs(eigen_value[1]),abs(eigen_valuep1[1])),
						max(abs(eigen_value[2]),abs(eigen_valuep1[2])) ]

			alpham1 = [ max(abs(eigen_value[0]),abs(eigen_valuem1[0])),
						max(abs(eigen_value[1]),abs(eigen_valuem1[1])),
						max(abs(eigen_value[2]),abs(eigen_valuem1[2])) ]

			deltaUp = U[:,j+1] - U[:,j]
			deltaUm = U[:,j] - U[:,j-1]
			gammap = dot(inv(R),deltaUp)
			gammam = dot(inv(R),deltaUm)
			# fp1 = (f[:,j]+f[:,j+1])/2 - 0.5*(alphap1[0]*gammap[0]*R[:,0] + alphap1[1]*gammap[1]*R[:,1] + alphap1[2]*gammap[2]*R[:,2])
			# fm1 = (f[:,j]+f[:,j-1])/2 - 0.5*(alpham1[0]*gammam[0]*R[:,0] + alpham1[1]*gammam[1]*R[:,1] + alpham1[2]*gammam[2]*R[:,2])

			fp1 = (f[:,j]+f[:,j+1])/2 - 0.5*(amax(alphap1)*deltaUp)
			fm1 = (f[:,j]+f[:,j-1])/2 - 0.5*(amax(alpham1)*deltaUm)


			Un[:,j] = (U[:,j-1]+U[:,j+1])/2 - (0.5*dt/dx)*(fp1-fm1)

			# if j==index:
			# 	H = gam*vect[1,j]/((gam-1)*vect[2,j]) + vect[0,j]**2/2
			# 	H2 = gam*vect[1,j+1]/((gam-1)*vect[2,j+1]) + vect[0,j+1]**2/2
			# 	print "a:", a
			# 	print "H:", H
			# 	print "H1:", H2
			# 	print "R:",R
			# 	print "R1:",Rp1
			# 	print "eigen_value:", eigen_value
			# 	print "eigen_value1:", eigen_valuep1
			# 	print "Un:", Un[:,j]
			# 	print "U:", U[:,j]
			# 	print "deltaUp:",deltaUp
			# 	print "deltaUm", deltaUm
			# 	print "gammap", gammap
			# 	print "gammam", gammam
			# 	print "fp1:", fp1
			# 	print "fm1:", fm1
			# 	print "vect:", vect[:,j] 
			# 	print "u^2/2:", vect[0,j]**2/2
			# 	print "u:", vect[0,j]
			# 	print "R_inv:", inv(R)
			# 	print "alphap1:", alphap1
			# 	print "alpham1:", alpham1


	
		# elif j==0:
		# 	(R,eigen_value) = calc_eigen(vect[:,j])
		# 	(Rp1,eigen_valuep1) = calc_eigen(vect[:,j+1])
		# 	alphap1 = [ max(abs(eigen_value[0]),abs(eigen_valuep1[0])),
		# 				max(abs(eigen_value[1]),abs(eigen_valuep1[1])),
		# 				max(abs(eigen_value[2]),abs(eigen_valuep1[2])) ]

		# 	deltaUp = U[:,j+1] - U[:,j]
		# 	gammap = dot(inv(R),deltaUp)
		# 	# fp1 = (f[:,j]+f[:,j+1])/2 - 0.5*(alphap1[0]*gammap[0]*R[:,0] + alphap1[1]*gammap[1]*R[:,1] + alphap1[2]*gammap[2]*R[:,2])
		# 	fp1 = (f[:,j]+f[:,j+1])/2 - 0.5*(amax(alphap1)*(U[:,j+1]-U[:,j]))

		# 	Un[:,j] = (0+U[:,j+1])/2 - (a*dt/dx)*(fp1-0)
	
		# else:
		# 	(R,eigen_value) = calc_eigen(vect[:,j])
		# 	(Rm1,eigen_valuem1) = calc_eigen(vect[:,j-1])
		# 	alpham1 = [ max(abs(eigen_value[0]),abs(eigen_valuem1[0])),
		# 				max(abs(eigen_value[1]),abs(eigen_valuem1[1])),
		# 				max(abs(eigen_value[2]),abs(eigen_valuem1[2])) ]

		# 	deltaUm = U[:,j] - U[:,j-1]
		# 	gammam = dot(inv(R),deltaUm)
		# 	fm1 = (f[:,j]+f[:,j-1])/2 - 0.5*(alpham1[0]*gammam[0]*R[:,0] + alpham1[1]*gammam[1]*R[:,1] + alpham1[2]*gammam[2]*R[:,2])

		# 	Un[:,j] = (0+U[:,j-1])/2 - (a*dt/dx)*(0-fm1)

	T = T+dt
	print T
	U = Un
	vect[2,:] = Un[0,:]								# updating rho in vect
	vect[0,:] = Un[1,:]/Un[0,:]						# updating u in vect (divide the first 2nd row of Un by 1st row)
	vect[1,:] = gm1*(Un[2,:] - 0.5*Un[1,:]*vect[0,:]) # updating p in vect
	f[0,:] = vect[0,:]*vect[2,:]					# rho*u
	f[1,:] = vect[1,:] + vect[2,:]*vect[0,:]**2		# p + rho*u^2
	f[2,:] = vect[0,:]*(vect[1,:]+Un[2,:])			# use Un[2,:] instead of energy
	n=n+1



			# Hp1 = gam*p[j+1]/((gam-1)*rho[j+1]) + u[j+1]**2/2
			
# 			uhatp1 = (sqrt(rho[j])*u[j] + sqrt(rho[j+1])*u[j+1])/(sqrt(rho[j]) + sqrt(rho[j+1]))
# 			Hhatp1 = (sqrt(rho[j])*H + sqrt(rho[j+1])*Hp1)/(sqrt(rho[j])+sqrt(rho[j+1]))
			
			
# 			rhat[:,0] = [1, uhatp1-a, Hhatp1 - uhatp1*a]
# 			lambdahat[0] = uhatp1-a

# 			rhat[:,1] = [1, uhatp1, uhatp1**2/2]
# 			lambdahat[1] = uhatp1

# 			rhat[:,2] = [1, uhatp1+a, Hhatp1 - uhatp1*a]
# 			lambdahat[2] = uhatp1+a

# 		maxlambda = amax(abs(lambdahat))
		
# 		Un[:,j] = (U[:,j]+U[:,j+1])/2 - 



plot(x,vect[0,:])
show()
