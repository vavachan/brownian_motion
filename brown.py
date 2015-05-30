import numpy as np
#gamma=1
m=1 #mass
k=2 #Strenght of harmonic potential
eta=10 #Viscosity coefficient
a=.01  #radius of brownian particle
beta=1 #!/K_B*Temp
delta=.01 #time step (delta T)
D=3 #Dimension
l=.2 # The Boundary is X=-l plane
F=np.zeros(shape=2*D) 
X=np.zeros(shape=2*D) #The first half of this array are position co-ordinates, Second half velcity co-ordinates in the same order.
def gamma(X): # this function gives the value of gamma, the friciton coefficient.  
	return 6*3.14*eta*a 
	#return 6*3.14*eta*a*(1+9/8*(a/(X[0]+l)))
def Func(X):# This function 
	W=np.random.normal(0,1,D)# the Gaussian-white noise, independent noise for each co-ordinate.
	for i in range (0,D): # The funtion in the runge-kutta scheme.
		F[i+D]=(-gamma(X)*X[i+D]-k/m*X[i])*delta+((2*gamma(X)/beta*m)**2)*W[i]*delta
		F[i]=X[i+D]*delta
	return F
K1=np.zeros(shape=2*D)
K2=np.zeros(shape=2*D)
K3=np.zeros(shape=2*D)
K4=np.zeros(shape=2*D)
def RK(X):#Runge-kutta
	K1[:]=Func(X)
	K2[:]=Func(X+0.5*K1)
	K3[:]=Func(X+0.5*K2)
	K4[:]=Func(X+K3)
def update(X):
	RK(X)
	X[:]=X+1.0/6*(K1+2*K2+2*K3+K4)
N=5000# the number of time steps 
gap=100 # the number of intervals while sampling for position or velocity distributions
auto=np.zeros(shape=(N,2*D))# This array will store all the informations , ie the vector X at each time step.
hist=np.zeros(shape=(gap,2*D))# For the histogram
infile=open("X_t.dat","w")
infile2=open("Corr.dat","w")
infile3=open("hist.dat","w")
MIN=np.zeros(shape=2*D)
for i in range (0,N):
	update(X)# Runge-Kutta 
	auto[i][:]=X
	infile.write(str(i)+'\t') # the following 3 lines writes all the information to a file.
	for p in range (0,2*D):
		infile.write(str(X[p])+'\t')
	infile.write('\n')
	continue # MAke this 'continue' a comment to run the code which will find histogram
	for k in range (0,2*D):
		if X[k]<MIN[k]:
			MIN[k]=X[k]
		for j in range (0,gap):
			if X[k]<MIN[k]+(j+1)*-2*MIN[k]/gap and X[k]>MIN[k]+j*-2*MIN[k]/gap:
				hist[j][k]=hist[j][k]+1
				break
k=0 # change this from 1 to 2D the write the corresponding distribution. k=0,1,2 for x,y,z. and k=3,4,5 for vx,vy,vz
for j in range (0,gap): # this would write the array hist to a file 
	infile3.write(str(MIN[k]+(j+1)*-2*MIN[k]/gap)+'\t'+str(hist[j][k])+'\n')
sum1=0
#for i in range (0,gap):
#	sum1=sum1+hist[i]
#print sum1/N	
#exit()

# to find correlation functions (same k holds here )
for T in range (0,N/2,2):
	corr=0
	if T*1.0%1000 == 0.0:
		print T
	for i in range (0,N-T):
		corr=corr+(auto[i][k]*auto[i+T][k])
	corr=corr/(N-T)
	infile2.write(str(T)+'\t'+str(corr)+'\n')
	
