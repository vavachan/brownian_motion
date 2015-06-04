m=1
omega=1
eta=1
a=1
KB=1
T=1
dt=.005
M=2
D=2
N=1000000
Xzero=0

import numpy as np
def gamma():
	return .5
	#return 6*3.14*eta*a
def drift(X):
        temp=np.zeros(shape=D)
	for i in range (0,D/2):
                temp[i+D/2]=-1*gamma()*X[i+D/2]-omega**2*X[i]
                temp[i]=X[i+D/2]

	return temp
def diffusion(X):
        temp=np.zeros(shape=(D,M))
	temp=[[0,0],[0,(2*gamma()*KB*T)**(0.5)/m]]
	return temp
def jac_diff(X):
        temp=np.zeros(shape=(D,M,D))
	return temp
import sde

X_ini=[0.0,0.0]
solver = sde.ito(N,dt,drift,diffusion,jac_diff,D,M,2) 
solver.set_initial_condition(X_ini)
x=solver.solve()
X=np.zeros(shape=N)
V=np.zeros(shape=N)
infile=open("t_X'1_X'2_X1_X2.dat","w")
infile2=open("Histogram.dat","w")
infile3=open("Mean_square_displacement_underdamped_em.dat","w")
infile4=open("Velocity_autocorrelation_underamped_em.dat","w")
infile5=open("Position_autocorrelation_underdamped_em.dat","w")
infile6=open("Position_velocity_autocorrelation_underdamped_em.dat","w")
for i in range (0,N):
	X[i]=x[i][0]
	V[i]=x[i][1]
hist=np.histogram(X,1000)
#print hist	
for i in range (0,N):
        infile.write(str(i*dt)+'\t') # the following 3 lines writes all the information to a file.
        infile.write(str(solver.X_all[i][0])+'\t'+str(solver.X_all[i][1]))
        infile.write('\n')
#exit()
for i in range (0,1000):
	infile2.write(str(hist[0][i]*1.0/N)+'\t'+str(hist[1][i])+'\n') 
for T in range (0,N/2,100):
        corrv=0
	corrp=0
	corrpv=0
	msd=0
	if T/1000 > 3:
		exit()
        if T%1000 == 0:
                print T/1000
        for i in range (0,N-5000):
                corrv=corrv+(solver.X_all[i][1]*solver.X_all[i+T][1])
		corrp=corrp+(solver.X_all[i][0]*solver.X_all[i+T][0])
		corrpv=corrpv+(solver.X_all[i][0]*solver.X_all[i+T][1])
		msd=msd+(solver.X_all[i][0]-solver.X_all[i+T][0])**2
        infile4.write(str(T*dt)+'\t'+str(corrv/(N-5000))+'\n')
	infile5.write(str(T*dt)+'\t'+str(corrp/(N-5000))+'\n')
	infile3.write(str(T*dt)+'\t'+str(msd/(N-5000))+'\n')
        infile6.write(str(T*dt)+'\t'+str(corrpv/(N-5000))+'\n')
exit()

