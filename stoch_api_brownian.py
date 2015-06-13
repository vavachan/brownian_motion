m=1
omega=0
omegax=0
omegay=15
eta=.1
a=1
KB=1
Temperature=1
dt=.0001
M=4
D=4
N=10000000
Xzero=0
import sys
Z_ini=float(sys.argv[1])
import numpy as np
def gamma():
	#return .5
	return 6*3.14*eta*a
def gamma_par(X):
        return gamma()/(1-(9.0/16)*(a/(X[1]+a))+1.0/8*(a/(X[1]+a))**3-(45.0/256)*(a/(X[1]+a))**4-(1.0/16)*(a/(X[1]+a))**5)
def gamma_per(X):
        return gamma()/(1-(9.0/8)*(a/(X[1]+a))+1.0/2*(a/(X[1]+a))**3-(57.0/100)*(a/(X[1]+a))**4-(1.0/5)*(a/(X[1]+a))**5+(7.0/100)*(a/(X[1]+a))**11-(1.0/25)*(a/(X[1]+a))**12)

def true (i):
	return 2*KB*Temperature/(m*gamma_par([0,Z_ini,0,0])**2)*(gamma_par([0,Z_ini,0,0])*i*dt-1+np.exp(-gamma_par([0,Z_ini,0,0])*i*dt))

omegas=(omegay**2-gamma_per([0,Z_ini,0,0])**2/4)**(.5)

def true_cvv(i):
        X=[0,Z_ini,0,0]
        return KB*Temperature/m*np.exp(-gamma_per(X)*(dt*i)/2)*(np.cos(omegas*dt*i)-gamma_per(X)/(2*omegas)*np.sin(omegas*dt*i))
def drift(X):
        temp=np.zeros(shape=D)
	visco=[[gamma_par(X),0],[0,gamma_per(X)]]
#        print gamma_par(X)/gamma(),' ',gamma_per(X)/gamma(),' ',X[1]/a,X[3]
        temp[2:4]=-1*np.dot(visco,X[2:4])-[(omegax**2)*X[0],omegay**2*(X[1]-Z_ini)]
   	temp[0:2]=X[2:4]
	return temp
def diffusion(X):
        temp=np.zeros(shape=(D,M))
	temp=[[0,0,0,0],[0,0,0,0],[0,0,(2*gamma_par(X)*KB*Temperature)**(0.5)/m,0],[0,0,0,(2*gamma_per(X)*KB*Temperature)**(0.5)/m]]
	return temp
def jac_diff(X):
        temp=np.zeros(shape=(D,M,D))
	return temp
import sde
print Z_ini
X_ini=[0.0,Z_ini,0.0,0.0]
solver = sde.ito(N,dt,drift,diffusion,jac_diff,D,M,2) 
solver.set_initial_condition(X_ini)
x=solver.solve()
X=np.zeros(shape=N)
V=np.zeros(shape=N)
#infile=open("t_X'1_X'2_X1_X2.dat","w")
l=float(sys.argv[1])
infile2=open("Histogram.dat","w")
exec "infile3=open('Mean_square_displacement_(%2f)','w')" %l 
exec "infile4=open('Velocity_autocorrelation_(%2f)','w')" %l
#infile5=open("Position_autocorrelation_underdamped_em.dat","w")
#nfile6=open("Position_velocity_autocorrelation_underdamped_em.dat","w")
for i in range (0,N):
	X[i]=x[i][1]
	V[i]=x[i][0]
hist=np.histogram(X,1000)
#print hist	
#for i in range (0,N):
#        infile.write(str(i*dt)+'\t') # the following 3 lines writes all the information to a file.
#        infile.write(str(solver.X_all[i][0])+'\t'+str(solver.X_all[i][1]))
#        infile.write('\n')
#exit()
for i in range (0,1000):
	infile2.write(str(hist[0][i]*1.0/N)+'\t'+str(hist[1][i])+'\n') 
#exit()
for T in range (0,N/2,100):
        corrv=0
	corrp=0
	corrpv=0
	msd=0
	if T/1000 > 10:
		exit()
        if T%1000 == 0:
                print T/1000
        for i in range (0,N-11000):
                corrv=corrv+(solver.X_all[i][3]*solver.X_all[i+T][3])
#		corrp=corrp+(solver.X_all[i][0]*solver.X_all[i+T][0])
#		corrpv=corrpv+(solver.X_all[i][0]*solver.X_all[i+T][1])
		msd=msd+(solver.X_all[i][0]-solver.X_all[i+T][0])**2
        infile4.write(str(T*dt)+'\t'+str(corrv/(N-11000))+'\t'+str(true_cvv(T))+'\n')
#	infile5.write(str(T*dt)+'\t'+str(corrp/(N-5000))+'\n')
	infile3.write(str(T*dt)+'\t'+str(msd/(N-11000))+'\t'+str(true(T))+'\n')
 #       infile6.write(str(T*dt)+'\t'+str(corrpv/(N-5000))+'\n')
exit()

