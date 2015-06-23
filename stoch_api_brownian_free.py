m=1
omega=[0,0]
g=[0,.98]
eta=.1
a=1
KB=1
Temperature=25.0
dt=.001
M=4
D=4
N=40000000
Xzero=0
import sys
Z_ini=float(sys.argv[1])
import numpy as np
import math 
#import scipy
#from scipy import optimize 
#from optimize import curve_fit
print Z_ini
X_ini=[0.0,Z_ini,0.0,0.0]
def gamma():
	return 6*3.14*eta*a
def GA(X):
        GAM=np.zeros(shape=2)
        GAM[0]=gamma()*(1+(9.0/16)*(a/(X[1]+a))-1.0/8*(a/(X[1]+a))**3)#/(1-(9.0/16)*(a/(X[1]+a))+1.0/8*(a/(X[1]+a))**3-(45.0/256)*(a/(X[1]+a))**4-(1.0/16)*(a/(X[1]+a))**5)
        GAM[1]=gamma()*(1+(9.0/8)*(a/(X[1]+a))-1.0/2*(a/(X[1]+a))**3)#/(1-(9.0/8)*(a/(X[1]+a))+1.0/2*(a/(X[1]+a))**3-(57.0/100)*(a/(X[1]+a))**4-(1.0/5)*(a/(X[1]+a))**5+(7.0/100)*(a/(X[1]+a))**11-(1.0/25)*(a/(X[1]+a))**12)
        return GAM
def true (i):
	return (2*KB*Temperature/(m*(gamma())**2))*(gamma()*i*dt-1+np.exp(-gamma()*i*dt))

#gammas=(GA(X_ini)[1]**2-4*omega[1]**2)**(.5)
#def true_cvv(i):
#        X=[0,Z_ini,0,0]
#        return KB*Temperature/m*np.exp(-GA(X)[1]*(dt*i)/2)*(np.cosh((gammas/2.0)*dt*i)-GA(X)[1]/(gammas)*np.sinh((gammas/2.0)*dt*i))

#omegas=(omega[1]**2-GA([0,Z_ini,0,0])[1]**2/4)**(.5)
def true_cvv(i):
#        X=[0,Z_ini,0,0]
	return 0
#        return KB*Temperature/m*np.exp(-GA(X)[1]*(dt*i)/2)*(np.cos(omegas*dt*i)-GA(X)[1]/(2*omegas)*np.sin(omegas*dt*i))

def drift(X):
        temp=np.zeros(shape=D)
        for i in range (0,D/2):
                temp[i+D/2]=-1*GA(X)[i]/m*X[i+D/2]-(omega[i]**2)*(X[i]-X_ini[i])-g[i]
                temp[i]=X[i+D/2]
        return temp
def diffusion(X):
        temp=np.zeros(shape=(D,M))
        temp=[[0,0,0,0],[0,0,0,0],[0,0,(2*GA(X)[0]*KB*Temperature/m)**(0.5),0],[0,0,0,(2*GA(X)[1]*KB*Temperature/m)**(0.5)]]
        return temp

def jac_diff(X):
        temp=np.zeros(shape=(D,M,D))
	return temp
import sde
solver = sde.ito(N,dt,drift,diffusion,jac_diff,D,M,2) 
solver.set_initial_condition(X_ini)
x=solver.solve()
X=np.zeros(shape=N)
V=np.zeros(shape=N)
#infile=open("t_X'1_X'2_X1_X2.dat","w")
l=float(sys.argv[1])
o=float(sys.argv[2])
infile2=open("Histogram_1.dat","w")
exec "infile3=open('Mean_square_displacement_proper_(%2f)_%d','w')" %(l,o) 
exec "infile4=open('Velocity_autocorrelation_proper_(%2f)','w')" %l
infile5=open("Position_autocorrelation_underdamped_em.dat","w")
infile6=open("Position_velocity_autocorrelation_underdamped_em.dat","w")
for i in range (0,N):
	X[i]=x[i][1]
	V[i]=x[i][0]
hist=np.histogram(X,1000)
#print hist	
#or i in range (0,N):
#      infile.write(str(i*dt)+'\t') # the following 3 lines writes all the information to a file.
#      infile.write(str(solver.X_all[i][1])+'\t'+str(solver.X_all[i][3]))
#      infile.write('\n')
#exit()
#popt, pcov = curve_fit(f, hist[1], hist[0])
#print popt,pcov
for i in range (0,1000):
	infile2.write(str(hist[0][i]*1.0/N)+'\t'+str(hist[1][i])+'\n') 
dz=.001
mean=0
z=0
for i in range (0,200000):
	mean=mean+(2*KB*Temperature)/(GA([0,z])[0])*(m*g[1])/(KB*Temperature)*np.exp(-m*g[1]/(KB*Temperature)*(z))*dz
	z=z+dz
print mean	
#exit()
for T in range (0,N/2,500):
        corrv=0
	corrp=0
	corrpv=0
	msd=0
	if T/1000 > 4:
		exit()
        if T%1000 == 0:
                print T/1000
        for i in range (0,N-5000):
                corrv=corrv+(solver.X_all[i][3]*solver.X_all[i+T][3])
#		corrp=corrp+(solver.X_all[i][0]*solver.X_all[i+T][0])
#		corrpv=corrpv+(solver.X_all[i][2]*solver.X_all[i+T][3])
		msd=msd+(solver.X_all[i][0]-solver.X_all[i+T][0])**2
        infile4.write(str(T*dt)+'\t'+str(corrv/(N-5000))+'\t'+str(true_cvv(T))+'\n')
	#infile5.write(str(T*dt)+'\t'+str(corrp/(N-3000))+'\n')
	infile3.write(str(T*dt)+'\t'+str(msd/(N-5000))+'\t'+str(true(T))+'\n')
 #       infile6.write(str(T*dt)+'\t'+str(corrpv/(N-5000))+'\n')
exit()

