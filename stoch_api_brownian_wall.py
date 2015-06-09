m=1
omegax=10
omegay=10
eta=1
a=.1
KB=1
Temperature=.5
dt=.0001
M=4
D=4
N=1000000
Xzero=0
shift=0
import numpy as np
def gamma():
        return 6*3.14*eta*a
omegas=(omegax**2-gamma()**2/4)**(.5)
def true(i):
        return KB*Temperature/m*np.exp(-gamma()*(dt*i)/2)*(np.cos(omegas*dt*i)-gamma()/(2*omegas)*np.sin(omegas*dt*i))

def gamma_par(X):
	#return .5
	return 6*3.14*eta*a*(1+(9.0/16)*(a*1.0/X[1])-(1.0/8)*(a/X[1])**3)
def gamma_per(X):
	#return .5
	return 6*3.14*eta*a*(1+(9.0/8)*(a*1.0/X[1])-(1.0/2)*(a/X[1])**3)
def drift(X):
        temp=np.zeros(shape=D)
	visco=[[gamma_par(X),0],[0,gamma_per(X)]]
	#print gamma_par(X),' ',gamma_per(X),' ',X[1],X[3]
	temp[2:4]=-1*np.dot(visco,X[2:4])-[(omegax**2)*X[0],omegay**2*(X[1]-Z_ini)]#-[0,m*g]#[0.0,force(X[1])]
	temp[0:2]=X[2:4]
	return temp
def diffusion(X):
        temp=np.zeros(shape=(D,M))
	temp=[[0,0,0,0],[0,0,0,0],[0,0,(2*gamma_par(X)*KB*Temperature)**(0.5)/m,0],[0,0,0,(2*gamma_per(X)*KB*Temperature)**(0.5)/m]]
	return temp
def jac_diff(X):
        temp=np.zeros(shape=(D,M,D))
	return temp
import sde_wall
Z_ini=1000
X_ini=[0.0,Z_ini,0.0,0.0]
solver = sde_wall.ito(N,dt,drift,diffusion,jac_diff,D,M,2) 
solver.set_initial_condition(X_ini)
x=solver.solve()
X=np.zeros(shape=N)
Z=np.zeros(shape=N)
l=0
exec "infile=open('Trajectory_%2d.dat','w')"%l 
exec "infile2=open('Histogram_%2d.dat','w')"%l
exec "infile3=open('msd_z_%2d','w')" %l
exec "infile4=open('msd_x_%2d','w')"%l
exec "infile5=open('czz_%2d','w')"%l
exec "infile6=open('cxx_%2d','w')"%l
exec "infile7=open('cvxvx_%2d','w')"%l
exec "infile8=open('cvzvz_%2d','w')"%l
for i in range (0,N):
	X[i]=x[i][0]
	Z[i]=x[i][1]
hist_z=np.histogram(Z,1000)
hist_x=np.histogram(X,1000)
#print hist
#exit()	
for i in range (0,N):
        infile.write(str(i*dt)+'\t') # the following 3 lines writes all the information to a file.
        infile.write(str(solver.X_all[i][0])+'\t'+str(solver.X_all[i][1]))
        infile.write('\n')
#exit()
for i in range (0,1000):
	infile2.write(str(hist_x[0][i]*1.0/N)+'\t'+str(hist_z[0][i]*1.0/N)+'\t'+str(hist_x[1][i])+'\t'+str(hist_z[1][i])+'\n') 
#exit()
A=10
for T in range (0,N,100):
        cvzvz=0
	cvxvx=0
	cxx=0
	czz=0
	msdx=0
	msdz=0
	if T/1000 > A:
		exit()
        if T%1000 == 0:
                print T/1000
        for i in range (0,N-(A+1)*1000):
		cvxvx=cvxvx+(solver.X_all[i][2]*solver.X_all[i+T][2])
                cvzvz=cvzvz+(solver.X_all[i][3]*solver.X_all[i+T][3])
		cxx=cxx+(solver.X_all[i][0]*solver.X_all[i+T][0])
		czz=czz+((solver.X_all[i][1]-solver.X_all[0][1])*((solver.X_all[i+T][1])-solver.X_all[0][1]))
		msdz=msdz+(solver.X_all[i][1]-solver.X_all[i+T][1])**2
		msdx=msdx+(solver.X_all[i][0]-solver.X_all[i+T][0])**2
        infile3.write(str(T*dt)+'\t'+str(msdz/(N-(A+1)*1000))+'\n')
	infile4.write(str(T*dt)+'\t'+str(msdx/(N-(A+1)*1000))+'\n')
	infile5.write(str(T*dt)+'\t'+str(czz/(N-(A+1)*1000))+'\n')
        infile6.write(str(T*dt)+'\t'+str(cxx/(N-(A+1)*1000))+'\n')
	infile7.write(str(T*dt)+'\t'+str(cvxvx/(N-(A+1)*1000))+'\n')
	infile8.write(str(T*dt)+'\t'+str(cvzvz/(N-(A+1)*1000))+'\t'+str(true(T))+'\n')
exit()

