#	The equation I am giving here as an example is 
#		dv=A_v(x,v)dt+B_v(x,v)dt
#		dx=A_x(x,v)dt+B_x(x,v)dt
#	where A_v(x,v)=-gamma/m*v-omega^2*x
#	      B_v(x,v)=sqrt(2*gamma*KB*T)/m
#	      A_x(x,v)=v
#	      B_x(x,v)=0
# choose 1 for Runge-kutta and 2 for Euler-Maruyama 
m=0.5
omega=1
eta=1
a=1
KB=1.38*10**(-2)
T=300#
dt=.0005
M=3
D=2
N=2000
Xzero=0

import numpy as np
def true(i,w1,w2,w3):
	X=np.zeros(shape=D)
	X[0]=np.exp(-2*i*dt+w1-w2)*np.cos(w3)
	X[1]=np.exp(-2*i*dt+w1-w2)*np.sin(w3)
	return X 
def gamma():
	return 100
	#return 6*3.14*eta*a
def drift(X):
        temp=np.zeros(shape=D)
        temp[0]=-(3.0/2)*X[0]
	temp[1]=-(3.0/2)*X[1]
	return temp
def diffusion(X):
        temp=np.zeros(shape=(D,M))
        temp[0][:]=[X[0],-X[0],-X[1]]
	temp[1][:]=[X[1],-X[1],X[0]]
	#print X,temp
	return temp
def jac_diff(X):
        temp=np.zeros(shape=(D,M,D))
        temp[0][:][:]=[[1,0],[-1,0],[0,-1]]
	temp[1][:][:]=[[0,1],[0,-1],[1,0]]
	return temp
import sde


X_ini=[1.0,0.0]
solver = sde.ito(N,dt,drift,diffusion,jac_diff,D,M,1) 

solver.set_initial_condition(X_ini)
x=solver.solve()
X=np.zeros(shape=N)
V=np.zeros(shape=N)
infile=open("t_X'1_X'2_X1_X2.dat","w")
infile2=open("Histogram.dat","w")
infile3=open("Correlation.dat","w")
for i in range (0,N):
	X[i]=x[i][0]
	V[i]=x[i][1]
hist=np.histogram(X,1000)
#print hist	
for i in range (0,N):
        infile.write(str(i*dt)+'\t') # the following 3 lines writes all the information to a file.
        infile.write(str(solver.X_all[i][0])+'\t'+str(solver.X_all[i][1])+'\t'+str(true(i,solver.W[i][0],solver.W[i][1],solver.W[i][2])[0])+'\t'+str(true(i,solver.W[i][0],solver.W[i][1],solver.W[i][2])[1]))
        infile.write('\n')
exit()
for i in range (0,1000):
	infile2.write(str(hist[0][i]*1.0/N)+'\t'+str(hist[1][i])+'\n') 
for T in range (0,N/2,100):
        corr=0
#	if T/1000 > 20:
#		exit()
        if T%1000 == 0:
                print T/1000
        for i in range (0,N/2):
                corr=corr+(solver.X_all[i][0]-solver.X_all[i+T][0])**2
        infile3.write(str(T*dt)+'\t'+str(corr/(N/2))+'\n')
exit()

