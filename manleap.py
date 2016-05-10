import numpy as np
import sys
import time
start_time=time.time()
n=float(sys.argv[1])
m=float(sys.argv[2]) 
omega=0
K=1
Temp=30
eta=.1
epsilon=(2*eta*K*Temp/m)**(0.5)
dT=.001
g=9.8
N=7000000
a=1
def gamma():
    return 6*3.14*eta*a
def gamma_par(X):
    return gamma()*(1+(9.0/16)*(a/(X[1]+a))-1.0/8*(a/(X[1]+a))**3)
def gamma_per(X):
    return gamma()*(1+(9.0/8)*(a/(X[1]+a))-1.0/2*(a/(X[1]+a))**3)
#this is the true solution for contant gamma case, from balky's book .
def true (i):
        return (2*K*Temp/((gamma()**2)/m))*((gamma()/m)*i*dT-1+np.exp(-gamma()/m*i*dT))
#
def S2(X):
	return (gamma_per(X)/(m*eta))**(0.5)
def S1(X):
	return (gamma_par(X)/(m*eta))**(0.5)
def f1(X):
	return 0
def f2(X):
	return -g
X=np.zeros(shape=(N,2))
V=np.zeros(shape=(N,2))
#Initial conditions, first index is time and second (x,y)
X[0][:]=[0,100.0]
V[0][:]=[0,0]
for i in range (0,N-1):
	W=np.random.normal(0,1,2)*dT**(0.5)
	X_star=X[i][:]+0.5*V[i][:]*dT
	S1_star=S1(X_star)
	S2_star=S2(X_star)
	V[i+1][0]=((1+0.5*S1_star**2*dT)**(-1))*((1-0.5*S1_star**2*dT)*V[i][0]+f1(X_star)*dT+epsilon*S1_star*W[0])
	V[i+1][1]=((1+0.5*S2_star**2*dT)**(-1))*((1-0.5*S2_star**2*dT)*V[i][1]+f2(X_star)*dT+epsilon*S2_star*W[1])
	X[i+1][:]=X_star+0.5*V[i+1][:]*dT
	#If the particle hits the wall which is y=0
	if X[i+1][1]<0:
		X[i+1][1]=-1*X[i+1][1]
		V[i+1][1]=-1*V[i+1][1]
infile3=open('Mean_square_displacement_proper_%d_%d'%(n,m),'w') 
infile=open('varghese_%d'%n,'w')
infile2=open('Histogram_%d'%n,'w') 
x=np.zeros(shape=N)
for i in range (0,N):
       x[i]=X[i][1]
hist=np.histogram(x,1000)
mx=max(x)
mi=min(x)
for i in range (0,1000):
       infile2.write(str(hist[0][i]*1.0/(N*(mx-mi))*1000)+'\t'+str(hist[1][i])+'\n')
dz=.001
mean=0
z=0
for i in range (0,200000):
       mean=mean+(2*K*Temp)/(gamma_par([0,z]))*(m*g)/(K*Temp)*np.exp(-m*g/(K*Temp)*(z))*dz
       z=z+dz
print(mean)    
#exit()
A=10
for T in range (0,int(N/2),1000):
        msd=0
        if T/1000 > A:
                exit()
        if T%1000 == 0:
                print (T/1000)
        for i in range (0,N-(A+1)*1000):
                msd=msd+(X[i][0]-X[i+T][0])**2
        infile3.write(str(T*dT)+'\t'+str(msd/(N-(A+1)*1000))+'\t'+str(true(T))+'\n')#+str(msd_em/(N-(A+1)*1000))+'\n')
	
