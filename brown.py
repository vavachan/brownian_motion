import numpy as np
#gamma=1
m=1
k=2
eta=10
a=.01
beta=1
delta=.01
D=3
l=.2
F=np.zeros(shape=2*D)
X=np.zeros(shape=2*D)
def gamma(X):
	return 6*3.14*eta*a 
	#return 6*3.14*eta*a*(1+9/8*(a/(X[0]+l)))
def Func(X):
	W=np.random.normal(0,1,D)
	for i in range (0,D):
		F[i+D]=(-gamma(X)*X[i+D]-k/m*X[i])*delta+((2*gamma(X)/beta*m)**2)*W[i]*delta
		F[i]=X[i+D]*delta
	return F
K1=np.zeros(shape=2*D)
K2=np.zeros(shape=2*D)
K3=np.zeros(shape=2*D)
K4=np.zeros(shape=2*D)
def RK(X):
	K1[:]=Func(X)
	K2[:]=Func(X+0.5*K1)
	K3[:]=Func(X+0.5*K2)
	K4[:]=Func(X+K3)
def update(X):
	RK(X)
	X[:]=X+1.0/6*(K1+2*K2+2*K3+K4)
N=5000
gap=100
auto=np.zeros(shape=(N,2*D))
hist=np.zeros(shape=(gap,2*D))
infile=open("X_t.dat","w")
infile2=open("Corr.dat","w")
infile3=open("hist.dat","w")
MIN=np.zeros(shape=2*D)
for i in range (0,N):
	update(X)
	auto[i][:]=X
	infile.write(str(i)+'\t')
	for p in range (0,2*D):
		infile.write(str(X[p])+'\t')
	infile.write('\n')
	continue
	for k in range (0,2*D):
		if X[k]<MIN[k]:
			MIN[k]=X[k]
		for j in range (0,gap):
			if X[k]<MIN[k]+(j+1)*-2*MIN[k]/gap and X[k]>MIN[k]+j*-2*MIN[k]/gap:
				hist[j][k]=hist[j][k]+1
				break
k=0
for j in range (0,gap):
	infile3.write(str(MIN[k]+(j+1)*-2*MIN[k]/gap)+'\t'+str(hist[j][k])+'\n')
sum1=0
#for i in range (0,gap):
#	sum1=sum1+hist[i]
#print sum1/N	
#exit()
for T in range (0,N/2,2):
	corr=0
	if T*1.0%1000 == 0.0:
		print T
	for i in range (0,N-T):
		corr=corr+(auto[i][k]*auto[i+T][k])
	corr=corr/(N-T)
	infile2.write(str(T)+'\t'+str(corr)+'\n')
	
