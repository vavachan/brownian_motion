import numpy as np
gamma=1
m=1
k=2
beta=1
delta=.01
F=np.zeros(shape=2)

X=np.zeros(shape=2)
def Func(X):
	W=np.random.normal(0,1)
	fun=(-gamma*X[1]-k/m*X[0])*delta+((2*gamma/beta*m)**2)*W*delta
	#print X
	F[1]=fun
	F[0]=X[1]*delta
	#print F
	return F
K1=np.zeros(shape=2)
K2=np.zeros(shape=2)
K3=np.zeros(shape=2)
K4=np.zeros(shape=2)
def RK(X):
	K1[:]=Func(X)
	K2[:]=Func(X+0.5*K1)
	K3[:]=Func(X+0.5*K2)
	K4[:]=Func(X+K3)
def update(X):
	RK(X)
	X[:]=X+1.0/6*(K1+2*K2+2*K3+K4)
N=200000
auto=np.zeros(shape=N)
hist=np.zeros(shape=1000)
infile=open("X_t.dat","w")
infile2=open("Corr.dat","w")
infile3=open("hist.dat","w")
MIN=0
for i in range (0,N):
	update(X)
	auto[i]=X[1]
	if X[1]<MIN:
		MIN=X[1]
	infile.write(str(i)+'\t'+str(X[0])+'\t'+str(X[1])+'\n')
	for j in range (0,1000):
		if X[0]<MIN+(j+1)*-2*MIN/1000 and X[0]>MIN+j*-2*MIN/1000:
			hist[j]=hist[j]+1
			break
for j in range (0,1000):
	infile3.write(str(MIN+(j+1)*-2*MIN/1000)+'\t'+str(hist[j]/N)+'\n')
sum1=0
for i in range (0,1000):
	sum1=sum1+hist[i]
print sum1/N	
exit()
for T in range (0,N/2,2):
	corr=0
	print T
	for i in range (0,N-T):
		corr=corr+(auto[i]*auto[i+T])**2
	corr=corr/(N-T)
	infile2.write(str(T)+'\t'+str(corr)+'\n')
	
