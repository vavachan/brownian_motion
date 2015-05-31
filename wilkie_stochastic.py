import numpy as np
import const 
#import acor
#import matplotlib.pyplot as plt
from function import F
D=const.D # dimension
N=const.N
class wilkie_stoch():
	def __init__(self,iniX,N):
		self.iniX=iniX
		self.N=N
		self.X_all=np.zeros(shape=(self.N,2*D))
	def RK(self):
		K1=np.zeros(shape=2*D)
		K2=np.zeros(shape=2*D)
		K3=np.zeros(shape=2*D)
		K4=np.zeros(shape=2*D)
		X=self.iniX
		for i in range (0,self.N):
			K1[:]=F(X)
			K2[:]=F(X+0.5*K1)
			K3[:]=F(X+0.5*K2)
			K4[:]=F(X+K3)
			X[:]=X+1.0/6*(K1+2*K2+2*K3+K4)
			self.X_all[i][:]=X
X=np.zeros(shape=2*D)		
B=wilkie_stoch(X,N)
B.RK()
infile=open("all_data.dat","w")
for i in range (0,N):
	infile.write(str(i)+'\t') # the following 3 lines writes all the information to a file.
	for p in range (0,2*D):
		infile.write(str(B.X_all[i][p])+'\t')
	infile.write('\n')	
A=np.zeros(shape=N)
for i in range (0,N):
	A[i]=B.X_all[i][0]
#tau,mean,sigma=acor.acor(A)
#print tau
#plt.plot(A)
#plt.show()
