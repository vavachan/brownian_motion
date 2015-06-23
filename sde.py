import numpy as np

class ito ():
	def __init__(self,N,dt,A,B,B_der,D,M,choice=1):
		self.choice=choice
		self.M=M
		self.D=D
                self.A=A
		self.B=B
		self.B_der=B_der
                self.N=N
		self.dt=dt
                self.X_all=np.zeros(shape=(self.N,self.D))
                self.dW=np.zeros(shape=(self.N,self.M))
                self.W=np.zeros(shape=(self.N+1,M))
	def multiply(self,A,B):
		temp=np.zeros(shape=self.D)
		for j in range (0,self.D):
			ten=(B[j][:][:])
			ten1=np.transpose(ten)
			sum2=0
			for k in range (0,self.M):
				sum1=0
				for i in range (0,self.D):
					sum1=sum1+A[k][i]*ten1[i][k]
				sum2=sum2+sum1
			temp[j]=sum2
		return temp	
			
	def F(self,X,W):
		D=self.D
		M=self.M
	        B1=np.zeros(shape=(D,M))
		K=np.zeros(shape=(D,M))
	        K=self.B(X)
	        B1=np.dot(K,W)
		F1=(self.A(X)-0.5*self.multiply((np.transpose(self.B(X))),self.B_der(X)))*self.dt+B1
		return F1
	def weiner(self):
		M=self.M
                for i in range (0,self.N):
                        self.dW[i][:]=self.dt**(.5)*np.random.normal(0,1,M)
                for j in range (0,M):
                        sum1=0
                        for i in range (1,self.N+1):
                                sum1=sum1+self.dW[i-1][j]
                                self.W[i][j]=sum1
        def RK(self):
		self.weiner()
		D=self.D
                K1=np.zeros(shape=D)
                K2=np.zeros(shape=D)
                K3=np.zeros(shape=D)
                K4=np.zeros(shape=D)
                X=np.zeros(shape=D)
                X[:]=self.X_all[0][:]
                for i in range (1,self.N):
                        K1[:]=self.F(X,self.dW[i-1][:])
                        K2[:]=self.F(X+0.5*K1,self.dW[i-1][:])
                        K3[:]=self.F(X+0.5*K2,self.dW[i-1][:])
                        K4[:]=self.F(X+K3,self.dW[i-1][:])
                        X[:]=X+1.0/6*(K1+2*K2+2*K3+K4)
                        #print self.dW[i][:]
                        self.X_all[i][:]=X
       
	def em(self):
		self.weiner()
		B1=np.zeros(shape=self.D)
		winc=np.zeros(shape=self.M)
		K=np.zeros(shape=(self.D,self.M))
		for i in range (1,self.N):
			K=self.B(self.X_all[i-1][:])
		        winc[:]=self.W[i][:]-self.W[i-1][:]
		        B1=np.dot(K,winc)
			self.X_all[i][:] = self.X_all[i-1][:] + self.dt*self.A(self.X_all[i-1][:]) + B1
			if self.X_all[i][1]<.2:
				self.X_all[i][1]=.4-self.X_all[i][1]
				self.X_all[i][3]=-1*self.X_all[i][3]
	def solve(self):
		if self.choice == 1:
			print 'RK'
			self.RK()
		if self.choice == 2:
			print 'Euler-Maruyama'
			self.em()
		if self.choice != 1 and self.choice != 2: 
			print 'wrong choice',self.choice
			exit()
		return self.X_all

	def set_initial_condition(self,X_ini):
		self.X_all[0][:]=X_ini
