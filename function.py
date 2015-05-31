import const
m=const.m
omega=const.omega
eta=const.eta
a=const.a
KB=const.KB
T=const.T
dt=const.dt
D=const.D
import numpy as np
def gamma():
	return 6*3.14*eta*a 
def A(X):
	temp=np.zeros(shape=2*D)
	for i in range (0,D):
		temp[i+D]=-1*gamma()*X[i+D]-omega**2*X[i]
		temp[i]=X[i+D]
	return temp
def B(X):
	temp=np.zeros(shape=2*D)
	for i in range (0,D):
		temp[i+D]=((2*gamma()*KB*T)/m)**(.5)
	return temp
def B_der(X):
	temp=np.zeros(shape=(2*D,2*D))
	return temp
def F(X):
	B1=B(X)
	for i in range (0,2*D):
		w=np.random.normal(0,1,1)
		B1[i]=B1[i]*w[0]
	F1=(A(X)-np.dot(B_der(X),B(X)))*dt+B1*dt
	return F1
