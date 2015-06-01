#	The equation I am giving here as an example is 
#		dv=A_v(x,v)dt+B_v(x,v)dt
#		dx=A_x(x,v)dt+B_x(x,v)dt
#	where A_v(x,v)=-gamma/m*v-omega^2*x
#	      B_v(x,v)=sqrt(2*gamma*KB*T)/m
#	      A_x(x,v)=v
#	      B_x(x,v)=0
# choose 1 for Runge-kutta and 2 for Euler-Maruyama 
m=1
omega=0
eta=1
a=1
KB=1#.38*10**(-2)
T=1#
dt=.0001
D=2
N=20000
Xzero=0
def gamma():
	return 6*3.14*eta*a
def drift(X):
        temp=np.zeros(shape=D)
        for i in range (0,D/2):
                temp[i+D/2]=-1*gamma()*X[i+D/2]-omega**2*X[i]
                temp[i]=X[i+D/2]
        return temp
def diffusion(X):
        temp=np.zeros(shape=D)
        for i in range (0,D/2):
                temp[i+D/2]=((2*gamma()*KB*T)/m)**(.5)
        return temp
def jac_diff(X):
        temp=np.zeros(shape=(D,D))
        return temp
import sde
import numpy as np

X_ini=[0.0,0.0]
solver = sde.ito(N,dt,drift,diffusion,jac_diff,D,1) #last entry is the choice
solver.set_initial_condition(X_ini)
x=solver.solve()
infile=open("Time_position_velocity","w")
for i in range (0,N):
        infile.write(str(i*dt)+'\t') # the following 3 lines writes all the information to a file.
        infile.write(str(solver.X_all[i][0])+'\t'+str(solver.X_all[i][1])+'\t')
        infile.write('\n')
exit()
