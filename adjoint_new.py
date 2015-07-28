import numpy as np
import plot_tools as pt
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b


def F_newton(unp1,un,param):
	# function used for the newton iterations. It is basically the discretized burgers equations, with 0 as the rigth part.
	dx = param.dx
	dt = param.dt
	nu = param.nu

	F = np.zeros(un.size)
	F[1:-1] = (unp1[1:-1]-un[1:-1])/dt + unp1[1:-1]*(unp1[1:-1]-unp1[0:-2])/dx - nu * (unp1[0:-2] - 2*unp1[1:-1] + unp1[2:])/(dx*dx)
	F[0] = (unp1[0]-un[0])/dt + unp1[0]*(unp1[0]-unp1[-1])/dx - nu * (unp1[-1] - 2*unp1[0] + unp1[1])/(dx*dx)
	F[-1] = (unp1[-1]-un[-1])/dt + unp1[-1]*(unp1[-1]-unp1[-2])/dx - nu * (unp1[-2] - 2*unp1[-1] + unp1[0])/(dx*dx)
	return F

def Jacob_F_newton(unp1,param):
	# Jacobian of F_newton. Used for the newton iterations
	dx = param.dx
	dt = param.dt
	nu = param.nu

	ip1 = np.array(range(1,un.size))
	ip1 = np.hstack((ip1,0))
	im1 = np.array(range(0,un.size-1))
	im1 = np.hstack((un.size-1,im1))

	d = np.zeros(un.size)
	d += 1.0/dt - 2 * nu / (dx*dx)
	d[1:] += ( 2.0 * unp1[1:] - unp1[0:-1] ) / dx
	d[0] += ( 2.0 * unp1[0] - unp1[-1] ) / dx

	dp1 = np.zeros(un.size-1)
	dp1 += -nu / (dx * dx)

#	dm1 = np.zeros(un.size-1)
	dm1 = - unp1[1:] / dx
	dm1 += - nu / (dx * dx)

	J = np.diag(d) + np.diag(dp1,1) + np.diag(dm1,-1)
	J[0,-1] = - unp1[0]/dx - nu / (dx * dx)
	J[-1,0] = - nu / (dx * dx)
#	print(J[0:5,0:5])
#	print(J[-5:,-5:])
	return J

def Newton_step(un,param):
# time integration of burgers, with newton iterations
	unp1 = np.copy(un)
	for i in range(5):
		J = Jacob_F_newton(unp1,param)
		F = F_newton(unp1,un,param)
		unp1 = np.linalg.solve(J,np.dot(J,un)-F)
#		unp1[0] = unp1[-1]
	return unp1

def Plot(unp1,un):
	pt.multiplot1(False,(unp1,un),('-r','-k'))

def Observable(u):
	y = np.copy(u) + 0.01*np.random.normal(0,0.1,u.size)
	return y

class Parameters():
	def __init__(self,dx,dt,nx,nt,psize = 0, obs=False):
		self.dx = dx
		self.dt = dt
		self.nu = nu
		self.nx = nx
		self.nt = nt
		self.obs = obs
		self.psize = psize


def DGDP(param):
# g is the function which relates the parameters to the initial conditions. here is computed partial_p g 
	return np.eye(param.psize)

def DHDP(param):
# h is the function so that the cost function J = int h dt. 
#Here is computed partial_p h 
	return np.zeros(param.psize)

def DFDP(param):
# f is the function so that the cost function J = int h dt. 
#Here is computed partial_p h 
	return np.zeros((param.psize,param.psize))

def Compute_mu(lamda,param):
# here is computed the second lagrangian parameter
	return lamda[:,0]

def Compute_U(p,param):
	u0 = p
	U = np.zeros((nx,nt))
	U[:,0] = u0
	unp1 = np.copy(U[:,0])
	for i in range(1,nt):
		un = np.copy(unp1)
		unp1 = Newton_step(un,param)
		U[:,i] = unp1
	return U

#def DlDt(l,dfdup,dtdfdup,dfdu,dhdu):
#	dldt = np.dot(np.linalg.inv(dfdup.T),np.dot(dfdu.T - dtdfdup.T)*l +dhdu)
#	return dldt

def F_newton_l(lm1,l,dtdfdup,dfdup,dfdu,dhdu,param):
	dt = param.dt
	F = (l-lm1)/dt + np.dot(np.linalg.inv(dfdup.T),np.dot(dtdfdup.T,lm1) - np.dot(dfdu.T,lm1) - dhdu.T)
	return F

def Jacob_F_newton_l(dtdfdup,dfdu,param):
	dt = param.dt
	J = - np.diag(np.ones(param.nx)/dt) + dtdfdup.T - dfdu.T
	return J

def DHDU(u,y):
# here is the partial derivative of the functional with respect to u
	return 2*(u-y)

def DTDFDUP(param):
# here is the time derivation of the partial derivative of the main equations (burgers) with respect to dudt
	return np.zeros((param.nx,param.nx))

def DFDUP(param):
	return np.eye(param.nx)

def DFDU(u,param):
# here is the partial derivative of the main equations (burgers) with respect to u
	dx = param.dx
	dt = param.dt
	nu = param.nu

	dfdu = np.diag(u/dx) + 2 * np.diag(np.ones(u.size)/(dx*dx))
	dfdu += np.diag(np.ones(u.size-1)/dx,-1) - np.diag(nu * np.ones(u.size-1)/(dx*dx),-1) - np.diag(nu * np.ones(u.size-1)/(dx*dx),1)
	dfdu[0,-1] = 1/dx - nu / (dx*dx)
	dfdu[-1,0] = - nu / (dx*dx)
	return dfdu

def Compute_lambda(U,param):
# The adjoint equations are solved in this function.
	nt = param.nt
	Y = param.obs

	dfdup = DFDUP(param)
	dtdfdup = DTDFDUP(param)
	lamda = np.zeros((param.psize,nt))
	for j in range(nt-1,0,-1):
		dhdu = DHDU(U[:,j],Y[:,j])
		dfdu = DFDU(U[:,j],param)
		J = Jacob_F_newton_l(dtdfdup,dfdu,param)
		i = j-1
		lm1 = np.copy(lamda[:,j])
		for k in range(5):
			F = F_newton_l(lm1,lamda[:,j],dtdfdup,dfdup,dfdu,dhdu,param)
			
			try:
				lm1 = np.linalg.solve(J,np.dot(J,lamda[:,j])-F)
			except:
				print('l fail')
				
				jlmf = np.dot(J,lamda[:,j])-F
				np.save('J',J)
				np.save('jlmf',jlmf)
				print(J)
		lamda[:,i] = lm1
		s = 'lamda[:,' + str(i) + '] value: ' + str(lamda[:,i])
#		print(s)
	return lamda

def Grad_Functional(p,param):
# Here we compute finaly the gradient
	U = Compute_U(p,param)
	lamda = Compute_lambda(U,param)
	mu = Compute_mu(lamda,param)
	dJ = np.zeros(p.size)
	dgdp = DGDP(param)
	dJ += np.dot(mu.T,dgdp)
	for i in range(nt):
		dhdq = DHDP(param)
		dfdq = DFDP(param)
		dJ += dt*(dhdq + np.dot(lamda[:,i].T,dfdq))
		
		
	s = 'Gradiant value: ' + str(dJ)
#	print(s)
	return -10.0*dJ # WHY IT IS MUCH MORE EFFICIENT WITH A x10 ??

def Functionnal_p(p,param):
# here is the functionnal
	nt = param.nt
	dt = param.dt
	U = Compute_U(p,param)
	Y = param.obs
	J = 0.0
	for i in range(nt):
		J+= np.linalg.norm(U[:,i] - Y[:,i])*dt
	return J

if __name__ == '__main__':
	print "This only executes when %s is executed rather than imported" % __file__
	nx = 35
	nu = 0.0001
	xl = 0.0
	xr = 1.0
	x = np.linspace(xl,xr,nx)
	dx = (xr-xl)/nx
	dt = 0.05
	nt = 50
	T = dt*nt
	param = Parameters(dx,dt,nx,nt,nx)

	u0 = 0.2*( 1.05 + -np.sin( -np.pi + 2.0*np.pi*(x-xl)/(xr-xl) ) ) + 1.0*np.random.normal(0,0.1,x.size)
	U = np.zeros((nx,nt))
	U[:,0] = u0

	y0 = Observable(u0)
	Y = np.zeros((y0.size,nt))
	Y[:,0] = y0

	unp1 = np.copy(U[:,0])

	for i in range(1,nt):
		un = np.copy(unp1)
		unp1 = Newton_step(un,param)
	#	if np.mod(i,1) == 0:
	#		Plot(unp1,un)
		U[:,i] = unp1
		Y[:,i] = Observable(unp1)

	param.obs = Y
	p0 = 0.1*( 0.25 + -np.sin( -np.pi + 1.05*np.pi*(x-xl)/(xr-xl) ) ) + 0.3*np.random.normal(0,0.1,x.size)
#	p0 = 0.2*( 1.05 + -np.sin( -np.pi + 2.0*np.pi*(x-xl)/(xr-xl) ) ) + 0.2*np.random.normal(0,0.1,x.size)
	# brute force
#	ic_opt = fmin_bfgs(Functionnal_p,p0,None,(param,))

	pt.closeall()
	ic_opt = fmin_l_bfgs_b(Functionnal_p,p0,Grad_Functional,(param,),iprint=1)
	Plot(u0,p0)
	Plot(u0,p0)
	Plot(u0,ic_opt[0])
	#Plot(u0,p0)

