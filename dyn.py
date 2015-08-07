import numpy as np
import plot_tools as pt
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b
import pickle as pk


def Noise_dyn(param):
	noise = param.dt * param.is_g * param.g * np.random.normal(0,1,param.nx)
	return noise

def Newton_step(un,param):
# time integration of burgers, with newton iterations
	unp1 = np.copy(un)
	for i in range(5):
		J = Jacob_F_newton(unp1,param)
		F = F_newton(unp1,un,param)
		noise = Noise_dyn(param)
		unp1 = np.linalg.solve(J,np.dot(J,un)-F) + noise
#		unp1[0] = unp1[-1]
	return unp1

def Compute_U(p,param):
	u0 = p
	nx = param.nx
	nt = param.nt
	
	U = np.zeros((nx,nt))
	U[:,0] = u0
	unp1 = np.copy(U[:,0])
	for i in range(1,nt):
		un = np.copy(unp1)
		unp1 = Newton_step(un,param)
		U[:,i] = unp1
	return U

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

	ip1 = np.array(range(1,unp1.size))
	ip1 = np.hstack((ip1,0))
	im1 = np.array(range(0,unp1.size-1))
	im1 = np.hstack((unp1.size-1,im1))

	d = np.zeros(unp1.size)
	d += 1.0/dt - 2 * nu / (dx*dx)
	d[1:] += ( 2.0 * unp1[1:] - unp1[0:-1] ) / dx
	d[0] += ( 2.0 * unp1[0] - unp1[-1] ) / dx

	dp1 = np.zeros(unp1.size-1)
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

def Noise_obs(param):
	noise = param.s * np.random.normal(0,1,param.nx)
	return noise

def Observable(u,param):
	noise = Noise_obs(param)
	y = np.copy(u) + noise
	return y


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

#def DlDt(l,dfdup,dtdfdup,dfdu,dhdu):
#	dldt = np.dot(np.linalg.inv(dfdup.T),np.dot(dfdu.T - dtdfdup.T)*l +dhdu)
#	return dldt

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

def Grad_dyn(X,param):
	
	pass


