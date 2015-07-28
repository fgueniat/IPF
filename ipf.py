import numpy as np
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b
from dyn import *


def P_int(f,X,Xp1,param):
# proba that (Xp05,Xp1) corresponds to transport of X by f, after a time dt.
	g=param.g
	dt = param.dt
	if np.abs(g) < 1.0e-12:
		g=1.0

	P = np.linalg.norm( Xp1 - X - dt * f(X,param) )**2.0 / (2.0 * dt * g * g)
#	print(P)
	return P
#	return 0.0

def P_obs(h,X,y,param):
# proba that the measure of X corresponds to observation y
	s=param.s
	if np.abs(s) < 1.0e-12:
		s=1.0
	P = np.linalg.norm(h(X,param) - y)**2.0 / (2.0  * s * s)
	return P

def F_intermediary(u,param):
# cancel noise for that part
	g = param.g
	param.g = 0.0
	unp1 = Newton_step(u,param)
	param.g = g
	return (unp1-u)/param.dt

def Observable_intermediary(u,param):
# cancel noise for that part
	s = param.s
	param.s = 0.0
	y = Observable(u,param)
	param.s = s
	return y

def Functionnal_IPF(p,param):
# Functional to minimize
	F = 0.0
	h = Observable_intermediary
	f = F_intermediary
	nt = param.nt
	X = p.reshape(param.nx,param.nt)
	for i in range(nt-1): #DA: t=0 included
		F += param.imp_int[i]*P_int(F_intermediary,X[:,i],X[:,i+1],param)
		F += param.imp_obs[i]*P_obs(h,X[:,i],param.obs[:,i],param)
	return F

def Ersatz(l,param):
# Ersatz for solving the algebraic equation
	X = param.xmin + l*param.Lnu
	X = X.reshape(param.nx,param.nt)
	return X

def Functionnal_IPF_lin(l,param):
# 1D functional, when X is replaced with the ersatz
	F = 0.0
	h = Observable
	f = F_intermediary
	nt = param.nt
	X = Ersatz(l,param)
	for i in range(nt-1): #DA: t=0 included
		F += param.imp_int[i]*P_int(F_intermediary,X[:,i],X[:,i+1],param)
		F += param.imp_obs[i]*P_obs(h,X[:,i],param.obs[:,i],param)
	return np.abs(F - param.fmin - param.epseps2)

def Id_Weight(lamda,param):
# log (to minimize risks of overflow) of weights associated with particles are identified here
		rho = np.dot(param.eps.T,param.eps)
		#numeric estimation of d lamda / d rho
		dl = 0.01
		Fp = Functionnal_IPF(Ersatz(lamda+dl,param),param)
		Fm = Functionnal_IPF(Ersatz(lamda,param),param)
		dldp = dl/(2.0*(Fp-Fm))
		# jacobian of the map
		(sm,logdet) = np.linalg.slogdet(param.L)
		logJ = np.log(2.0) + logdet + np.log(rho) * (1 - param.xmin.size) + np.log(  np.abs( ( lamda**(param.nt-1.0) ) * dldp)  )
#			s = 'J: ' + str(np.abs(np.linalg.det(self.L))) + '  ' + str(rho) + '    ' + str(np.abs( (lamda**(2.0*self.ndim-1.0)))) + '  ' + str(dldp)
#			print(s)
		#weight
		s = "log det L :" + str(logdet) + " rho: "+ str(rho) + " dldp: " + str(dldp)  
		weight = -param.fmin +  logJ
		print(s)
		return weight

def Weight_rescale(w):
# log of weights are rescaled, so a log representation is no necessary anymore
	w2 = np.zeros(w.size)
	w2 = np.exp( w-np.min(w) )/np.sum(  np.exp( w-np.min(w) )  )
	return w2


def Hessianm1(X0,param,dx=0.01):
# Hessian approximation. This is really dirty...
	f = Functionnal_IPF
	X0 = X0.flatten().reshape(param.nx,param.nt)
	ndim = X0.size
	H = np.zeros((ndim,ndim))
	for i1 in range(0,ndim):
		s = "Hessian computations ......... " + str(np.round(100.0*i1**2/ndim**2)) + "%"
#		print(s)
		dX1=np.zeros(ndim)
		dX1[i1] = dx
		dX1 = dX1.reshape(param.nx,param.nt)
		for i2 in range(0,i1+1):
			dX2=np.zeros(ndim)
			dX2[i2] = dx
			dX2 = dX2.reshape(param.nx,param.nt)
			
			fpp = f(X0+(dX1+dX2)/2.0,param)
			fpm = f(X0+(dX1-dX2)/2.0,param)
			fmp = f(X0+(-dX1+dX2)/2.0,param)
			fmm = f(X0-(dX1+dX2)/2.0,param)
			d2fdx1dx2 = (fpp-fmp -(fpm-fmm))/(dx*dx)
			H[i1,i2]=d2fdx1dx2
			H[i2,i1]=d2fdx1dx2

	return np.linalg.pinv(H)

def Grad_Functional_IPF(X,param):
# TO DO: derive directly the hessian and the gradient of the functional ?
	nu = param.nu
	nx = param.nx
	nt = param.nt
	DJ = np.zeros(X.size)
	for ix in range(nx):
		for it in range(nt):
			pass
def Hess_Functional_IPF(X,param):
# TO DO: derive directly the hessian and the gradient of the functional ?
	nu = param.nu
	nx = param.nx
	nt = param.nt
	H = np.zeros(X.size,X.size)
	for ix in range(nx):
		for it in range(nt):
			pass

