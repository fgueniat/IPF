import numpy as np
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b
from dyn import *
from la_tools import *


def P_int(f,X,Xp1,param):
# proba that (Xp05,Xp1) corresponds to transport of X by f, after a time dt.
	g=param.g
	dt = param.dt
	if np.abs(g) < 1.0e-12:
		g=1.0
	fX = f(X,param)
	P = mdot((Xp1-X - dt*fX, param.Gm1,Xp1[:,np.newaxis]-X[:,np.newaxis] - dt*fX[:,np.newaxis]))
#	P = np.linalg.norm( Xp1 - X - dt * f(X,param) )**2.0 / (2.0 * dt * g * g)
#	print(P)
	return P
#	return 0.0

def P_obs(h,X,y,param):
# proba that the measure of X corresponds to observation y
	s=param.s
	if np.abs(s) < 1.0e-12:
		s=1.0
#	P = np.linalg.norm(h(X,param) - y)**2.0 / (2.0  * s * s)
	hX = h(X,param)
	P = mdot((hX - y, param.Sm1,hX[:,np.newaxis]-y[:,np.newaxis]))
	return P

def F_intermediary(u,param):
# cancel noise for that part
	g = param.g
	param.is_g = 0.0
	unp1 = Newton_step(u,param)
	param.is_g = 1.0
	return (unp1-u)/param.dt

def Observable_intermediary(u,param):
# cancel noise for that part
	s = param.s
	param.is_s = 0.0
	y = Observable(u,param)
	param.is_s = 1.0
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

def Jacob_Functionnal_IPF(p,param):
# Jacob of Functional to minimize
	nt = param.nt
	nx = param.nx
	F = np.zeros((nx,nt))
	h = Observable_intermediary
	f = F_intermediary

	
	X = p.reshape(param.nx,param.nt)
	
	for i in range(nt-1): #DA: t=0 included
		#J = 2*(df/dxn+1 + df/dx) * (GG^T)-1 * f
		JFGG = np.dot(JF_int(X[:,i],param), param.GG) * (2*param.dt * param.g * param.g)
		F[:,i] = 2.0 * np.dot( JFGG, F_int(f,X[:,i],X[:,i+1],param) )
		JFSS = np.dot( JF_obs(X[:,i],param) , param.SS) * (2* param.s * param.s)
		F[:,i] += 2.0 * np.dot( JFSS,  F_obs(h,X[:,i],param.obs[:,i],param)   )
	JFSS = np.dot( JF_obs(X[:,-1],param) , param.SS) * (2* param.s * param.s)
	F[:,-1] += 2.0 * np.dot( JFSS,  F_obs(h,X[:,-1],param.obs[:,-1],param)   )
	F = F.flatten()
#	print(np.max(F))
	return 1.0 * F

def F_int(f,X,Xp1,param):
#Xnp1-f(Xn)
	g = param.g
	if np.abs(g) < 1.0e-12:
		g=1.0
	dt = param.dt
	F = Xp1 - X - dt * f(X,param)
#	print(P)
	return F


def JF_int(X,param):
# D(Xnp1-f(Xn))
	g = param.g
	if np.abs(g) < 1.0e-12:
		g=1.0
	dt = param.dt
	dx = param.dx
	nu = param.nu
	xm1 = np.append(X[-1],X[0:-1])
	JF = np.diag( 2.0*X - xm1 ) / dx - nu * (np.diag( np.ones(  X.size-1  ),1 ) - 2.0*np.diag( np.ones(  X.size  ) + np.diag(  np.ones(   X.size-1   ),-1  ) )) / (dx*dx)
	JF[0,-1] = -nu/(dx * dx)
	JF[-1,0] = -nu/(dx * dx)
	JF = np.zeros((X.size,X.size))
#	JF += np.eye(Xp1.size) - np.eye(X.size) - dt * JF
	JF = - dt * JF
#	print(P)1
	return JF

def F_obs(h,X,y,param):
	g = param.g
	if np.abs(g) < 1.0e-12:
		g=1.0
	dt = param.dt
	F = h(X,param) - y
#	print(P)
	return F


def JF_obs(X,param):
	JF = np.eye(X.size)
#	print(P)1
	return JF

def Hessian_IPF(p,param):
	H = np.zeros(2)
	print('approx with Jacob ?')
	pass
	
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

def Hess_deriv(X,param):
	pass



def Hess_gradiant(X,param):
# TO DO: derive directly the hessian and the gradient of the functional ?
	nu = param.nu
	nx = param.nx
	nt = param.nt
	H = np.eye(X.size)
#	H = np.zeros((X.size,X.size))
	Id = np.eye(X.size)
	Xp = X.flatten()  #+ 0.1*np.random.normal(0,1,X.size)
	Jp = -1.0*Jacob_Functionnal_IPF(X,param)
	for ih in range(100):
		nJ = np.linalg.norm(Jp)
		s = 'log func: ' +str(np.log(Functionnal_IPF(Xp,param))) + ' norm grad: ' + str(nJ)
		print(s)
		if nJ<1e-8:
			break
		X = 1.0*Xp
		J = 1.0*Jp
		if ih ==0:
			dX = -1.0*0.001*J
		else:
			dX = -1.0*0.000001*mdot((H,J))
		Xp = X + dX
		Jp = -1.0*Jacob_Functionnal_IPF(Xp,param)
		Y = Jp-J

		H = H + mdot((  dX[:,np.newaxis] - mdot((H,Y))[:,np.newaxis]  ,  dX[np.newaxis,:] - mdot((H,Y))[np.newaxis,:]  )) / mdot((  dX - mdot((H,Y))  ,  Y  ))
#		print( mdot((  temp ,  Y  )))
	return H,dX,Y



