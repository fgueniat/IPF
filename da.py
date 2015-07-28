import numpy as np
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b
from dyn import *

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
	nt = param.nt
	dt = param.dt

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
	return -1.0*dJ # WHY IT IS MUCH MORE EFFICIENT WITH A x10 ??

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



