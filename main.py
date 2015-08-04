import numpy as np
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b

import plot_tools as pt
from la_tools import *
from dyn import *
from da import *
from ipf import *


class Parameters():
	def __init__(self,dx,dt,nx,nt,nu,psize = 0, obs=False,g=0,s=0,imp_int = False,imp_obs = False,eps = np.eye(2),fmin = 0, hessianm1 = np.eye(2),xmin = np.ones(2)):
		self.dx = dx
		self.dt = dt
		self.nx = nx
		self.nt = nt
		self.m = 10
		self.obs = obs
		self.psize = psize
		self.g = g
		self.s = s
		self.fmin = fmin
		self.xmin = xmin
		if imp_int is False:
			self.imp_int = 1.0*np.linspace(10*nt,1.0,nt)
		else:
			self.imp_int = np.ones(nt)
		if imp_obs is False:
			self.imp_obs = 5.0*np.linspace(10*nt,1.0,nt)
		else:
			self.imp_obs = np.ones(nt)
		self.Set_eps(eps)
		self.Set_hess(hessianm1)

	def Set_eps(self,eps):
		self.eps = eps
		self.epseps2 = 0.5 * np.dot(eps.T,eps)
		try:
			self.Lnu = np.dot(self.L.T,self.eps)/np.dot(eps.T,eps)
		except:
			print "L not defined"

	def Set_hess(self,Hm1):
		self.hessianm1 = Hm1
		try:
			self.L = np.linalg.cholesky(self.hessianm1)
		except:
			print("gen cholesky")
			self.L = gen_cholesky(self.hessianm1)

	def GG(m):
		

if __name__ == '__main__':
	print "This only executes when %s is executed rather than imported" % __file__
	
	ntest = 10
	ng = 1
	ns = 1
	err_da = np.zeros((ntest,2))
	err_f = np.zeros((ntest,2))
	lg = np.logspace(-5,-2,ns)
	ls = np.logspace(-5,-2,ns)
	for ig in range(lg.size):
		for i_s in range(ls.size):
			for itest in range(ntest):
				s = "loop: "+str(itest) + " g: "+str(lg[ig]) + " s: "+str(ls[i_s])
				print(s)

				print("Definition of parameters")
				nx = 50
				nu = 0.0001
				xl = 0.0
				xr = 1.0
				x = np.linspace(xl,xr,nx)
				dx = (xr-xl)/nx
				dt = 0.05
				nt = 3
				g = lg[ig]
				s = ls[i_s]
				T = dt*nt
				param = Parameters(dx,dt,nx,nt,nu,nx,[],g,s)

				print("            ***               ")
				print("            ***               ")
				print("Definition of the real initial conditions")
				u0 = 0.2*( 1.05 + -np.sin( -np.pi + 2.0*np.pi*(x-xl)/(xr-xl) ) ) + 0.01*np.random.normal(0,1,x.size)
			#	u0 = 0.2*( 1.05 + -np.sin( -np.pi + 2.0*np.pi*(x-xl)/(xr-xl) ) )


				print("            ***               ")
				print("Integration of the system, and computations of the observables")
				U = np.zeros((nx,nt))
				U[:,0] = u0
				y0 = Observable(u0,param)
				Y = np.zeros((y0.size,nt))
				Y[:,0] = y0

				unp1 = np.copy(U[:,0])
				for i in range(1,nt):
					un = np.copy(unp1)
					unp1 = Newton_step(un,param)
				#	if np.mod(i,1) == 0:
				#		Plot(unp1,un)
					U[:,i] = unp1
					Y[:,i] = Observable(unp1,param)

				param.obs = Y


				print("            ***               ")
				print("            ***               ")
				print("Approximate initial conditions")
			#	p0 = 0.2*( 1.0 + -np.sin( -np.pi + 2.05*np.pi*(x-xl)/(xr-xl) ) ) + 0.5*np.random.normal(0,0.1,x.size) 
				p0 = np.mean(u0) + 0.2*np.random.normal(0,0.1,x.size) 
				p0 = np.mean(u0) * np.ones(x.size) 
			#	p0 = 0.2*( 1.05 + -np.sin( -np.pi + 2.0*np.pi*(x-xl)/(xr-xl) ) ) + 0.2*np.random.normal(0,0.1,x.size)



			#	# brute force
			#	print("Minimisation over all the parameters")
			#	ic_opt = fmin_bfgs(Functionnal_p,p0,None,(param,))
				print("            ***               ")
				print("Adjoint/4D-VAR")
			#	ic_opt = fmin_l_bfgs_b(Functionnal_p,p0,Grad_Functional,(param,),iprint=1,factr = 0.0)
				res = fmin_bfgs(Functionnal_p,p0,Grad_Functional,(param,),full_output = 1, disp = 0)
				ic_4DVAR = res[0]
				param.umin = ic_4DVAR
			#	Plot(u0,p0)
			#	Plot(u0,p0)
			#	Plot(u0,ic_4DVAR)
				#Plot(u0,p0)
			#	import pdb #@@@
			#	pdb.set_trace() #@@@


				print("            ***               ")
				print("            ***               ")
				print("Start from 4D-VAR for the minimisation of probabilities")
				Usmart = np.zeros((nx,nt))
				Usmart[:,0] = ic_4DVAR
				unp1 = np.copy(Usmart[:,0])

				for i in range(1,nt):
					un = np.copy(unp1)
					unp1 = Newton_step(un,param)
				#	if np.mod(i,1) == 0:
				#		Plot(unp1,un)
					Usmart[:,i] = unp1

				print("            ***               ")
				print("Minimisation with weigths on first times")
			#	p02 = U.flatten() + .01*np.random.normal(0,1,U.size)
				p02 = Usmart.flatten()
				ic_smart = fmin_l_bfgs_b(Functionnal_IPF,p02,None,(param,),approx_grad=True ,iprint=0,factr = 0.0)
				param.imp_int = np.ones(nt)
				param.imp_obs = np.ones(nt)

				print("            ***               ")
				print("Regular Minimisation")
			#	ic_opt2 = fmin_l_bfgs_b(Functionnal_IPF,ic_smart[0],None,(param,),approx_grad=True ,iprint=1,factr = 0.0)
				res = fmin_bfgs(Functionnal_IPF,ic_smart[0],None,(param,),full_output = 1, disp = 0)
				ic_opt_IPF = res[0]
				param.fmin = res[1]
				param.xmin = ic_opt_IPF
				print("            ***               ")
				print("Hessian")
				hesstype = "approx"
				if hesstype == "bfgs":
					hess = res[3]
				else:
					hess = Hessianm1(ic_opt_IPF,param,0.01)

				param.Set_hess(hess)

				ic_opt_IPF = ic_opt_IPF.reshape(param.nx,param.nt)
			#	Plot(u0,ic_opt_IPF[:,0])


				print("            ***               ")
				print("            ***               ")
				print("Particules")
				n_particle = 15
				ICs = np.zeros((param.nx,n_particle))
				TCs = np.zeros((param.nx,n_particle))
				w = np.zeros(n_particle)
				for ip in range(n_particle):
					print("            ***               ")
					s = "Particules no: " + str(ip)
					print(s)
					eps = np.random.normal(0,1,ic_opt_IPF.size)
					param.Set_eps(eps)
					rho = np.sqrt( np.dot(eps.T,eps) )
					res = fmin_bfgs(Functionnal_IPF_lin,rho,None,(param,),full_output = 1, disp = 1)
					ui = Ersatz(res[0],param)
					ICs[:,ip] = ui[:,0]
					TCs[:,ip] = ui[:,-1]
					w[ip] = Id_Weight(res[0],param)
				w2 = Weight_rescale(w)
				u0_IPF = np.zeros(nx)
				uf_IPF = np.zeros(nx)
				for ip in range(n_particle):
					u0_IPF += w2[ip] * ICs[:,ip]
					uf_IPF += w2[ip] * TCs[:,ip]


				U4dvar = np.zeros((nx,nt))
				U4dvar[:,0] = ic_4DVAR

				unp1 = np.copy(U4dvar[:,0])
				for i in range(1,nt):
					un = np.copy(unp1)
					unp1 = Newton_step(un,param)
				#	if np.mod(i,1) == 0:
				#		Plot(unp1,un)
					U4dvar[:,i] = unp1

				s = "diff u0 u4dvar: " + str(np.linalg.norm(u0 - ic_4DVAR)/np.linalg.norm(u0)) + " & diff u0 IPF: " + str(np.linalg.norm(u0 - u0_IPF)/np.linalg.norm(u0))
				print(s)
				s = "diff u0 u4dvar: " + str(np.linalg.norm(U[:,-1] - U4dvar[:,-1])/np.linalg.norm(u0)) + " & diff u0 IPF: " + str(np.linalg.norm(U[:,-1] - uf_IPF)/np.linalg.norm(u0))
				print(s)
				err_da[itest,0] =  np.linalg.norm(u0 - ic_4DVAR)/np.linalg.norm(u0)
				err_da[itest,1] =  np.linalg.norm(u0 - u0_IPF)/np.linalg.norm(u0)
		
				err_f[itest,0] =  np.linalg.norm(U[:,-1] - U4dvar[:,-1])/np.linalg.norm(u0)
				err_f[itest,1] =  np.linalg.norm(U[:,-1] - uf_IPF)/np.linalg.norm(u0)
				pt.multiplot1(False,(u0,ic_4DVAR))
				pt.multiplot1(False,(u0,u0_IPF))

#				s = '../' + 'loop_nx_' + str(nx) + '_g_' + str(ig) + '_s_' + str(i_s) + '.dat'
#				data = (u0, U , param ,err_da ,err_f , u0_IPF, ic_4DVAR )
#				PIK = s
#				with open(PIK, "wb") as f:
#					pk.dump(len(data), f)
#					for value in data:
#						pk.dump(value, f)
				#load:
#				data2 = []
#				with open(PIK, "rb") as f:
#					for _ in range(pk.load(f)):
#						data2.append(pk.load(f))


