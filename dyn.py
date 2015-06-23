# bibliotheque dyn
import numpy as np
import scipy.optimize as opt

def f_lorenz(X,t):

	beta = 8.0/3.0
	rho = 28.0
	sigma = 10.0
	
	dXdt = np.zeros(X.shape)
	dXdt[0] = sigma*(X[1] - X[0])
	dXdt[1] = X[0]*(rho-X[2]) - X[1]
	dXdt[2] = X[0]*X[1] - beta*X[2]
	
	return dXdt

def h_lorenz(X,t):

#	hX = X[1]
	hX = np.abs(np.sum(X))
#	hX = X
	return hX

def f_burgers(u,t):
	# burgers. Upwind scheme + central diff
	nu = 0.01

	dx = 1./u.size
# advection
	up = np.copy(u)
	up[up<0.] = 0.
	um = np.copy(u)
	um[um>0.] = 0.

	dudtp = np.zeros(u.size)
	dudtp[0:-2] = -3.*u[0:-2] + 4.*u[1:-1] - u[2:]
	dudtp[-2] = -3.*u[-2] + 4.*u[-1] - u[0]
	dudtp[-1] = -3.*u[-1] + 4.*u[0] - u[1]

	dudtm = np.zeros(u.size)
	dudtm[2:] = +3.*u[2:] - 4.*u[1:-1] + u[:-2]
	dudtm[1] = +3.*u[1] - 4.*u[0] + u[-1]
	dudtm[0] = +3.*u[0] - 4.*u[-1] + u[-2]

	dudt = -u*(up*dudtm + um*dudtp)/dx

	dudt[1:-1] = dudt[1:-1] + nu * (u[2:] - 2 * u[1:-1] + u[:-2]) / dx**2
	dudt[0] = dudt[0] + nu * (u[1] - 2 * u[0] + u[-1]) / dx**2
	dudt[-1] = dudt[-1] + nu * (u[0] - 2 * u[-1] + u[-2]) / dx**2

	return dudt

def h_burgers(U,t):

	n = np.int_(np.linspace(0,U.size-1,5))
#	hU = U
#	print(np.sum(U))
#	print(np.max(U))
#	print(np.min(U))
#	print(U[n])
	
	hU = np.concatenate((np.array([ np.sum(U),np.max(U),np.min(U)]),U[n] ))
	return hU



def mrand(mean,cov):
	if mean.size == 1:
		monrand = np.random.uniform(mean,cov)
	else:
		monrand = np.random.multivariate_normal(mean,cov)
	#	monrand = np.reshape(monrand,mean.shape)

	return monrand

def kp(f,point,t,dt,g=np.sqrt(2)):
#Klauder perterson scheme
	mean = np.zeros(point.shape)
	cov = g*dt*np.eye(point.size)
	fp = f(point,t)
	pstar = point + dt*fp + mrand(mean,cov)
	fps = f(pstar,t)

	return point+dt*(fps + fp)/2.0 + mrand(mean,cov),pstar,fp,fps

def P_int(f,t,dt,g,X,Xp05,Xp1):
# proba that (Xp05,Xp1) corresponds to transport of X by f, after a time dt.
#Klauder perterson scheme

	fX = f(X,t)
	P1 = norme_vec( Xp05 - X - dt * fX )**2.0 / (2.0 * dt * g * g)

	P2 = norme_vec( Xp1 - X - dt * ( fX + f(Xp05,t) )/2.0 )**2.0 / (2.0 * dt * g * g)

	return P1+P2


def P_int_min(f,t,dt,g,X,Xp1):
# proba that (Xp05,Xp1) corresponds to transport of X by f, after a time dt.
#Klauder perterson scheme

	fX = f(X,t)
	Xs = X + dt * fX
	
	P = norme_vec( Xp1 - X - dt * ( fX + f(Xs,t) )/2.0 )**2.0 / (2.0 * dt * g * g)

	return P

def P_obs(h,t,s,X,y):
# proba that the measure of X corresponds to observation y

	P = norme_vec(h(X,t) - y)**2.0 / (2.0  * s * s)
	return P


def ressample(X,w):
# from M.S. Arulampalam et al., A tutorial on particle filters for online nonlinear / non-Gaussian Bayesian tracking, IEEE Trans. Sign. Proc.
	Ns = w.size
	wrs = np.ones(Ns)/Ns
	cdf = np.zeros(Ns)
	Xrs = X
	for i in range(1,Ns):
		cdf[i] = cdf[i-1] + w[i]
	i=0
	u0 = np.random.uniform(0,1.0/Ns)
	u = np.zeros(Ns)
	perm = np.zeros(Ns)
	for j in range(Ns):
		u[j] = u0 + 1.0*j/Ns
		sortie = False
		while u[j]>cdf[i]:
			if i<Ns-1:
				i=i+1
			else:
				break

		Xrs[:,j] = X[:,i]
		perm[j]=i

	return Xrs,wrs,perm



def Hessian(f,X0,dx=0.1):
	ndim = X0.size
	H = np.zeros((ndim,ndim))
	for i1 in range(0,ndim):
		dX1=np.zeros(ndim)
		dX1[i1] = dx
		for i2 in range(0,i1+1):
			dX2=np.zeros(ndim)
			dX2[i2] = dx
			
			fpp = f(X0+(dX1+dX2)/2.0)
			fpm = f(X0+(dX1-dX2)/2.0)
			fmp = f(X0+(-dX1+dX2)/2.0)
			fmm = f(X0-(dX1+dX2)/2.0)
			d2fdx1dx2 = (fpp-fmp -(fpm-fmm))/(dx*dx)
			H[i1,i2]=d2fdx1dx2
			H[i2,i1]=d2fdx1dx2

	return H

def gen_cholesky(H, eps_machine=1e-15, print_prefix=0, print_flag=0): 
#The Gill, Murray, and Wright modified Cholesky algorithm, from Practical Optimization, Academic Press, London, 1981
 
	ndim = len(H) 
	I = np.eye(ndim) 
	# Calculate gamma(A) and xi(A). 
	gamma = 0.0 
	xi = 0.0 
	for i in range(ndim): 
		gamma = np.max(np.abs(H[i, i]), gamma) 
		for j in range(i+1, ndim): 
			xi = max(np.abs(H[i, j]), xi) 

# Identify delta and beta. 
	try:
		delta = eps_machine * np.max(gamma + xi, 1.0) 
	except:
		print(gamma)
		print(xi)
		print(eps_machine)
		delta = eps_machine * np.max(gamma + xi, 1.0) 

	if ndim == 1: 
		beta = np.sqrt(np.max(gamma, eps_machine)) 
	else: 
		beta = np.sqrt(np.max(gamma, xi / np.sqrt(n**2 - 1.0), eps_machine)) 
 
# Initialise data structures. 
	a = 1.0 * H 
	r = 0.0 * H 
	P = 1.0 * I 

# Main loop. 
	for j in range(ndim): 
# Row and column swapping, find the index > j of the largest diagonal element. 
		q = j 
		for i in range(j+1, ndim): 
			if np.abs(a[i, i]) >= np.abs(a[q, q]): 
				q = i 

# swap row and column j and q (if j != q). 
		if q != j: 
# Temporary permutation matrix for swaping 2 rows or columns. 
			p = 1.0 * I 

# Modify the permutation matrix P by swaping columns. 
			row_P = 1.0*P[:, q] 
			P[:, q] = P[:, j] 
			P[:, j] = row_P 

# Modify the permutation matrix p by swaping rows (same as columns because p = pT). 
			row_p = 1.0*p[q] 
			p[q] = p[j] 
			p[j] = row_p 

# Permute a and r (p = pT). 
		a = np.dot(p, dot(a, p)) 
		r = np.dot(r, p) 
 
# Calculate dj. 
		theta_j = 0.0 
		if j < ndim-1: 
			for i in range(j+1, ndim): 
				theta_j = np.max(theta_j, np.abs(a[j, i])) 
		dj = np.max(np.abs(a[j, j]), (theta_j/beta)**2, delta) 

# Calculate row j of r and update a. 
		r[j, j] = sqrt(dj)     # Damned sqrt introduces roundoff error. 
		for i in range(j+1, ndim): 
			r[j, i] = a[j, i] / r[j, j] 
			for k in range(j+1, i+1): 
				a[i, k] = a[k, i] = a[k, i] - r[j, i] * r[j, k]     # Keep matrix a symmetric. 
 
# Finally, the Cholesky decomposition of H. 
	L = np.dot(P, transpose(r))
	return L

def norme_vec(a):
	return np.sum(a*a)**0.5


class particle_parameters:

	def __init__(self,ptype = False,X_init = np.array([ -4.4090617 ,   0.94099541,  31.65011104]),verbose = False, fdyn = f_lorenz,h_obs = h_lorenz,dt = 0.01,t0=0.0,s_obs = np.sqrt(0.1),g_int=np.sqrt(2.0),objective = 'filter'):
		self.ptype = ptype # true if real state, false is particle
		self.objective = objective # DA for data assimilation, filter for filtering
		self.fdyn = fdyn # dynamical function
		self.h_obs = h_obs # observable function
		self.dt = dt
		self.t = t0

		self.s = s_obs
		self.g = g_int
		self.verbose = verbose
		
		if ptype:
			self.X = X_init
		else:
			self.X = X_init + np.random.uniform(0,10*self.s,X_init.shape)
#			self.X[:] = 0

		self.Y = self.h_obs(self.X,self.t)

	def __len__(self):
		return 1

	def set_type(self,ptype):
		self.ptype = ptype
	def get_type(self):
		return self.ptype

	def set_fdyn(self,fdyn):
		self.fdyn = fdyn
	def get_fdyn(self):
		return self.fdyn

	def set_h(self,h):
		self.h_obs = h_obs
	def get_h(self):
		return self.h_obs

	def set_t(self,t):
		self.t = t
	def get_t(self):
		return self.t

	def set_dt(self,dt):
		self.dt = dt
	def get_dt(self):
		return self.dt

	def set_var_obs(self,s):
		self.s = s
	def get_var_obs(self):
		return self.s

	def set_var_eq(self,g):
		self.g = g
	def get_var_eq(self):
		return self.g

	def set_position(self,X):
		self.X = X
	def get_position(self):
		return self.X	

#	def set_precond(self,cond):
#		self.precond = cond
#	def get_precond(self):
#		return self.precond

	def set_verbose(self,v):
		self.verbose = v
	def get_verbose(self):
		return self.verbose

	def get_dim(self):
		return self.X.size
	def get_dim_obs(self):
		return self.Y.size

	def get_objective(self):
		return self.objective

class particle:

	def __init__(self,param):

		self.ptype = param.get_type()
		self.objective = param.get_objective()
		
		self.ndim = param.get_dim()

		self.g = param.get_var_eq()
		self.s = param.get_var_obs()


		self.fdyn = param.get_fdyn()
		self.h_obs = param.get_h()

		self.X = param.get_position()
		self.Xp1 = np.copy(self.X)
		self.X0 = np.copy(self.X) # for assimilation
		self.Xp05 = np.copy(self.X)

#		self.precond = param.get_precond()

		self.verbose = param.get_verbose()

		self.cov_fact = 1

		self.t = param.get_t()
		self.dt = param.get_dt()

		if self.ptype is True:
			self.path = np.copy(self.X)
			self.Yp1 = self.compute_obs(self.Xp1,self.t)
		else:
			self.old_weight = 1.0
			self.weight = 1.0
			self.eps = np.zeros(2.0*self.ndim)

		self.isready = False # True as soon as a new observation is available
		self.n_obs = 0 # number of time steps between two observations
		self.intermediate_steps = []

	def compute_obs(self,X,t):

		hX = self.h_obs(X,t)
		mean = np.zeros(hX.shape)
		cov = np.eye(hX.size)
		Yp1 = hX + self.s*mrand(mean,cov)

		return Yp1



	def next_step(self):

		self.X = np.copy(self.Xp1)
		self.n_obs = self.n_obs + 1

		if self.ptype is True:
			self.Xp05,self.Xp1 = self.integration(self.X,self.t)
			self.path = np.vstack((self.path,self.Xp1))
			self.Yp1 = self.compute_obs(self.Xp1,self.t)
			
		else:
			if self.objective == 'filter':
				if self.isready is True:
					self.F_min()
					Xp1,self.weight = self.sample()
					self.X0 = np.copy(Xp1[0:self.ndim])
					self.Xp1 = np.copy(Xp1[-self.ndim:])
					if self.n_obs>1:
						self.intermediate_steps = np.array([Xp1[2*i*self.ndim:(2*i+1)*self.ndim] for i in np.array(range(0,self.n_obs))])
					self.isready = False
					self.n_obs = 0
					if self.verbose is True:
						s = 'old point ' + str(self.X) + ' to new point ' + str(self.Xp1)
						print(s)
					self.t = self.t + self.dt
	#			else:
	#				print('to do: store path when estimated')

			elif self.objective == 'DA':
				self.F_min()
				Xp1,self.weight = self.sample()
				self.X0 = np.copy(Xp1[0:self.ndim])
				self.Xp1 = np.copy(Xp1[-self.ndim:])
				self.intermediate_steps = np.array( [ Xp1[i*self.ndim:(i+1)*self.ndim] for i in range(1,Xp1.shape[0]-2) ] )


	def integration(self,X,t):
	#Klauder perterson scheme
		mean = np.zeros(X.shape)
		cov = np.eye(X.size)
		fX = self.fdyn(X,t)
		Xp05 = X + self.dt*fX + self.g*np.sqrt(self.dt)*mrand(mean,cov)
		fXp05 = self.fdyn(Xp05,t)
		Xp1 = X+self.dt*(fXp05 + fX)/2.0 + self.g*np.sqrt(self.dt)*mrand(mean,cov)

		return Xp05,Xp1


	def F(self,X):
		""" Xp05 = X[0:self.ndim]
		Xp1 = X[self.ndim:]
		=> the Klauder perterson scheme is used, and is a two-steps scheme
		"""
		if self.ptype is True:
			print('Reference particle!')
		else:

			if self.objective == 'filter':
				Fint= 0 
				for i in range(self.n_obs):
					if i==0:
#						Fint = Fint + P_int(self.fdyn,self.t,self.dt,self.g,self.X,X[2*i*self.ndim:(2*i+1)*self.ndim],X[(2*i+1)*self.ndim:2*(i+1)*self.ndim])
						Fint = Fint + P_int_min(self.fdyn,self.t,self.dt,self.g,self.X,X[(i+0)*self.ndim:(i+1)*self.ndim])
#						print(X[2*i*self.ndim:(2*i+1)*self.ndim])
#						print(X[(2*i+1)*self.ndim:2*(i+1)*self.ndim])
					else:
						Fint = Fint + P_int(self.fdyn,self.t,self.dt,self.g,X[(2*i-1)*self.ndim:2*(i)*self.ndim],X[2*i*self.ndim:(2*i+1)*self.ndim],X[(2*i+1)*self.ndim:2*(i+1)*self.ndim])

				Fobs = P_obs(self.h_obs,self.t,self.s,X[-self.ndim:],self.Yp1)

				Fx = Fint+Fobs
#				s = 'X = ' + str(X) + ' and registered X is ' + str(self.X)
#				print(s)
#				s = 'F(x) = ' + str(Fx)
#				print(s)
#				s = 'F = ' + str(Fint) + '   '+ str(Fobs)
#				print(s)

			elif self.objective == 'DA':
				Fx = 0
				Fint= 0
				Fobs = 0
				for i in range(self.n_obs.size):

					Fint = Fint + P_int(self.fdyn,self.t[i],self.dt,self.g,X[(2*i-1)*self.ndim:2*(i)*self.ndim],X[2*i*self.ndim:(2*i+1)*self.ndim],X[(2*i+1)*self.ndim:2*(i+1)*self.ndim])

				for i in range(self.n_obs.size):
					iobs = self.n_obs[i]
					Fobs = Fobs + P_obs(self.h_obs,self.t,self.s,X[2*iobs*self.ndim:(2*iobs+1)*iself.ndim],self.Yp1[i])
				Fx = Fx + Fint+Fobs
#		print(Fx)
		return Fx

	def F_ersatz(self,lamda):
		X = self.Xmin + lamda*np.dot(self.L.T , self.eps) / norme_vec(self.eps)
		return np.abs(self.F(X) - self.Fmin - self.EE)

	def F_min(self):

		if self.ptype is True:
			print('Reference particle!')

		else:

			if self.objective == 'filter':
#construction of a good set of initial conditions
				X04min = np.zeros(2*self.ndim*self.n_obs)
				for i in range(self.n_obs):
					if i == 0:
						Xp05,Xp1 = self.integration(self.X,self.t-(self.n_obs-i-1)*self.dt)
					else:
						Xp05,Xp1 = self.integration(X04min[(2*(i-1))*self.ndim:(2*(i-1)+1)*self.ndim],self.t-(self.n_obs-i-1)*self.dt)
					X04min[2*i*self.ndim:(2*i+1)*self.ndim] = Xp05
					X04min[(2*i+1)*self.ndim:2*(i+1)*self.ndim] = Xp1
				self.debug = X04min

# minimisation
#				res = opt.minimize(self.F, X04min, method='BFGS',options={'gtol': 1e-4,'disp':self.verbose})
#				print(X04min)
				res = opt.minimize(self.F, X04min, method='BFGS')
				s = 'F = ' + str(self.F(res.x)) 
				print(s)

				
# if the mimimisation have failed
				if res.success is False:
					print('minimisation failed')
					print('try with ncg algorithm')
					print(res.message)
					self.cov_fact = 15
					res = opt.minimize(self.F, X04min, method='Nelder-Mead',jac=None, hess=None)
					if res.success is False:
						print('minimisation failed')
						print('try with cg algorithm')
						print(res.message)
						res = opt.minimize(self.F, X04min,jac=None, hess=None,method='CG')
						if res.success is False:
							print('minimisation failed again')
							print('Point will just move with the flow')
							print(res.message)

	#						Xp05,Xp1 = self.integration(self.X,self.t)
							res.x =  X04min
							res.fun = self.F(res.x)
				else:
					self.cov_fact = 1

				self.res = res
				self.Xmin = res.x
				self.Fmin = res.fun

				if res.success is True:
					self.H = Hessian(self.F,self.Xmin, 0.1)
	#				print(np.linalg.cond(self.H))
					try:
						self.L = np.linalg.cholesky(np.linalg.pinv(self.H))
					except np.linalg.LinAlgError:
# if the minimum is not that good, the (inverse of the) hessian might not be positive: we use a general cholesky decomposition:
						print('use of gmw')
						try:
							self.L = gen_cholesky(np.linalg.pinv(self.H))
						except:
							print('Singular matrix')
							self.H = np.eye(2.0*self.ndim)
							self.L = self.H

				else: # the minimum is badly estimated anyway
					self.H = np.eye(2.0*self.ndim)
					self.L = self.H


			elif self.objective == 'DA':
				X04min = np.concatenate( ( self.X0.flatten(),self.intermediate_steps.flatten(),self.Xp1.flatten() ) )

				res = opt.minimize(self.F, X04min, method='BFGS',options={'gtol': 1e-4,'disp':self.verbose})

				if res.success is False:
					print('minimisation failed')
					print('try with ncg algorithm')
					print(res.message)
					self.cov_fact = 15
					res = opt.minimize(self.F, X04min, method='Nelder-Mead',jac=None, hess=None)
					if res.success is False:
						print('minimisation failed')
						print('try with cg algorithm')
						print(res.message)
						res = opt.minimize(self.F, X04min,jac=None, hess=None,method='CG')
						if res.success is False:
							print('minimisation failed again')
							print('Point will just move with the flow')
							print(res.message)

							res.x =  X04min
							res.fun = self.F(res.x)
				else:
					self.cov_fact = 1

				self.res = res
				self.Xmin = res.x
				self.Fmin = res.fun
				if res.success is True:
					self.H = Hessian(self.F,self.Xmin, 0.1)
					try:
						self.L = np.linalg.cholesky(np.linalg.pinv(self.H))
					except np.linalg.LinAlgError:
						print('use of gmw')
						self.L = gen_cholesky(np.linalg.pinv(self.H))
				else:
					self.H = np.eye(2.0*self.ndim)
					self.L = self.H


	def sample(self):
		if self.ptype is True:
			print('Reference particle!')
		else:
			mean = np.zeros(self.Xmin.size)
			
#			if self.res.success is True:
#				self.cov_fact = 1.0
#			else: # need to seek away
#				self.cov_fact = self.cov_fact + 5 
#			cov = self.cov_fact*np.eye(self.Xp1.size+self.Xp05.size) 

			weightisok = False
			nit = 0
			cov = self.cov_fact*np.eye(self.Xmin.size) 
			
			while weightisok is False:
				
				self.eps = mrand(mean,cov)
				self.EE = 0.5 * norme_vec(self.eps)**2.0
				lamda = self.solveEqAlgebraic()
				Xp1 = self.Xmin + lamda*np.dot(self.L.T , self.eps) / norme_vec(self.eps)
				weight = self.Id_Weight(Xp1,lamda)
				nit = nit + 1
				if nit>2:
					weightisok = True
				if weight > 0:
					weightisok = True
#			self.cov_fact = 15.0

#			Xp1 = Xp1[self.ndim:]
#			Xp05 = Xp1[0:self.ndim]

			if self.verbose is True:
				s = 'Algebraic solution lambda: ' + str(lamda)
				print(s)
				s = 'Xmin: ' + str(self.Xmin[self.ndim:]) + ' to sample: ' + str(Xp1)
				print(s)

			return Xp1,weight

	def solveEqAlgebraic(self):
		if self.ptype is True:
			print('Reference particle!')
		else:
			rho = norme_vec(self.eps)
			res = opt.minimize(self.F_ersatz, np.sqrt(rho), method='BFGS',options={'gtol': 1e-4,'disp':self.verbose})
			lamda = res.x

			return lamda

	def Id_Weight(self,X,lamda):
		if self.ptype is True:
			print('Reference particle!')
		else:
			rho = norme_vec(self.eps)
			#numeric estimation of d lamda / d rho
			dl = 0.01
			Fp = self.F(self.Xmin + (lamda+dl) * np.dot(self.L.T , self.eps) / norme_vec(self.eps))
			Fm = self.F(X)
			dldp = dl/(2.0*(Fp-Fm))
			# jacobian of the map
			J = 2.0 * np.abs(np.linalg.det(self.L)) * (rho ** (1 - self.ndim)) * np.abs( (lamda**(self.ndim-1.0)) * dldp)
#			s = 'J: ' + str(np.abs(np.linalg.det(self.L))) + '  ' + str(rho) + '    ' + str(np.abs( (lamda**(2.0*self.ndim-1.0)))) + '  ' + str(dldp)
#			print(s)
			#weight
			weight = self.weight * np.exp(-self.Fmin) * J
#			s = 'wn: ' + str(self.weight) + ' expphi: ' + str(np.exp(-self.Fmin)) + ' J: ' + str(J)
#			print(s)

			if self.verbose is True:
				if np.isnan(weight):
					print('isnan')
					print(J)
					print(self.Fmin)
				if weight == 0:
					print('iszeros')
					print(J)
					print(self.Fmin)
			if np.isnan(weight):
				weight = 0
			return weight

	def set_weight(self,w):
		self.old_weight = self.weight
		self.weight = w

	def get_weight(self):
		if self.ptype is True:
			print('Reference particle!')
		else:
			return self.weight
	def get_old_weight(self):
		if self.ptype is True:
			print('Reference particle!')
		else:
			return self.old_weight

	def get_last_position(self):
		return self.X
	def set_last_position(self,X):
		self.X = X

	def get_init_position(self):
		return self.X0

	def get_current_position(self):
		return self.Xp1
	def set_current_position(self,Xp1):
		self.Xp1 = Xp1

	def get_intermediate_steps(self):
		return self.intermediate_steps

	def temporary_set_var(self,var):
		self.cov_fact = var

	def set_obs(self,Y):
		self.Yp1 = Y
		if self.objective == 'filter':
			self.isready = True

	def set_obs_time(self,nobs):
		self.n_obs = nobs

	def get_obs(self):
		return self.Yp1

	def set_init(self,XO):
		n = X0.size/self.ndim
		self.X0 = X0[0:self.ndim]
		self.intermediate_steps = np.array([X0[i*self.ndim:(i+1)*self.ndim] for i in np.array(range(1,n-1))])


class density_particle:
	def __init__(self,n_particle = 15, param_ref = particle_parameters(True),param_p = particle_parameters(False),isverbose = False):
		self.new_data = True
		self.isverbose = isverbose

		self.n_particle = n_particle
		self.p_ref = particle(param_ref)
		self.current_estimate_position = np.zeros(param_ref.get_dim())
		if len(param_p) == self.n_particle:
			self.liste_p = [particle(param_p[k]) for k in range(n_particle)]
		else:
			self.liste_p = [particle(param_p) for k in range(n_particle)]

		for i_p in range(self.n_particle):
			self.current_estimate_position = self.current_estimate_position + self.liste_p[i_p].get_current_position()
		self.current_estimate_position = self.current_estimate_position/n_particle

		self.X_ips = np.zeros((param_ref.get_dim(),n_particle))
		self.intermediate_steps = []

		for i_p in range(self.n_particle):
			self.X_ips[:,i_p] = self.liste_p[i_p].get_current_position()

		self.w = np.zeros(n_particle)
		self.n_obs = []

	def __len__(self):
		return self.n_particle

	def next_step(self):
	#it of main particle
		self.p_ref.next_step()
		obs = self.p_ref.get_obs()
		self.current_estimate_position[:] = 0

		if self.new_data is True:
			print('give obs')
# next step
		for i_p in range(0,self.n_particle):
#			print(i_p)
			if self.new_data is True:
				self.liste_p[i_p].set_obs(obs)
			self.liste_p[i_p].next_step()
			self.X_ips[:,i_p] = self.liste_p[i_p].get_current_position()
			if self.isverbose:
				s = 'particle no: ' + str(i_p)
				print(s)
			self.w[i_p] = self.liste_p[i_p].get_weight()
#		print(self.w)
# normalize weights
		sumw = np.sum(self.w)
		if sumw==0.0:
			print('lost')
			print(self.w)
			self.w=np.ones(self.w.size)/self.n_particle
		else:
			self.w = self.w/sumw

		for i_p in range(0,self.n_particle):
			self.liste_p[i_p].set_weight(self.w[i_p])
#		print(self.w)

# approximate position
		for i_p in range(0,self.n_particle):
			self.current_estimate_position = self.current_estimate_position + self.liste_p[i_p].get_current_position() * self.w[i_p]
			
			if sumw ==0:
				self.liste_p[i_p].temporary_set_var(15)

		try:
			self.intermediate_steps = self.liste_p[0].get_intermediate_steps() * self.w[0]
		except:
			self.intermediate_steps = self.liste_p[0].get_intermediate_steps() * 0.0
		for i_p in range(1,self.n_particle):
			try:
				self.intermediate_steps = self.intermediate_steps + self.liste_p[i_p].get_intermediate_steps() * self.w[i_p]
			except:
				self.intermediate_steps = self.intermediate_steps + self.liste_p[i_p].get_intermediate_steps() * 0.0
			
#		print(w)


#ressample	
		doressample = False
		if np.random.uniform()<0.05: # force ressample from time to time
			doressample = True
		if self.w.min() < 1e-8: # one particle is neglectible
			doressample = True
		if self.w.max() > 5.0/self.n_particle: # one particle is too important
			doressample = True

		if doressample is True:
			Xrs,wrs,perm = ressample(self.X_ips,self.w)
			if self.isverbose:
				s = "resampling: " + str(np.unique(perm))
				print(s)
#			s = "resampling: " + str(np.unique(perm)) + ' weights: ' + str(self.w)
#			print(s)
			for i_p in range(0,self.n_particle):

				self.liste_p[i_p].set_weight(wrs[i_p]/np.sum(wrs))
				self.liste_p[i_p].set_current_position(Xrs[:,i_p])

	def new_data_provided(self,boolean):
		self.new_data = boolean

	def get_p(self,i):
		return self.liste_p[i]
	def get_ref(self):
		return self.p_ref
	def get_current_position(self):
		return self.p_ref.get_current_position()

	def get_estimate_position(self):
		return self.current_estimate_position
	def get_current_positions(self):
	# position BEFORE ressampling
		return self.X_ips

	def get_estimate_intermediate_position(self):
		return self.intermediate_steps

	def get_estimate_intermediate_positions(self,i_p):
		return self.liste_p[i_p].get_intermediate_steps()
		
	def set_obs_and_time(self,Y,n_obs):
		if Y.size == n_obs.size:
			for i in range(self.n_particle):
				self.liste_p[i].set_time_obs(n_obs)
				self.liste_p[i].set_obs(Y)
		else:
			print('Error: size between observables and times do not match')


