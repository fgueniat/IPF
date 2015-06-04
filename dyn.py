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
		gamma = max(abs(H[i, i]), gamma) 
		for j in range(i+1, ndim): 
			xi = max(abs(H[i, j]), xi) 

# Identify delta and beta. 
	delta = eps_machine * max(gamma + xi, 1.0) 
	if n == 1: 
		beta = np.sqrt(np.max(gamma, eps_machine)) 
	else: 
		beta = np.sqrt(np.max(gamma, xi / sqrt(n**2 - 1.0), eps_machine)) 
 
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

	def __init__(self,ptype = False,X_init = np.array([ -4.4090617 ,   0.94099541,  31.65011104]),verbose = False,precond = np.array([15.0,15.0,25.0]), fdyn = f_lorenz,h_obs = h_lorenz,dt = 0.01,t0=0.0,s_obs = np.sqrt(0.1),g_int=np.sqrt(2.0)):
		self.ptype = ptype # true if real state, false is particle
		self.fdyn = fdyn # dynamical function
		self.h_obs = h_lorenz # observable function
		self.dt = dt
		self.t = t0
		self.s = s_obs
		self.g = g_int
		self.precond = np.diag(np.ones(precond.size)/precond)
		self.verbose = verbose
		if ptype:
			self.X = X_init
		else:
			self.X = X_init + np.random.uniform(0,1,X_init.shape)

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

	def set_precond(self,cond):
		self.precond = cond
	def get_precond(self):
		return self.precond

	def set_verbose(self,v):
		self.verbose = v
	def get_verbose(self):
		return self.verbose

	def get_dim(self):
		return self.X.size


class particle:

	def __init__(self,param):
		
		self.ptype = param.get_type()

		self.g = param.get_var_eq()
		self.s = param.get_var_obs()

		self.ndim = param.get_dim()
		self.fdyn = param.get_fdyn()
		self.h_obs = param.get_h()
		self.X = param.get_position()
		self.Xp1 = self.X
		self.Xp05 = self.X
		self.precond = param.get_precond()
		self.verbose = param.get_verbose()

		self.cov_fact = 1
		self.t = param.get_t()
		self.dt = param.get_dt()

		if self.ptype is True:
			self.path = self.X
			self.Yp1 = self.compute_obs(self.Xp1,self.t)
		else:
			self.weight = 1.0
			self.eps = np.zeros(2.0*self.ndim)


#		self.verbose = True

	def compute_obs(self,X,t):

		hX = self.h_obs(X,t)
		mean = np.zeros(hX.shape)
		cov = self.s*np.eye(hX.size)
		Yp1 = hX + mrand(mean,cov)

		return Yp1



	def next_step(self):

		self.X = self.Xp1

		if self.ptype is True:
			self.Xp05,self.Xp1 = self.integration(self.X,self.t)
			self.path = np.vstack((self.path,self.Xp1))
			self.Yp1 = self.compute_obs(self.Xp1,self.t)
		else:
			self.F_min()
			self.Xp05,self.Xp1,self.weight = self.sample()
		self.t = self.t + self.dt
		if self.verbose is True:
			s = 'old point ' + str(self.X) + ' to new point ' + str(self.Xp1)
			print(s)


	def integration(self,X,t):
	#Klauder perterson scheme
		mean = np.zeros(X.shape)
		cov = self.g*self.dt*np.eye(X.size)
		fX = self.fdyn(X,t)
		Xp05 = X + self.dt*fX + mrand(mean,cov)
		fXp05 = self.fdyn(Xp05,t)
		Xp1 = X+self.dt*(fXp05 + fX)/2.0 + mrand(mean,cov)

		return Xp05,Xp1

	
	def F(self,X):
		""" Xp05 = X[0:self.ndim]
		Xp1 = X[self.ndim:]
		=> the Klauder perterson scheme is used, and is a two-steps scheme
		"""
		if self.ptype is True:
			print('Reference particle!')
		else:

#			Xn = self.X
#			t = self.t
#			dt = self.dt

#			Xp05 = X[0:self.ndim]
#			Xp1 = X[self.ndim:]

			

#			fX = self.fdyn(Xn,t)
##			F1 = norme_vec( np.dot(self.precond,Xp05 - Xn - dt * fX) )**2.0 / (2.0 * dt * self.g * self.g)
#			F1 = norme_vec( Xp05 - Xn - dt * fX )**2.0 / (2.0 * dt * self.g * self.g)
#			fXs = self.fdyn(Xp05,t)
##			F2 = norme_vec( np.dot(self.precond,Xp1 - Xn - dt * (fX + fXs)/2.0) )**2.0 / (2.0 * dt * self.g * self.g)
#			F2 = norme_vec( Xp1 - Xn - dt * (fX + fXs)/2.0 )**2.0 / (2.0 * dt * self.g * self.g)
#			F3 = norme_vec(self.h_obs(Xp1,t) - self.Yp1)**2.0 / (2.0 * dt * self.s * self.s)
			Fint = P_int(self.fdyn,self.t,self.dt,self.g,self.X,X[0:self.ndim],Xp1 = X[self.ndim:])
			Fobs = P_obs(self.h_obs,self.t,self.s,X[self.ndim:],self.Yp1)
			Fx = Fint+Fobs
			s = 'F: ' + str(Fint) + '   ' + str(Fobs)
#			print(s)
#			if self.verbose is True:
#				s = 'X: ' + str(np.around(X,decimals=3)) + ' F(X): ' + str(np.around(Fx,decimals=3)) + ' F1: ' + str(np.around(F1,decimals=3)) + ' F2: ' + str(np.around(F2,decimals=3)) + ' F3: ' + str(np.around(F3,decimals=3))
#				print(s)

			return Fx

	def F_ersatz(self,lamda):
		X = self.Xmin + lamda*np.dot(self.L.T , self.eps) / norme_vec(self.eps)
		return np.abs(self.F(X) - self.Fmin - self.EE)

	def F_min(self):
		if self.ptype is True:
			print('Reference particle!')
		else:
#			res = opt.fmin_bfgs(self.F,np.concatenate((self.X,self.X)),full_output=True)
#			self.Xmin = res[0]
#			self.Fmin = res[1]
#			self.H = res[3] # inverse of the Hessian at the min point_init
#			# usefull for the estimation of L, hence the weight
#			if res[6]==2: #bfgs failed ?
#				res = opt.minimize(self.F, self.Xmin, method = 'Nelder-Mead')
#				self.Xmin = res.x
#				self.Fmin = res.fun
#				self.H = Hessian(self.F,self.Xmin, 0.01*np.abs(self.Xmin).min())

#			self.L = np.linalg.cholesky(self.H)
			dX = self.fdyn(self.X,self.t)
			X0 = np.concatenate((self.X+dX*self.dt/2,self.X+dX*self.dt))
#			print('x0 et res')
#			print(self.X)
#			print(X0)
			res = opt.minimize(self.F, X0, method='BFGS',options={'gtol': 1e-4,'disp':self.verbose})
#			print(res.x)
			if res.success is False:
				print('minimisation failed')
				print('try with ncg algorithm')
				print(res.message)
				self.cov_fact = 15
				res = opt.minimize(self.F, X0, method='Nelder-Mead',jac=None, hess=None)
				if res.success is False:
					print('minimisation failed')
					print('try with cg algorithm')
					print(res.message)
					res = opt.minimize(self.F, X0,jac=None, hess=None,method='CG')
					if res.success is False:
						print('minimisation failed again')
						print('Point will just move with the flow')
						print(res.message)
						Xp05,Xp1 = self.integration(self.X,self.t)
						res.x =  np.concatenate((Xp05,Xp1))
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
					self.L = np.linalg.cholesky(np.linalg.inv(self.H))
				except np.linalg.LinAlgError:
					print('use of gmw')
					self.L = gen_cholesky(np.linalg.inv(self.H))
#				print(self.H)
#				print(self.L)
			else:
				self.H = np.eye(2.0*self.ndim)
				self.L = self.H


	def sample(self):
		if self.ptype is True:
			print('Reference particle!')
		else:
			mean = np.zeros(self.Xp1.size+self.Xp05.size)
			
#			if self.res.success is True:
#				self.cov_fact = 1.0
#			else: # need to seek away
#				self.cov_fact = self.cov_fact + 5 
#			cov = self.cov_fact*np.eye(self.Xp1.size+self.Xp05.size) 

			weightisok = False
			nit = 0
			cov = self.cov_fact*np.eye(self.Xp1.size+self.Xp05.size) 
			
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

			Xp1 = Xp1[self.ndim:]
			Xp05 = Xp1[0:self.ndim]

			if self.verbose is True:
				s = 'Algebraic solution lambda: ' + str(lamda)
				print(s)
				s = 'Xmin: ' + str(self.Xmin[self.ndim:]) + ' to sample: ' + str(Xp1)
				print(s)

			return Xp05,Xp1,weight

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
			return weight

	def set_weight(self,w):
		self.weight = w
	def get_weight(self):
		if self.ptype is True:
			print('Reference particle!')
		else:
			return self.weight

	def get_last_position(self):
		return self.X
	def set_last_position(self,X):
		self.X = X

	def get_current_position(self):
		return self.Xp1
	def set_current_position(self,Xp1):
		self.Xp1 = Xp1

	def temporary_set_var(self,var):
		self.cov_fact = var

	def set_obs(self,Y):
		self.Yp1 = Y

	def get_obs(self):
		return self.Yp1


class density_particle:
	def __init__(self,n_particle = 15, param_ref = particle_parameters(True),param_p = particle_parameters(False),isverbose = False):

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
		for i_p in range(self.n_particle):
			self.X_ips[:,i_p] = self.liste_p[i_p].get_current_position()

		self.w = np.zeros(n_particle)

	def __len__(self):
		return self.n_particle

	def next_step(self):
	#it of main particle
		self.p_ref.next_step()
		obs = self.p_ref.get_obs()
		self.current_estimate_position[:] = 0
# next step
		for i_p in range(0,self.n_particle):
#			print(i_p)
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
#		print(w)


#ressample	
		doressample = False
		if np.random.uniform()<0.005: # force ressample from time to time
			doressample = True
		if self.w.min() < 1e-8: # one particle is neglectible
			doressample = True
		if self.w.max() > 0.3: # one particle is too important
			doressample = True

		if doressample is True:
			Xrs,wrs,perm = ressample(self.X_ips,self.w)
			if self.isverbose:
				s = "resampling: " + str(np.unique(perm))
				print(s)
#			s = "resampling: " + str(np.unique(perm)) + ' weights: ' + str(self.w)
#			print(s)
			for i_p in range(0,self.n_particle):
				self.liste_p[i_p].set_weight(wrs[i_p])
				self.liste_p[i_p].set_current_position(Xrs[:,i_p])

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


