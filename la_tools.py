import numpy as np

def mdot(ABC):
	A = ABC[0]
	for i in range(1,len(ABC)):
		A = np.dot(A,ABC[i])
	return A

def gen_cholesky(H, eps_machine=1.0e-15, print_prefix=0, print_flag=0): 
#The Gill, Murray, and Wright modified Cholesky algorithm, from Practical Optimization, Academic Press, London, 1981
 
	ndim = len(H) 
	I = np.eye(ndim) 
	# Calculate gamma(A) and xi(A). 
	gamma = 0.0 
	xi = 0.0 
	for i in range(ndim): 
		gamma = np.max((np.abs(H[i, i]), gamma)) 
		for j in range(i+1, ndim): 
			xi = max(np.abs(H[i, j]), xi) 

# Identify delta and beta. 
	try:
		delta = eps_machine * np.max((gamma + xi, 1.0)) 
	except:
		print(gamma)
		print(xi)
		print(eps_machine)
		delta = eps_machine * np.max((gamma + xi, 1.0)) 

	if ndim == 1: 
		beta = np.sqrt( np.max((gamma, eps_machine)) )
	else: 
		beta = np.sqrt(  np.max((gamma, xi / np.sqrt(ndim**2 - 1.0), eps_machine))  ) 
 
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
		a = np.dot(p, np.dot(a, p)) 
		r = np.dot(r, p) 
 
# Calculate dj. 
		theta_j = 0.0 
		if j < ndim-1: 
			for i in range(j+1, ndim): 
				theta_j = np.max(  (theta_j, np.abs(a[j, i]))  ) 
		dj = np.max(  (np.abs(a[j, j]), (theta_j/beta)**2, delta)  ) 

# Calculate row j of r and update a. 
		r[j, j] = np.sqrt(dj)     # Damned sqrt introduces roundoff error. 
		for i in range(j+1, ndim): 
			r[j, i] = a[j, i] / r[j, j] 
			for k in range(j+1, i+1): 
				a[i, k] = a[k, i] = a[k, i] - r[j, i] * r[j, k]     # Keep matrix a symmetric. 
 
# Finally, the Cholesky decomposition of H. 
	L = np.dot(P, np.transpose(r))
	return L
