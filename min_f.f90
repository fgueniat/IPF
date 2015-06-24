FUNCTION F_proba(XX,t)
C Probability that X is a solution of F_dyn and coherent with measure y
C XX: contains Xn+1 and X05
C n: dim of X
	COMMON n
	COMMON g
	COMMON s
	COMMON dt
	COMMON Xn
	COMMON yp1
	REAL :: P_int, P_obs
C F_dyn such as dxdt = F_dyn(x,t). 
c Something as Xnp1 = F_dyn(Xn) is also possible - then dt has to be removed in the numerator.

c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Probability than Xp1,X05 is a solution of the dynamical system at hand
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	fx = F_dyn(Xn,t)
	X05 = XX(0:n)
	Xp1 = XX(n+1:2*n)
c traduces the 1st step of Klauder perterson scheme
	P_int = NORM2( X05 - Xn - dt*fx )**2.0 / (2.0 * dt * g * g)
c traduces the 2nd step of Klauder perterson scheme
	P_int = P_int + NORM2(Xp1 - Xn - dt*( fx + F_dyn(X05) )/20 )**2.0 / (2.0 * dt * g * g)


c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Probability than Xp1 compatible with the measurement y
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	P_obs = NORM2( H_obs(Xp1,t) - yp1 )**2.0 / (2.0 * dt * s * s)

	F_proba = P_int + P_obs
	
END FUNCTION F_proba



FUNCTION F_1d(XX,t)
c find a sample using argmin F_proba and min(F_proba)
	COMMON Xn
	COMMON minF
	COMMON Xmin
	COMMON Eps

	F_1d = abs(F_proba(XX,t) - minF - 0.5* NORM2(Eps))

END FUNCTION F_1d


SUBROUTINE Ressample(X,w,Xrs,wrs):
C from M.S. Arulampalam et al., A tutorial on particle filters for online nonlinear / non-Gaussian Bayesian tracking, IEEE Trans. Sign. Proc.
	COMMON n_particle
	COMMON n

	REAL, DIMENSION(n_particle) :: u
	REAL, DIMENSION(n_particle) :: cdf
	REAL, DIMENSION(n_particle) :: perm
	INTEGER :: ir
	INTEGER :: ex

c	REAL, DIMENSION(n_particle) :: wrs
c	REAL, DIMENSION(n, n_particle) :: Xrs

	wrs(:) = 1.0/n_particle
	cdf(:) = 0.0d0
	u0 = random_number()/n_particle
	u(:) = 0.0d0
	perm(:) = 0.0d0

c	Xrs = X
	DO i = 1,n_particle:
		cdf(i) = cdf(i-1) + w(i)
	END DO
	
	ir=0
	DO j =0,n_particle
		u(j) = u0 + 1.0*j/n_particle
		ex = 0

		DO WHILE ( u(j) > cdf(ir) .OR. ex == 1 ):
			IF ( ir<n_particle-1) THEN
				ir=ir+1
			ELSE
				ex =1
		END DO

		Xrs(:,j) = X(:,ir)
		perm(j)=ir
	END DO

END SUBROUTINE Ressample




REAL FUNCTION F_dyn(x,t)
c	COMMON n
c	REAL, DIMENSION(n) :: F_dyn
	REAL, DIMENSION(3) :: F_dyn
	REAL, PARAMETER :: beta = 8.0/3.0
	REAL, PARAMETER :: rho = 28.0
	REAL, PARAMETER :: sigma = 10.0
	
	F_dyn(0) = sigma*(x(1) - x(0))
	F_dyn(1) = X(0)*(rho-x(2)) - x(1)
	F_dyn(2) = x(0)*x(1) - beta*x(2)

END FUNCTION F_dyn

REAL FUNCTION H_obs(x,t)
c	COMMON n
c	REAL, DIMENSION(n) :: H_obs
	REAL, DIMENSION(3) :: H_obs
	REAL, DIMENSION(3) :: y_scalar

	call random_number(y_scalar)
	H_obs = x + y_scalar

END FUNCTION H_obs

