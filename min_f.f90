C Probability that X is a solution of F_dyn and coherent with measure y
SUBROUTINE F_proba(Xp1,X05,X,tp1,yp1,P)
	COMMON n
	COMMON g
	COMMON dt
	REAL :: P_int, P_obs
C F_dyn such as dxdt = F_dyn(x,t). 
c Something as Xnp1 = F_dyn(Xn) is also possible - then dt has to be removed in the numerator.

c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Probability than Xp1,X05 is a solution of the dynamical system at hand
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	fx = F_dyn(Xn)

c traduces the 1st step of Klauder perterson scheme
	P_int = NORM2( X05 -X - dt*fx )**2.0 / (2.0 * dt * g * g)
c traduces the 2nd step of Klauder perterson scheme
	P_int = P_int + NORM2(Xp1 - X - dt*( fx + F_dyn(X05) )/20 )**2.0 / (2.0 * dt * g * g)


c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Probability than Xp1 compatible with the measurement y
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	P_obs = NORM2( H_obs(Xp1,t) - y )**2.0 / (2.0 * dt * s * s)

	F_proba = P_int + P_obs
	
END SUBROUTINE F_proba


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

