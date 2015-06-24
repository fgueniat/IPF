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

Function NORM2(X)
C already present in fortran 2008
	NORM2 = DOT_PRODUCT(X,X) ** 0.5
end FUNCTION NORM2


FUNCTION F_1D(l,t)
c Used for identifying sample using argmin F_proba and min(F_proba)
	COMMON Xmin
	COMMON Eps
	COMMON L

	F_1D = F_proba( Xmin + l * TRANSPOSE(L) * Eps / NORM2(Eps),t )

END FUNCTION F_1D

FUNCTION Weight(X,l,t)
c find the weight associated with a sample
	COMMON Xn
	COMMON minF
	COMMON Xmin
	COMMON Eps
	COMMON w
	COMMON L
	COMMON n
	REAL :: rho = NORM2(Eps)
	REAL :: dl = 0.01
	REAL :: dldp
	real :: det

	det = GETDET(L,n)

	Fp = F_1D( l+dl,t )
	Fm = F_1D( l,t )
	dldp = dl / ( 2.0 * (Fp - Fm) )

	J = 2.0 * ABS(det) * (rho ** (1-n) ) * ABS( (l ** (n-1)) * dldp )

	Weight = w*EXP(-minF)*J
END FUNCTION F_1D


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
	u0 = RANDOM_NUMBER()/n_particle
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

	call RANDOM_NUMBER(y_scalar)
	H_obs = x + y_scalar

END FUNCTION H_obs




REAL*8 FUNCTION GETDET(A, N)
C A function written in FORTRAN77 to calculate determinant of a square matrix

C Passed parameters:
C A = the matrix
C N = dimension of the square matrix
       
C A modification of a code originally written by Ashwith J. Rego, available from http://www.dreamincode.net/code/snippet1273.htm
C Modified by Syeilendra Pramuditya, available from http://wp.me/p61TQ-zb
C Last modified on January 13, 2011 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ELEM(N,N),A(N,N)
      REAL*8 M, TEMP
      INTEGER I, J, K, L
      LOGICAL DETEXISTS

      DO I=1,N
      DO J=1,N
      ELEM(I,J)=A(I,J)
      END DO
      END DO
        
      DETEXISTS = .TRUE.
      L = 1
      !CONVERT TO UPPER TRIANGULAR FORM
      DO K = 1, N-1
	  IF (DABS(ELEM(K,K)).LE.1.0D-20) THEN
	  DETEXISTS = .FALSE.
	  DO I = K+1, N
	  IF (ELEM(I,K).NE.0.0) THEN
      
	  DO J = 1, N
	  TEMP = ELEM(I,J)
	  ELEM(I,J)= ELEM(K,J)
	  ELEM(K,J) = TEMP
	  END DO
      
	  DETEXISTS = .TRUE.
	  L=-L
	  EXIT
      
	  END IF
      
	  END DO
	  IF (DETEXISTS .EQV. .FALSE.) THEN
	  GETDET = 0
	  RETURN
	  END IF
      END IF
      DO J = K+1, N
	  M = ELEM(J,K)/ELEM(K,K)
	  DO I = K+1, N
	  ELEM(J,I) = ELEM(J,I) - M*ELEM(K,I)
	  END DO
	  END DO
      END DO
	
      !CALCULATE DETERMINANT BY FINDING PRODUCT OF DIAGONAL ELEMENTS
      GETDET = L
      DO I = 1, N
	  GETDET = GETDET * ELEM(I,I)
      END DO
	
      END        
