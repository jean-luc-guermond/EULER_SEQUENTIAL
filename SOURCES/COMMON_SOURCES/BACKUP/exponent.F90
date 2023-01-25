MODULE exponent
PUBLIC :: power
REAL(KIND=8), PARAMETER :: tol=1.d-2
CONTAINS

  FUNCTION power(x,a) RESULT(y) !Computes x^(1/a), 0<x<1
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: x
    INTEGER,      INTENT(IN) :: a
    REAL(KIND=8)             :: y
    REAL(KIND=8)             :: p1, p2, phi1, phi2, phi11, phi22, phi12, phi112, phi221
    IF (x==0.d0) THEN
       y = 0.d0
       RETURN
    END IF
    p1 = SQRT(x)
    p1= sqrt(p1)
    p2= sqrt(p1)
    !p2 = 1.d0
    DO WHILE(.TRUE.) 

       IF (p2-p1 .LE. tol) THEN
          y = p1
          RETURN
       END IF
       phi1 = p1**a
       phi2 = p2**a
       phi11 = a*phi1/p1
       phi22 = a*phi2/p2
       phi1 = phi1 - x
       phi2 = phi2 - x
       IF (phi2<0.d0) THEN
          y = p2
          RETURN
       END IF
       IF (phi1>0.d0) THEN
          y = p1
          RETURN
       END IF
       phi12 = (phi2-phi1)/(p2-p1) 
       phi112 = (phi12-phi11)/(p2-p1)
       phi221 = (phi22-phi12)/(p2-p1)
       p1 = p1 - 2*phi1/(phi11 + SQRT(phi11**2 - 4*phi1*phi112))
       p2 = p2 - 2*phi2/(phi22 + SQRT(phi22**2 - 4*phi2*phi221))
    END DO
  END FUNCTION power

END MODULE exponent
 
PROGRAM check
  USE exponent
  IMPLICIT NONE
  INTEGER, PARAMETER :: a=7, nmax=1000000
  INTEGER      :: n
  REAL(KIND=8) :: x, y, z=-1.d0, t1, t2
  CALL CPU_TIME(t1)
  DO n = 1, nmax
     x = RAND(0)
     y = power(x,a)
  END DO
  CALL CPU_TIME(t2)
  WRITE(*,*) ' CPU ', t2-t1
  DO n = 1, nmax
     x = RAND(0)
     y = ABS(power(x,a) - x**(1.d0/a))
     z = MAX(y,z)
  END DO
  WRITE(*,*) 'ERROR = ', z, power(x,a), x**(1.d0/a)
  CALL CPU_TIME(t1)
  DO n = 1, nmax
     x = RAND(0)
     y = x**(1.d0/a)
  END DO
  CALL CPU_TIME(t2)
  WRITE(*,*) ' CPU ', t2-t1
END PROGRAM check
