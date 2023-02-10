PROGRAM t
  IMPLICIT NONE
  INTEGER,      PARAMETER :: N=300
  REAL(KIND=8), PARAMETER :: angle=6.d0/180.d0, theta=angle*ACOS(-1.D0), l_oat15a=2.30
  INTEGER      :: beg, i
  REAL(KIND=8) :: x(300), y(300), ratio
  
  OPEN(unit=10,file='data',form='formatted', status='unknown')
  OPEN(unit=11,file='data_o',form='formatted', status='unknown')
  WRITE(11,*) N+1
  DO i = 1, N
     READ(10,*) beg, x(i), y(i)
  END DO

  ratio = l_oat15a/(MAXVAL(x) - MINVAL(x))
  write(*,*) (MAXVAL(x) - MINVAL(x)), ratio
  x = x*ratio
  y = y*ratio
  DO i = N,1,-1
     WRITE(11,*) COS(theta)*x(i)+SIN(theta)*y(i), -SIN(theta)*x(i)+COS(theta)*y(i)
  END DO
  WRITE(11,*) COS(theta)*x(N)+SIN(theta)*y(N), -SIN(theta)*x(N)+COS(theta)*y(N)
END PROGRAM t
