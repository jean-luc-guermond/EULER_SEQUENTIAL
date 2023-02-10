PROGRAM test
  IMPLICIT NONE
  REAL(KIND=8) :: a, b, c, gamma, b_covolume, rhol, rhor, pl, pr, &
       ul, ur, rhoel, rhoer, gammal, gammar, lambdal, lambdar
  a = 1.d0
  b = 1.d0
  b_covolume=0.d0
  gamma=1.02d0
  
  rhol= 0.15d0
  rhor = 0.015d0
  pl = .5d0
  pr = .003d0
  ul = -4d0
  ur = -2d0
  rhoel = (pl+a*rhol**2)*(1-b*rhol)/(gamma-1)-a*rhol**2
  rhoer = (pr+a*rhor**2)*(1-b*rhor)/(gamma-1)-a*rhor**2

  write(*,*) 'rhol, ul, pl, el', rhol, ul, pl, rhoel/rhol
  write(*,*) 'rhor, ur, pr, er', rhor, ur, pr, rhoer/rhor
  
  gammal = 1 + pl*(1-b_covolume*rhol)/(rhoel)
  gammar= 1 + pr*(1-b_covolume*rhor)/(rhoer)

  lambdal = lbd(ul,gammal,pl,rhol,-1)
  lambdar = lbd(ul,gammal,pl,rhol,+1) 
  write(*,*) 'lambdal, lambdar', lambdal, lambdar
  
  lambdal = lbd(ul,gammal,pl,rhol,-1)
  lambdar = lbd(ur,gammar,pr,rhor,+1) 
  write(*,*) 'lambdal, lambdar', lambdal, lambdar

  lambdal = lbd(ur,gammar,pr,rhor,-1)
  lambdar = lbd(ur,gammar,pr,rhor,+1) 
  write(*,*) 'lambdal, lambdar', lambdal, lambdar
  
CONTAINS
  FUNCTION lbd(u,gamma,p,rho,i) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: u, gamma, p, rho, vv
    INTEGER :: i
    vv = u+i*sqrt(gamma*p/(rho*(1-b_covolume*rho)))
  END FUNCTION lbd
END PROGRAM test
