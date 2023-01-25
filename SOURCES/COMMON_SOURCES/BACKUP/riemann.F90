! Authors: Jean-Luc Guermond and Bojan Popov, Texas A&M, Jan 22, 2016
MODULE lambda_module
  PUBLIC :: lambda
  PRIVATE
  INTEGER,      PARAMETER :: Mgas=3, Mgasptwo = Mgas+2
  REAL(KIND=8), PARAMETER :: gamma=(Mgas+2.d0)/Mgas
  REAL(KIND=8), PARAMETER :: exponent=(gamma-1)/(2*gamma), one_m_exponent = 1-exponent
CONTAINS

  SUBROUTINE lambda(rhol,ul,pl,rhor,ur,pr,lambdal,lambdar)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: rhol, ul, pl, rhor, ur, pr
    REAL(KIND=8), INTENT(OUT):: lambdal, lambdar 
    REAL(KIND=8)             :: ratio, pstar
    REAL(KIND=8)             :: pmin, pmax
    REAL(KIND=8)             :: capAmin, capBmin, amin, amax, phimin, phimax
    REAL(KIND=8)             :: al, ar

    al=SQRT(gamma*pr/rhor)
    ar=SQRT(gamma*pl/rhol)
    pstar = al+ar-0.5d0*(gamma-1.0d0)*(ur-ul)
    if (pstar>0.0d0) THEN
       pmin=pr; pmax=pl; amin=ar; amax=al
       if (pl < pr) THEN 
          pmin=pl; pmax=pr; amin=al; amax=ar
       END if
       pstar = pstar/(amin + amax*(pmin/pmax)**exponent)
       pstar = pstar**(Mgasptwo)*pmin
    END if

    lambdal = ul - al*SQRT(1.d0+MAX((pstar-pl)/pl,0.d0)*one_m_exponent)
    lambdar = ur + ar*SQRT(1.d0+MAX((pstar-pr)/pr,0.d0)*one_m_exponent)

!!$    IF (pl .le. pr) THEN
!!$       capAmin = 2/((gamma+1)*rhol)
!!$       capBmin = pl*(gamma-1)/(gamma+1)
!!$       amin = SQRT(gamma*pl/rhol)
!!$       amax = SQRT(gamma*pr/rhor)
!!$       pmin = pl 
!!$       pmax = pr
!!$       al = amin
!!$       ar = amax
!!$    ELSE
!!$       capAmin = 2/((gamma+1)*rhor)
!!$       capBmin = pr*(gamma-1)/(gamma+1)
!!$       amin = SQRT(gamma*pr/rhor)
!!$       amax = SQRT(gamma*pl/rhol)
!!$       pmin = pr 
!!$       pmax = pl
!!$       ar = amin
!!$       al = amax
!!$    END IF
!!$    ratio= (pmin/pmax)**(exponent)
!!$    phimin = (2/(gamma-1.d0))*amax*(ratio-1.d0) + ur-ul
!!$    IF (phimin.GE.0) THEN
!!$       pstar = 0.d0
!!$       lambdal = ul-al
!!$       lambdar = ur+ar
!!$       RETURN
!!$    END IF
!!$    phimax = (pmax-pmin)*SQRT(capAmin/(pmax+capBmin)) + ur-ul
!!$    pstar = pmin*((amin+amax-(ur-ul)*(gamma-1)/2)/(amin + amax*ratio))**(Mgasptwo)
!!$    IF (phimax > 0.d0) THEN
!!$       pstar = MIN(pmax,pstar)
!!$    END IF
!!$    lambdal = ul - al*SQRT(1.d0+MAX((pstar-pl)/pl,0.d0)*one_m_exponent)
!!$    lambdar = ur + ar*SQRT(1.d0+MAX((pstar-pr)/pr,0.d0)*one_m_exponent)
!!$    RETURN
  END SUBROUTINE lambda

END MODULE lambda_module

