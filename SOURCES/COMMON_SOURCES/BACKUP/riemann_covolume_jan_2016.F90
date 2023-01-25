! Authors: Jean-Luc Guermond and Bojan Popov, Texas A&M, Nov 2, 2015
MODULE lambda_covolume_module
  PUBLIC :: lambda
  PRIVATE
  INTEGER,      PARAMETER :: Mgas=5
  REAL(KIND=8), PARAMETER :: gamma=(Mgas+2.d0)/Mgas, b=0.0d0
  REAL(KIND=8), PARAMETER :: exponent=(gamma-1)/(2*gamma), one_m_exponent = 1-exponent
  REAL(KIND=8)            :: al, capAl, capBl, covl, ar, capAr, capBr, covr
CONTAINS

  SUBROUTINE init(rhol,pl,rhor,pr)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: rhol, pl, rhor, pr
    al = SQRT(gamma*pl/rhol)
    capAl = 2/((gamma+1)*rhol)
    capBl = pl*(gamma-1)/(gamma+1)
    covl = SQRT(1-b*rhol)
    ar = SQRT(gamma*pr/rhor)
    capAr = 2/((gamma+1)*rhor)
    capBr = pr*(gamma-1)/(gamma+1)
    covr = SQRT(1-b*rhor)
    !exponent = (gamma-1)/(2*gamma)
    !one_m_exponent = 1-exponent !((gamma+1)/(2*gamma))
  END SUBROUTINE init

  SUBROUTINE lambda(tol,rhol,ul,pl,rhor,ur,pr,lambda_max,pstar,k)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: tol, rhol, ul, pl, rhor, ur, pr
    REAL(KIND=8), INTENT(OUT):: lambda_max, pstar
    INTEGER,      INTENT(OUT):: k
    REAL(KIND=8)             :: lambda_min, phimax, ptilde, num, denom, ratio
    REAL(KIND=8)             :: phi1, phi11, phi12, phi112, phi2, phi22, phi221
    REAL(KIND=8)             :: p1, p2 , pmin, pmax, v11, v12, v31, v32
    REAL(KIND=8)             :: capAmin, capBmin, amin, covmin, acovmin, phimin
    REAL(KIND=8)             :: capAmax, capBmax, amax, covmax, acovmax
    REAL(KIND=8)             :: aovercovl, aovercovr
    !===Initialization
    CALL init(rhol,pl,rhor,pr)
!return
    IF (pl .le. pr) THEN
       pmin = pl 
       capAmin = 2/((gamma+1)*rhol)
       capBmin = pl*(gamma-1)/(gamma+1)
       !acovmin = SQRT(gamma*pl*(1-b*rhol)/rhol)
       acovmin = SQRT(gamma*pl/rhol)
       pmax = pr
       !acovmax = SQRT(gamma*pr*(1-b*rhor)/rhor)
       acovmax = SQRT(gamma*pr/rhor)
    ELSE
       pmin = pr 
       capAmin = 2/((gamma+1)*rhor)
       capBmin = pr*(gamma-1)/(gamma+1)
       !acovmin = SQRT(gamma*pr*(1-b*rhor)/rhor)
       acovmin = SQRT(gamma*pr/rhor)
       pmax = pl
       !acovmax = SQRT(gamma*pl*(1-b*rhol)/rhol)
       acovmax = SQRT(gamma*pl/rhol)
    END IF
    ratio = (pmin/pmax)**exponent
    phimin = (2/(gamma-1.d0))*acovmax*(ratio-1.d0) + ur-ul
    IF (phimin.GE.0) THEN
       pstar = 0.d0
       lambda_max = MAX(MAX(-(ul-al),0.d0), MAX(ur+ar,0.d0))
       RETURN
    END IF
    phimax = (pmax-pmin)*SQRT(capAmin/(pmax+capBmin)) + ur-ul
    pstar = pmin*((acovmin+acovmax - (ur-ul)*(gamma-1)/2)&
         /(acovmin + acovmax*ratio))**(Mgas+2)
    IF (phimax < 0.d0) THEN
       p1 = pmax
       p2 = pstar
    ELSE
       p1=pmin
       p2 = MIN(pmax,ptilde)
    END IF
    IF (pl.LE.pr) THEN
       aovercovl = acovmin
       aovercovr = acovmax
    ELSE
       aovercovr = acovmin
       aovercovl = acovmax
    END IF
    !aovercovl = SQRT(gamma*pl/((1-b*rhol)*rhol))
    !aovercovr = SQRT(gamma*pr/((1-b*rhor)*rhor))
    lambda_max = MAX(MAX(-lambdaz(ul,pl,aovercovl,p2,-1),0.d0), &
         MAX(lambdaz(ur,pr,aovercovr,p2,1),0.d0))

    num = al*covl+ar*covr+(ul-ur)*(gamma-1)/2
    denom = al*covl/pl**(exponent)+ar*covr/pr**(exponent)
    ptilde = (num/denom)**(Mgas+2)
    IF (ABS(pstar-ptilde)/ptilde .GE. 1.d-10) THEN
       WRITE(*,*) ' BUG ' 
       STOP
    END IF

    RETURN


    pmin = MIN(pl,pr)
    pmax = MAX(pl,pr)
    k = 0
    IF (phi(pmin,ul,pl,ur,pr).GE.0) THEN
       pstar = 0.d0
       lambda_max = MAX(MAX(-(ul-al),0.d0), MAX(ur+ar,0.d0))
       RETURN
    END IF
    phimax= phi(pmax,ul,pl,ur,pr)
    IF (phimax==0) THEN
       pstar = pmax
       lambda_max = MAX(MAX(-lambdaz(ul,pl,al/covl,pstar,-1),0.d0), &
            MAX(lambdaz(ur,pr,ar/covr,pstar,1),0.d0))
       RETURN
    END IF
    num = al*covl+ar*covr+(ul-ur)*(gamma-1)/2
    denom = al*covl/pl**(exponent)+ar*covr/pr**(exponent)
    !ptilde = (num/denom)**(1/exponent)
    ptilde = (num/denom)**(7)

    write(*,*) ABS(pstar-ptilde)/ptilde
    RETURN

    !ptilde = exp(log(num/denom)*7)
    IF (phimax < 0.d0) THEN
       p1 = pmax
       p2 = ptilde
    ELSE
       p1=pmin
       p2 = MIN(pmax,ptilde)
    END IF
    p1 = MAX(p1,p2-phi(p2,ul,pl,ur,pr)/phi_prime(p2,pl,pr))
    !TEST
    v11 = lambdaz(ul,pl,al/covl,p2,-1)
    v32 = lambdaz(ur,pr,ar/covr,p2,1)
    lambda_max = MAX(MAX(v32,0.d0),MAX(-v11,0.d0))
    RETURN
    !TEST
    !===Iterations
    DO WHILE(.TRUE.) 
       v11 = lambdaz(ul,pl,al/covl,p2,-1)
       v12 = lambdaz(ul,pl,al/covl,p1,-1)
       v31 = lambdaz(ur,pr,ar/covr,p1,1)
       v32 = lambdaz(ur,pr,ar/covr,p2,1)
       lambda_max = MAX(MAX(v32,0.d0),MAX(-v11,0.d0))
       lambda_min = MAX(MAX(MAX(v31,0.d0),MAX(-v12,0.d0)),0.d0)
       IF (lambda_min>0.d0) THEN
          IF (lambda_max/lambda_min -1.d0 .LE. tol) THEN
             pstar = p2
             RETURN
          END IF
       END IF
       phi1 =  phi(p1,ul,pl,ur,pr)
       phi11 = phi_prime(p1,pl,pr)
       phi2 =  phi(p2,ul,pl,ur,pr)
       phi22 = phi_prime(p2,pl,pr)
       IF (phi1>0.d0) THEN
          lambda_max = lambda_min
          RETURN
       END IF
       IF (phi2<0.d0) RETURN
       phi12 = (phi2-phi1)/(p2-p1) 
       phi112 = (phi12-phi11)/(p2-p1)
       phi221 = (phi22-phi12)/(p2-p1)
       p1 = p1 - 2*phi1/(phi11 + SQRT(phi11**2 - 4*phi1*phi112))
       p2 = p2 - 2*phi2/(phi22 + SQRT(phi22**2 - 4*phi2*phi221))
       k = k+1
    END DO
  END SUBROUTINE lambda

  FUNCTION lambdaz(uz,pz,az,pstar,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: uz,pz,az,pstar
    INTEGER,      INTENT(IN) :: z
    REAL(KIND=8)             :: vv
    vv = uz + z*az*SQRT(1+MAX((pstar-pz)/pz,0.d0)*one_m_exponent)
  END FUNCTION lambdaz

  FUNCTION phi(p,ul,pl,ur,pr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: p, ul, pl, ur, pr
    REAL(KIND=8)             :: vv, fl, fr
    IF (p>pl) THEN
       fl = (p-pl)*SQRT(capAl/(p+capBl))
    ELSE
       !fl = 2*al/(gamma-1)*((p/pl)**exponent-1)
       fl = 2*al/(gamma-1)*(exp(log(p/pl)*exponent)-1)
    END IF
    IF (p>pr) THEN
       fr = (p-pr)*SQRT(capAr/(p+capBr))
    ELSE
       !fr = 2*ar/(gamma-1)*((p/pr)**exponent-1)
       fr = 2*ar/(gamma-1)*(exp(log(p/pr)*exponent)-1)
    END IF
    vv = fl*covl + fr*covr + ur - ul
  END FUNCTION phi

  FUNCTION phi_prime(p,pl,pr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: p, pl, pr
    REAL(KIND=8)             :: vv, fl, fr
    IF (p>pl) THEN
       fl = SQRT(capAl/(p+capBl))*(1-(p-pl)/(2*(capBl+p)))
    ELSE
       !fl = (al/(gamma*pl))*(pl/p)**((gamma+1)/(2*gamma))
       !fl = (al/(gamma*pl))*(pl/p)**one_m_exponent
       fl = (al/(gamma*pl))*exp(log(pl/p)*one_m_exponent)
    END IF
    IF (p>pr) THEN
       fr = SQRT(capAr/(p+capBr))*(1-(p-pr)/(2*(capBr+p)))
    ELSE
       !fr = (ar/(gamma*pr))*(pr/p)**((gamma+1)/(2*gamma))
       !fr = (ar/(gamma*pr))*(pr/p)**one_m_exponent
       fr = (ar/(gamma*pr))*exp(log(pr/p)*one_m_exponent)
    END IF
    vv = fl*covl + fr*covr
  END FUNCTION phi_prime
END MODULE lambda_covolume_module

