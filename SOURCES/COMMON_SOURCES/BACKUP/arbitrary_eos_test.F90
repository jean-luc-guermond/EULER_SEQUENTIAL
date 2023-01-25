! Authors: Bennett Clayton, Jean-Luc Guermond and Bojan Popov, Texas A&M, April 5, 2021
MODULE arbitrary_eos_lambda_module
  use vdw !=hack
  IMPLICIT NONE
  PUBLIC :: lambda_arbitrary_eos
  REAL(KIND=8), PUBLIC :: b_covolume = .5d0
  PRIVATE :: phi
  PRIVATE
  INTEGER, PARAMETER:: NUMBER = KIND(1.d0)
  REAL(KIND=NUMBER), PARAMETER :: zero = 0
  REAL(KIND=NUMBER), PARAMETER :: one = 1
  REAL(KIND=NUMBER), PARAMETER :: half = 0.5d0
  REAL(KIND=NUMBER), PARAMETER :: five_third = 5.d0/3.d0
  REAL(KIND=NUMBER) :: rhol, ul, pl, el
  REAL(KIND=NUMBER) :: rhor, ur, pr, er
  REAL(KIND=NUMBER) :: gammal, al, alphal, capAl, capBl, capCl, expol
  REAL(KIND=NUMBER) :: gammar, ar, alphar, capAr, capBr, capCr, expor
  REAL(KIND=NUMBER) :: p_min, rho_min, gamma_min, alpha_min
  REAL(KIND=NUMBER) :: p_max, rho_max, gamma_max, alpha_max
  REAL(KIND=NUMBER) :: gamma_lm, expo_lm
  REAL(KIND=NUMBER) :: gamma_uM, expo_uM
  REAL(KIND=NUMBER) :: numerator, vacuum
  CHARACTER(LEN=1)  :: gamma_min_index, gamma_lm_index 

CONTAINS

  SUBROUTINE lambda_arbitrary_eos(in_rhol,in_ul,in_el,in_pl,in_rhor,in_ur,in_er,in_pr,in_tol,no_iter,&
       lambda_maxl_out,lambda_maxr_out,pstar,k)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: in_rhol, in_el, in_rhor, in_er, in_tol
    REAL(KIND=8), INTENT(IN), TARGET :: in_ul, in_pl, in_ur, in_pr
    LOGICAL,      INTENT(IN) :: no_iter
    REAL(KIND=8), INTENT(OUT):: lambda_maxl_out, lambda_maxr_out, pstar
    INTEGER,      INTENT(OUT):: k
    REAL(KIND=NUMBER)        :: lambda_maxl, lambda_maxr, tol
    REAL(KIND=NUMBER)        :: csl,  csr, p1, phi1, phi11, p2, phi2, phi22, phi12, phi112, phi221
    LOGICAL                  :: check
    !===Initialization
    rhol= in_rhol
    ul = in_ul
    pl = in_pl
    el = in_el
    rhor =in_rhor
    ur = in_ur
    pr = in_pr
    er = in_er
    k = 0
    CALL init(rhol,el,pl,gammal,al,alphal,capAl,capBl,capCl,expol)
    CALL init(rhor,er,pr,gammar,ar,alphar,capAr,capBr,capCr,expor)
!===
    !lambda_maxl = ul-sqrt(gammal*pl/(rhol*(1-b_covolume*rhol)))
    !lambda_maxr = ur+sqrt(gammar*pr/(rhor*(1-b_covolume*rhor)))
    !===CHECK POSITIVITY OF SOUND SPEED
    !csl =(gamma_vdw*(pl+avdw*rhol**2)/(rhol*(1-bvdw*rhol)) - 2*avdw*rhol)
    !csr =(gamma_vdw*(pr+avdw*rhor**2)/(rhor*(1-bvdw*rhor)) - 2*avdw*rhor)
    !IF (MIN(csl,csr)<0.d0) THEN
    !   WRITE(*,*) ' COMPLEX sound speed', p1, p2, ul, ur
    !   STOP
    !END IF
    !lambda_maxl_out = ul - SQRT(csl)
    !lambda_maxr_out = ur + SQRT(csr)
    !RETURN
!===
    IF (pl.LE.pr) THEN
       p_min     = pl
       rho_min   = rhol
       gamma_min = gammal
       gamma_min_index = 'l'
       alpha_min = alphal
       p_max     = pr
       rho_max   = rhor
       gamma_max = gammar
       alpha_max = alphar
    ELSE
       p_min     = pr
       rho_min   = rhor
       gamma_min = gammar
       gamma_min_index = 'r'
       alpha_min = alphar
       p_max     = pl
       rho_max   = rhol
       gamma_max = gammal
       alpha_max = alphal
    END IF
    IF (gammal.LE.gammar) THEN
       gamma_lm = gammal
       gamma_lm_index = 'l'
       gamma_uM = gammar 
    ELSE
       gamma_lm = gammar
       gamma_lm_index = 'r'
       gamma_uM = gammal
    END IF
    expo_lm = (gamma_lm-1)/(2*gamma_lm)  
    expo_uM = (gamma_uM-1)/(2*gamma_uM) 
    numerator = alphal+alphar-ur+ul
    vacuum = capCl+capCr+ul-ur

    !===Initialize p1 and p2
    CALL initialize_p1_p2(p1,p2)
 
    IF (no_iter) THEN
       CALL no_iter_update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p2,lambda_maxl,lambda_maxr)
       pstar = p2
       lambda_maxl_out = lambda_maxl
       lambda_maxr_out = lambda_maxr
    ELSE
       !===Iterations
       p1 = MAX(p1,p2-phi(p2)/phi_prime(p2))
       tol = in_tol
       DO WHILE(.TRUE.)
          CALL update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p1,p2,tol,&
               lambda_maxl,lambda_maxr,check)
          pstar = p2
          IF (check) RETURN
          phi1 =  phi(p1)
          phi11 = phi_prime(p1)
          phi2 =  phi(p2)
          phi22 = phi_prime(p2)
          IF (phi1>0.d0) THEN
             lambda_maxl = lambdaz(ul,pl,al,gammal,p1,-1)
             lambda_maxr = lambdaz(ur,pr,ar,gammar,p1, 1)
             pstar = p1
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
    END IF
  END SUBROUTINE lambda_arbitrary_eos

  SUBROUTINE init(rho,e,p,gamma,a,alpha,capA,capB,capC,expo)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: rho, e, p
    REAL(KIND=NUMBER), INTENT(OUT) :: gamma, a, alpha, capA, capB, capC, expo
    REAL(KIND=NUMBER) :: x
    x = 1-b_covolume*rho
    gamma = 1 + p*x/(rho*e) !=== e is the specific internal energy
    a = SQRT(gamma*p/(rho*x))
    capC = 2*a*x/(gamma-1)
    alpha = cc(gamma)*capC
    capA = 2*x/((gamma+1)*rho)
    capB = p*(gamma-1)/(gamma+1)
    expo = (gamma-1)/(2*gamma)
  CONTAINS
    function cc(gamma) result(vv)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN) :: gamma
      REAL(KIND=NUMBER)             :: vv
      IF (gamma.LE.1) THEN
         WRITE(*,*) "BUG: gamma .LE. 1", gamma, p
         STOP
      ELSE IF (gamma .LE. five_third) THEN
         vv = 1.d0 
      ELSE IF (gamma .LE. 3) THEN
         vv = SQRT((3*gamma+11)/(6*(gamma+1)))
      ELSE
         expo = (4-2*gamma)/(gamma-1)
         vv = SQRT(0.5d0+2*3**expo/(gamma-1))
      END IF
    END function cc
  END SUBROUTINE init

  SUBROUTINE no_iter_update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p2,lambda_maxl,lambda_maxr)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: ul, pl, al, gammal, ur, pr, ar, gammar
    REAL(KIND=NUMBER), INTENT(IN)  :: p2
    REAL(KIND=NUMBER), INTENT(OUT) :: lambda_maxl, lambda_maxr
    REAL(KIND=NUMBER) :: v11, v32, lambda_max
    v11 = lambdaz(ul,pl,al,gammal,p2,-1)
    v32 = lambdaz(ur,pr,ar,gammar,p2,1)
    lambda_maxl = MAX(-v11,0.d0)
    lambda_maxr = MAX(v32,0.d0)
    lambda_max = MAX(lambda_maxl,lambda_maxr)
  END SUBROUTINE no_iter_update_lambda
  
  SUBROUTINE update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p1,p2,tol,lambda_maxl,lambda_maxr,check)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: ul, pl, al, gammal, ur, pr, ar, gammar
    REAL(KIND=NUMBER), INTENT(IN)  :: p1, p2, tol
    REAL(KIND=NUMBER), INTENT(OUT) :: lambda_maxl, lambda_maxr
    LOGICAL,      INTENT(OUT) :: check
    REAL(KIND=NUMBER) :: v11, v12, v31, v32, lambda_max, err1, err3
    v11 = lambdaz(ul,pl,al,gammal,p2,-1)
    v12 = lambdaz(ul,pl,al,gammal,p1,-1)
    v31 = lambdaz(ur,pr,ar,gammar,p1,1)
    v32 = lambdaz(ur,pr,ar,gammar,p2,1)
    lambda_maxl = MAX(-v11,0.d0)
    lambda_maxr = MAX(v32,0.d0)
    lambda_max = MAX(lambda_maxl,lambda_maxr)
    err3 =  abs(v32 - v31)/lambda_max
    err1 =  abs(v12 - v11)/lambda_max
    IF (MAX(err1,err3).LE.tol) THEN
       check = .TRUE.
    ELSE
       check = .FALSE.
    END IF
  END SUBROUTINE update_lambda

  FUNCTION lambdaz(uz,pz,az,gammaz,pstar,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: uz,pz,az,gammaz,pstar
    INTEGER,      INTENT(IN) :: z
    REAL(KIND=NUMBER)             :: vv
    vv = uz + z*az*SQRT(1+MAX((pstar-pz)/pz,0.d0)*(gammaz+1)/(2*gammaz))
  END FUNCTION lambdaz

  FUNCTION phi(p) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: p
    REAL(KIND=NUMBER)             :: vv
    vv = f(p,pl,capAl,capBl,capCl,expol) + f(p,pr,capAr,capBr,capCr,expor)+ur-ul
  END FUNCTION phi

  function f(p,pz,capAz,capBz,capCz,expoz) result(ff)
    REAL(KIND=NUMBER), INTENT(IN) :: p, pz, capAz, capBz, capCz, expoz
    REAL(KIND=NUMBER)             :: ff
    IF (p.LE.pz) THEN
       ff = capCz*((p/pz)**expoz-1)
    ELSE
       ff = (p-pz)*SQRT(capAz/(p+capBz))
    END IF
  END function f

  FUNCTION phi_prime(p) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: p
    REAL(KIND=NUMBER)             :: vv
    vv = fp(p,pl,capAl,capBl,capCl,expol) + fp(p,pr,capAr,capBr,capCr,expor)
  CONTAINS
    function fp(p,pz,capAz,capBz,capCz,expoz) result(ff)
      REAL(KIND=NUMBER), INTENT(IN) :: p, pz, capAz, capBz, capCz, expoz
      REAL(KIND=NUMBER)             :: ff
      IF (p.LE.pz) THEN
         ff = capCz*expoz*(p/pz)**(expoz-1)/pz
      ELSE
         ff = SQRT(capAz/(p+capBz))*(1-(p-pz)/(2*(capBz+p)))
      END IF
    END function fp
  END FUNCTION phi_prime

  SUBROUTINE initialize_p1_p2(p1,p2)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(OUT) :: p1, p2
    REAL(KIND=NUMBER) :: phat1, phat2, r, rp
    REAL(KIND=NUMBER) :: xl, xr, a, b, c
    IF (vacuum.LE.0.d0) THEN
       p1=0.d0
       p2=0.d0
    ELSE IF (0.d0.LE.phi(p_min)) THEN
       p1 = 0.d0 
       phat1 = p_min*(numerator/(alpha_min +alpha_max*(p_min/p_max)**expo_uM))**(1/expo_uM)
       p2 = min(p_min,phat1)
    ELSE IF (0.d0.LE.phi(p_max)) THEN
       p1 = p_min
       rp = p_min/p_max
       r = (rp)**((gamma_uM-gamma_lm)/(2*gamma_lm*gamma_uM))
       IF (gamma_min_index==gamma_lm_index) THEN
          phat1 = p_min*(numerator/(r*alpha_min +alpha_max*(rp)**expo_uM))**(1/expo_uM)
          phat2 = p_min*(numerator/(alpha_min +r*alpha_max*(rp)**expo_lm))**(1/expo_lm)
       ELSE
          phat1 = p_min*(numerator/(alpha_min +alpha_max*(rp)**expo_lm))**(1/expo_lm)
          phat2 = p_min*(numerator/(alpha_min +alpha_max*(rp)**expo_uM))**(1/expo_uM)
       END IF
       p2 = min(p_max,phat1,phat2)
    ELSE
       p1 = p_max
       p2 = p_min*(numerator/(alpha_min +alpha_max*(p_min/p_max)**expo_lm))**(1/expo_lm)
       xl = SQRT(capAl/(1+capBl/p_max))
       xr = SQRT(capAr/(1+capBr/p_max))
       a = xl+xr
       b = ur-ul
       c = -pl*xl-pr*xr
       phat2 = ((-b+sqrt(b*b-4*a*c))/(2*a))**2
       p2 = min(p2,phat2)
    END IF
  END SUBROUTINE initialize_p1_p2
  
END MODULE arbitrary_eos_lambda_module
