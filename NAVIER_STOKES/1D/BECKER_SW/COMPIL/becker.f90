MODULE newton
  PUBLIC :: newton_iter
  PRIVATE
CONTAINS
  SUBROUTINE newton_iter(v,vmin,vmax,psi_func,dpsi_func,tol)
    IMPLICIT NONE
    INTERFACE
       FUNCTION psi_func(vin) RESULT(vout)
         IMPLICIT NONE
         REAL(KIND=8), INTENT(IN) :: vin
         REAL(KIND=8) :: vout
       END FUNCTION psi_func
       FUNCTION dpsi_func(vin) RESULT(vout)
         IMPLICIT NONE
         REAL(KIND=8), INTENT(IN) :: vin
         REAL(KIND=8) :: vout
       END FUNCTION dpsi_func
    END INTERFACE
    REAL(KIND=8), INTENT(IN) :: vmin, vmax, tol
    REAL(KIND=8) :: v
    REAL(KIND=8) :: norm, psi, vnext, dp
    INTEGER :: iter
    norm = MAX(ABS(psi_func(0.25*vmin+0.75*vmax)), ABS(psi_func(0.75*vmin+0.25*vmax)))
    !===Test trivial cases
    IF (psi_func(vmin+tol).LE.tol) THEN
       v = vmin+tol
       RETURN
    ELSE IF (psi_func(vmax-tol).GE.tol) THEN
       v=vmax-tol
       RETURN
    END IF
    
    iter = 0
    psi =  psi_func(v)
    DO WHILE (ABS(psi).GE.tol*norm)
       !dp = (psi_func(v+sqrt(tol)) - psi_func(v-sqrt(tol)))/(2*sqrt(tol))
       dp = dpsi_func(v)
       vnext = v - psi/dp
       IF (ABS(vnext-v).LT.tol*(vmax+vmin)/2) THEN
          RETURN
       ELSE IF (vnext<vmin) THEN
          v = (vmin+v)/2
       ELSE IF (vnext>vmax) THEN
          v = (vmax+v)/2
       ELSE
          v = vnext
       END IF
       psi = psi_func(v)
       iter = iter+1
    END DO
  END SUBROUTINE newton_iter
END MODULE newton

MODULE becker
  USE NEWTON
  REAL(KIND=8), PUBLIC :: becker_vmin, becker_vmax, becker_const, &
       becker_expomax, becker_expomin, becker_x, becker_vorigin
  PUBLIC :: set_becker_parameters, psi, dpsi
  PRIVATE
CONTAINS
  SUBROUTINE set_becker_parameters
    becker_expomax = becker_vmax/(becker_vmax-becker_vmin)
    becker_expomin = -becker_vmin/(becker_vmax-becker_vmin)
  END SUBROUTINE set_becker_parameters
  FUNCTION psi(v) RESULT(vout)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: v
    REAL(KIND=8) :: vout, ratio
    !ratio = ((becker_vmax-v)/(becker_vmax-becker_vorigin))**becker_expomax &
    !     * ((v-becker_vmin)/(becker_vorigin-becker_vmin))**becker_expomin
    !vout = becker_const*LOG(ratio)
    vout = becker_const*(becker_expomax*LOG(becker_vmax-v)              + becker_expomin*LOG(v-becker_vmin) &
                       - becker_expomax*LOG(becker_vmax-becker_vorigin) - becker_expomin*LOG(becker_vorigin-becker_vmin)) - becker_x
  END FUNCTION  psi
  FUNCTION dpsi(v) RESULT(vout)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: v
    REAL(KIND=8) :: vout
    vout = becker_const*(-becker_expomax/(becker_vmax-v) + becker_expomin/(v-becker_vmin))
  END FUNCTION  dpsi
END MODULE becker
