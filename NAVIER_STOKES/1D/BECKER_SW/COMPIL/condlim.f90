MODULE boundary_conditions
  USE input_data
  USE space_dim
  USE becker
  PUBLIC :: flux, sol_anal, pressure, init, rho_anal, press_anal, mt_anal, E_anal, vit_anal
  REAL(KIND=8), PUBLIC :: gamma, uninfty
  PRIVATE
  REAL(KIND=8) :: Mach, rho0, v0, m0, cv, cp, Pr, tol, Rinfty
CONTAINS

  FUNCTION flux(comp,un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    INTEGER,                            INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(k_dim,SIZE(un,2))      :: vv
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: H, u
    SELECT CASE(comp)
    CASE(1) 
       vv(1,:) = un(2,:)
    CASE(2)
       u = un(2,:)/un(1,:)
       vv(1,:) = un(2,:)*u + pressure(un)
    CASE(3) 
       H = pressure(un) + un(3,:)
       vv(1,:) = (un(2,:)/un(1,:))*H
    CASE DEFAULT
       WRITE(*,*) ' BUG in flux'
       STOP
    END SELECT
  END FUNCTION flux
  
  SUBROUTINE init(un,rr,time)
    USE input_data_ns
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),                INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(k_dim+2,SIZE(rr,2)), INTENT(OUT):: un
    REAL(KIND=8) :: time
    !===Ref J. Fluid Mech. (2013), vol. 726, R4
    uninfty = 0.2d0
    tol = 1.d-10
    Pr = 0.75d0
    gamma = 1.4
    cv = 1.d0/(gamma-1.d0)
    cp =  gamma*cv
    Mach = 3.d0
    rho0 = 1.0
    v0 = 1.d0
    m0 = rho0*v0
    Rinfty = (gamma+1)/(gamma-1) !===(3.8)
    inputs_ns%thermal_diff = cp*inputs_ns%mu_visc/Pr
    becker_vmax=v0
    becker_vmin=becker_vmax*(gamma-1+2/Mach**2)/(gamma+1) !===(2.10)
    becker_vorigin = sqrt(becker_vmax*becker_vmin) !===Origin at adiabatic sonic point ((3.6) to (3.7))
    becker_const=(2.d0/(gamma+1))*inputs_ns%thermal_diff/(m0*cv) !===(3.4) and (3.6)
    CALL set_becker_parameters
    write(*,*) becker_vmax, becker_vmin 
    un(1,:) = rho_anal(rr,time)
    un(2,:) = mt_anal(1,rr,time)
    un(3,:) = E_anal(rr,time)
  END SUBROUTINE init
  
  FUNCTION rho_anal(rr,time) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv
    REAL(KIND=8) :: time
    REAL(KIND=8) :: xi, cl, cr, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    vv = m0/(vit_anal(1,rr,time)-uninfty)
  END FUNCTION rho_anal

  FUNCTION press_anal(rr,time) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv, zz
    REAL(KIND=8) :: time
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    zz = vit_anal(1,rr,time)-uninfty
    vv = (m0/zz)*(Rinfty*becker_vmin*becker_vmax - zz**2)/(2*cp) !===(3.7)
  END FUNCTION press_anal

  FUNCTION vit_anal(comp,rr,time) RESULT(vv)
    USE newton
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    REAL(KIND=8) :: time
    INTEGER :: n
    REAL(KIND=8) :: vold
    IF (SIZE(vv)==0) RETURN
    vv(1) = 0.99*becker_vmax+0.01*becker_vmin
    DO n = 1, SIZE(rr,2)
       IF (n>1) THEN
          vv(n) = vold
       END IF
       becker_x = rr(1,n)-uninfty*time
       CALL newton_iter(vv(n),becker_vmin,becker_vmax,psi,dpsi,tol)
       vold = vv(n)
       vv(n) = vv(n)+uninfty
    END DO
  END FUNCTION vit_anal

  FUNCTION E_anal(rr,time) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    REAL(KIND=8) :: time
    vv = press_anal(rr,time)/(gamma-1.d0) &
         + rho_anal(rr,time)*(vit_anal(1,rr,time)**2)/2
  END FUNCTION E_anal

  FUNCTION mt_anal(comp,rr,time) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    REAL(KIND=8) :: time
    vv = rho_anal(rr,time)*vit_anal(comp,rr,time)
  END FUNCTION mt_anal

  FUNCTION pressure(un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: vv
       vv=(un(3,:)-0.5d0*(un(2,:)**2)/un(1,:))*(gamma-1)
  END FUNCTION pressure
  
  FUNCTION sol_anal(comp,rr,time)RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                            INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))            :: vv
    REAL(KIND=8) :: time
    SELECT CASE(comp)
       CASE(1)
          vv = rho_anal(rr,time)
       CASE(2)
          vv = mt_anal(1,rr,time)
       CASE(3)
          vv = E_anal(rr,time)
           CASE DEFAULT
       WRITE(*,*) ' BUG in sol_anal'
       STOP
    END SELECT
  END FUNCTION sol_anal
END MODULE boundary_conditions
