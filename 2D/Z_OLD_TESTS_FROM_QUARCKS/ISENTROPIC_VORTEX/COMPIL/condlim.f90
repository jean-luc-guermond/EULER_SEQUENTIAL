MODULE boundary_conditions
  USE input_data
  REAL(KIND=8), PARAMETER :: gamma=1.4d0
  REAL(KIND=8), PARAMETER :: r0=0.15d0, x0=0d0, y0=0.0d0
  REAL(KIND=8), PARAMETER :: u_infty=1.d0, rho_infty=1.d0
  REAL(KIND=8) :: Mach, sound_spd, ubar_over_uinfty
  REAL(KIND=8) :: p_infty 
  REAL(KIND=8) :: ubar, chi
CONTAINS

  FUNCTION flux(comp,un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    INTEGER,                            INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(2,SIZE(un,2))          :: vv
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: H, u
    SELECT CASE(comp)
       CASE(1) 
          vv(1,:) = un(2,:)
          vv(2,:) = un(3,:)
       CASE(2)
          u = un(2,:)/un(1,:)
          vv(1,:) = un(2,:)*u + pressure(un)
          vv(2,:) = un(3,:)*u
       CASE(3)
          u = un(3,:)/un(1,:)
          vv(1,:) = un(2,:)*u 
          vv(2,:) = un(3,:)*u + pressure(un)
       CASE(4) 
          H = pressure(un) + un(4,:)
          vv(1,:) = (un(2,:)/un(1,:))*H
          vv(2,:) = (un(3,:)/un(1,:))*H
    CASE DEFAULT
       WRITE(*,*) ' BUG in flux'
       STOP
    END SELECT
  END FUNCTION flux

  FUNCTION rho_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv, radius_sq, z
    REAL(KIND=8) :: rtsq, phi, rsq
    INTEGER :: n, i

    DO n = 1, SIZE(rr,2)
       rsq = (rr(1,n)-x0-u_infty*inputs%time)**2 + (rr(2,n)-y0)**2
       z(n) = exp(1-rsq/(r0**2))
    END DO

    IF (inputs%type_test=='vrtx') THEN
       vv = rho_infty*(1-chi*z)**(1.d0/(gamma-1.d0))
    END IF
  END FUNCTION rho_anal

  FUNCTION press_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    INTEGER :: n
    vv = p_infty*(rho_anal(rr)/rho_infty)**gamma
  END FUNCTION press_anal

  FUNCTION vit_anal(comp,rr) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv, radius_sq, z
    REAL(KIND=8) :: rtsq, sqrt_mphip_over_r, rsq
    INTEGER :: n, i

    DO n = 1, SIZE(rr,2)
       rsq = (rr(1,n)-x0-u_infty*inputs%time)**2 + (rr(2,n)-y0)**2
       z(n) = exp(0.5d0*(1-rsq/(r0**2)))
    END DO
    
    IF (comp==1) THEN
       IF (inputs%type_test=='vrtx') THEN
          vv = u_infty - ubar*z*(rr(2,:)-y0)/r0
       END IF
    ELSE
       IF (inputs%type_test=='vrtx') THEN
          vv = ubar*z*(rr(1,:)-x0-u_infty*inputs%time)/r0
       END IF
    END IF
  END FUNCTION vit_anal

  FUNCTION E_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    vv = press_anal(rr)/(gamma-1.d0) &
         + rho_anal(rr)*(vit_anal(1,rr)**2+vit_anal(2,rr)**2)/2 
  END FUNCTION E_anal

  FUNCTION mt_anal(comp,rr) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    vv = rho_anal(rr)*vit_anal(comp,rr)
  END FUNCTION mt_anal
  
  FUNCTION pressure(un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: vv
    vv=(un(4,:)-0.5d0*(un(2,:)**2+un(3,:)**2)/un(1,:))*(gamma-1)
  END FUNCTION pressure
  
  FUNCTION sound_speed(un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: vv
    vv= SQRT(gamma*pressure(un)/un(1,:))
  END FUNCTION sound_speed

  SUBROUTINE init(un,rr)
    USE character_strings
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(4,SIZE(rr,2)),   INTENT(OUT):: un
    INTEGER, PARAMETER :: in_unit=21
    LOGICAL :: okay
  
    OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(in_unit, "===Mach number & ubar_over_uinfty===") 
    READ(in_unit,*) Mach, ubar_over_uinfty
    CLOSE(in_unit)
    
    sound_spd = u_infty/Mach 
    p_infty = (rho_infty/gamma)*sound_spd**2
    ubar = ubar_over_uinfty*u_infty
    chi=((gamma-1)/(2*gamma))*(rho_infty/p_infty)*ubar**2
    
    IF (inputs%if_restart) THEN
       OPEN(unit = 10, file = 'restart.'//inputs%file_name, form = 'unformatted', status = 'unknown')
       READ(10) inputs%time, un
       CLOSE(10)
       WRITE(*,*) ' inputs%time at restart', inputs%time 
    ELSE
       un(1,:) = rho_anal(rr)
       un(2,:) = mt_anal(1,rr)
       un(3,:) = mt_anal(2,rr)
       un(4,:) = E_anal(rr)
    END IF
  END SUBROUTINE init

  FUNCTION sol_anal(comp,rr,time,opt_out) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                                 INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    LOGICAL,      OPTIONAL,                  INTENT(IN) :: opt_out
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                 :: vv
    REAL(KIND=8) :: time
    IF (PRESENT(opt_out)) THEN
       SELECT CASE(comp)
       CASE(1)
          vv = rho_infty
       CASE(2)
          vv = rho_infty*u_infty
       CASE(3)
          vv = 0.d0
       CASE(4)
          vv = P_infty/(gamma-1) + 0.5d0*rho_infty*u_infty**2
       CASE DEFAULT
          WRITE(*,*) ' BUG in sol_anal'
          STOP
       END SELECT
       RETURN
    END IF
    SELECT CASE(comp)
    CASE(1)
       vv = rho_anal(rr)
    CASE(2)
       vv = mt_anal(1,rr)
    CASE(3)
       vv = mt_anal(2,rr)
    CASE(4)
       vv = E_anal(rr)
    CASE DEFAULT
       WRITE(*,*) ' BUG in sol_anal'
       STOP
    END SELECT
  END FUNCTION sol_anal


END MODULE boundary_conditions
