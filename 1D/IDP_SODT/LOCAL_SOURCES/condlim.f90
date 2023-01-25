MODULE boundary_conditions
  USE space_dim
  PUBLIC :: flux, sol_anal, pressure, init, rho_anal, press_anal, mt_anal, E_anal
  REAL(KIND=8), PUBLIC :: gamma
  PRIVATE
  REAL(KIND=8) :: rhol, pl, ul
  REAL(KIND=8) :: rhor, pr, ur
  REAL(KIND=8) :: long, x0, x1
CONTAINS

  FUNCTION flux(comp,un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    INTEGER,                            INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(SIZE(un,1),k_dim)      :: vv
    REAL(KIND=8), DIMENSION(SIZE(un,1))            :: H, u
    SELECT CASE(comp)
    CASE(1) 
       vv(:,1) = un(:,2)
    CASE(2)
       u = un(:,2)/un(:,1)
       vv(:,1) = un(:,2)*u + pressure(un)
    CASE(3) 
       H = pressure(un) + un(:,3)
       vv(:,1) = (un(:,2)/un(:,1))*H
    CASE DEFAULT
       WRITE(*,*) ' BUG in flux'
       STOP
    END SELECT
  END FUNCTION flux

  SUBROUTINE init(un,rr)
    USE def_of_gamma
    USE lambda_module
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),                INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2),k_dim+2), INTENT(OUT):: un
    REAL(KIND=8) :: cl, cr, rho_plus
    REAL(KIND=8) :: in_state(2), in_data(3), out_state(4)

    gamma = 1.4d0                                                                
    long=1.d0                                                                    
    x0 = long*0.5d0                                                              
    rhol=1.d0                                                                    
    rhor=0.125d0                                                                 
    pl = 1.d0                                                                    
    pr = 0.1d0                                                                   
    ul = 0.d0                                                                    
    ur = 0.d0
    
    un(:,1) = rho_anal(rr)
    un(:,2) = mt_anal(1,rr)
    un(:,3) = E_anal(rr)
    CALL set_gamma_for_riemann_solver(gamma)
  END SUBROUTINE init

  FUNCTION rho_anal(rr) RESULT(vv)
    USE input_data
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv
    REAL(KIND=8) :: xi, cl, cr, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    DO n = 1, SIZE(vv)
       IF (rr(1,n)<x0) THEN
          vv(n) = rhol  
       ELSE
          vv(n) = rhor
       END IF
    END DO
  END FUNCTION rho_anal

  FUNCTION press_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    REAL(KIND=8) :: xi, cl, cr, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    DO n = 1, SIZE(vv)
       IF (rr(1,n)<x0) THEN
          vv(n) = pl
       ELSE
          vv(n) = pr
       END IF
    END DO
  END FUNCTION press_anal

  FUNCTION vit_anal(comp,rr) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    REAL(KIND=8) :: xi, cl, cr, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    DO n = 1, SIZE(vv)
       IF (rr(1,n)<x0) THEN
          vv(n) = ul 
       ELSE
          vv(n) = ur
       END IF
    END DO
  END FUNCTION vit_anal

  FUNCTION E_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    vv = press_anal(rr)/(gamma-1.d0) &
         + rho_anal(rr)*(vit_anal(1,rr)**2)/2
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
    REAL(KIND=8), DIMENSION(SIZE(un,1))            :: vv
    vv=(un(:,3)-0.5d0*(un(:,2)**2)/un(:,1))*(gamma-1)
  END FUNCTION pressure

  FUNCTION sol_anal(comp,rr,time) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                                 INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                 :: vv
    REAL(KIND=8) :: time
    SELECT CASE(comp)
    CASE(1)
       vv = rho_anal(rr)
    CASE(2)
       vv = mt_anal(1,rr)
    CASE(3)
       vv = E_anal(rr)
    CASE DEFAULT
       WRITE(*,*) ' BUG in sol_anal'
       STOP
    END SELECT
  END FUNCTION sol_anal
END MODULE boundary_conditions
