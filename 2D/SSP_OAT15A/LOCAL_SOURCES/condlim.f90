MODULE boundary_conditions
  USE input_data
  REAL(KIND=8) :: gamma=1.4d0, rho_infty=0.413509d0, p_infty=26.5d3, Mach=0.73
  REAL(KIND=8) :: u_infty
CONTAINS
  FUNCTION flux(comp,un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    INTEGER,                            INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(2,SIZE(un,2))            :: vv
    REAL(KIND=8), DIMENSION(SIZE(un,2))              :: H, u
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
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    vv = rho_infty
  END FUNCTION rho_anal

  FUNCTION press_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    INTEGER :: n
    vv = p_infty
  END FUNCTION press_anal

  FUNCTION vit_anal(k,rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    INTEGER,                             INTENT(IN) :: k
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    INTEGER :: n
    IF (k==1) THEN
       vv = u_infty
    ELSE IF (k==2) THEN
       vv = 0.d0
    ELSE 
       WRITE(*,*) ' BUG '
       STOP
    END IF

  END FUNCTION vit_anal

  FUNCTION E_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv, radius
    vv = press_anal(rr)/(gamma-1.d0) &
         + rho_anal(rr)*(vit_anal(1,rr)**2+vit_anal(2,rr)**2)/2 
  END FUNCTION E_anal

  FUNCTION mt_anal(k,rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    INTEGER,                             INTENT(IN) :: k
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    vv(:) = rho_anal(rr)*vit_anal(k,rr)
  END FUNCTION mt_anal

  FUNCTION pressure(un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: vv
    vv=(un(4,:)-0.5d0*(un(2,:)**2+un(3,:)**2)/un(1,:))*(gamma-1)
  END FUNCTION pressure

  SUBROUTINE init(un,rr)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(4,SIZE(rr,2)),   INTENT(OUT):: un

    u_infty = Mach*sqrt(gamma*p_infty/rho_infty)
    !===Initialization of un
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

  FUNCTION sol_anal(comp,rr,time,opt_out)RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: time
    INTEGER,                                 INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    LOGICAL,      OPTIONAL,                  INTENT(IN) :: opt_out
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
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
