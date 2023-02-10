MODULE boundary_conditions
  USE input_data
  USE space_dim
  REAL(KIND=8), PUBLIC :: gamma
  !===Sonic test
  !REAL(KIND=8), PRIVATE :: ratio=0.95d0
  !REAL(KIND=8), PRIVATE :: rhol=1.d0, pl=0.01d0, ul=-1.d0
  !REAL(KIND=8), PRIVATE :: rhor=1.d0, pr=0.01d0, ur=1.d0
  REAL(KIND=8), PRIVATE :: rhol=3.857d0, pl=10.333d0, ul=0.92d0
  REAL(KIND=8), PRIVATE :: rhor=1.d0, pr=1.d0, ur=3.55d0
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

  FUNCTION rho_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv
    INTEGER :: n
    REAL(KIND=8) :: x0, long
    SELECT CASE(inputs%type_test)
    CASE ('sodt') ! Sod
       long = 1.d0
       x0 = long*0.5d0
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = 1.d0
          ELSE
             vv(n) = 0.125d0
          END IF
       END DO
    CASE('cont') ! Contact
       long = 1.d0
       x0 = long*0.15
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = 2.d0 
          ELSE
             vv(n) = 1.d0
          END IF
       END DO
    CASE('lebl') ! Leblanc
       long = 1.d0
       x0 = long*0.5d0
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = 1.d0 
          ELSE
             vv(n) = 0.001d0
          END IF
       END DO
    CASE ('soni')
       long = 1.d0
       x0 = long*0.5d0
       !rhor=(ratio+(gamma-1)/(gamma+1))/(1+ ratio*(gamma-1)/(gamma+1))*rhol
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = rhol 
          ELSE
             vv(n) = rhor
          END IF
       END DO
    CASE DEFAULT
       WRITE(*,*) ' BUG in rho_anal'
       STOP
    END SELECT
  END FUNCTION rho_anal

  FUNCTION press_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    INTEGER :: n
    REAL(KIND=8) :: x0, long
    SELECT CASE(inputs%type_test)
    CASE ('sodt') ! Sod
       long = 1.d0
       x0 = long*0.5d0
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = 1.d0
          ELSE
             vv(n) = 0.1d0
          END IF
       END DO
    CASE('lebl') ! Leblanc
       long = 1.d0
       x0 = long*0.5d0
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = (gamma-1.d0)*1.d-1
          ELSE
             vv(n) = (gamma-1.d0)*1.d-7
          END IF
       END DO
       !vv = vv*(gamma-1)
    CASE('cont') ! Contact
       vv = 1.d0
    CASE ('soni')
       long = 1.d0
       x0 = long*0.5d0
       !pr = ratio*pl
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = pl 
          ELSE
             vv(n) = pr
          END IF
       END DO
    CASE DEFAULT
       WRITE(*,*) ' BUG in press_anal'
       STOP
    END SELECT
  END FUNCTION press_anal

  FUNCTION vit_anal(comp,rr) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    INTEGER :: n
    REAL(KIND=8) :: x0, long
    SELECT CASE(inputs%type_test)
    CASE ('sodt','lebl') ! Sod and Leblanc
       vv = 0.d0
    CASE('cont') ! Contact
       IF (comp==1) THEN
          vv = 1.d0
       END IF
    CASE ('soni')
       long = 1.d0
       x0 = long*0.5d0
       !ul=SQRT(0.5d0*(pl/rhol)*((gamma+1)*ratio+gamma-1))
       !rhor=(ratio+(gamma-1)/(gamma+1))/(1+ ratio*(gamma-1)/(gamma+1))*rhol
       !ur = rhol*ul/rhor
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = ul 
          ELSE
             vv(n) = ur
          END IF
       END DO
    CASE DEFAULT
       WRITE(*,*) ' BUG in vit_anal'
       STOP
    END SELECT
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
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: vv
    vv=(un(3,:)-0.5d0*(un(2,:)**2)/un(1,:))*(gamma-1)
  END FUNCTION pressure

  FUNCTION sound_speed(un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: vv
    vv= SQRT(gamma*pressure(un)/un(1,:))
  END FUNCTION sound_speed

  SUBROUTINE init(un,rr)
    USE def_of_gamma 
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(k_dim+2,SIZE(rr,2)),   INTENT(OUT):: un
    SELECT CASE(inputs%type_test)
    CASE ('sodt','cont') ! Sod
       gamma = 1.4d0
    CASE ('lebl')
       gamma = 5.d0/3.d0
    CASE ('soni')
       gamma = 1.4d0
    CASE DEFAULT
       WRITE(*,*) ' BUG in init'
       STOP
    END SELECT
    un(1,:) = rho_anal(rr)
    un(2,:) = mt_anal(1,rr)
    un(3,:) = E_anal(rr)
    CALL set_gamma_for_riemann_solver(gamma)
  END SUBROUTINE init

  FUNCTION sol_anal(comp,rr)RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                                 INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
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
