MODULE boundary_conditions
  USE input_data
  REAL(KIND=8) :: gamma=1.4d0
  REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0) , sqonethird=SQRT(1.d0/3.d0), &
       COSA=COS(pi/6), SINA=SIN(pi/6), x0=1/6.d0
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
    REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv
    REAL(KIND=8) :: xref, yref
    INTEGER :: n

    yref = 0.d0
    IF (inputs%time<1.d-8) THEN
       xref = x0
    ELSE 
       xref = x0+20*inputs%time*sqonethird
    END IF

    DO n = 1, SIZE(vv)
       IF (rr(1,n) .LE. xref + rr(2,n)*sqonethird) THEN
          vv(n) = 8.d0
       ELSE
          vv(n) = 1.4d0
       END IF
    END DO

  END FUNCTION rho_anal

  FUNCTION press_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
       REAL(KIND=8) :: xref, yref
    INTEGER :: n

    IF (inputs%time<1.d-8) THEN
       xref = x0
    ELSE 
       xref = x0+20*inputs%time*sqonethird
    END IF

    DO n = 1, SIZE(vv)
       IF (rr(1,n) .LE. xref + rr(2,n)*sqonethird)THEN
          vv(n) = 116.5
       ELSE
          vv(n) = 1.d0
       END IF
    END DO
  END FUNCTION press_anal

  FUNCTION vit_anal(comp,rr) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    REAL(KIND=8) :: xref, yref
    INTEGER :: n

    IF (inputs%time<1.d-8) THEN
       xref = x0
    ELSE 
       xref = x0+20*inputs%time*sqonethird
    END IF
 
    IF (comp==1) THEN
       DO n = 1, SIZE(vv)
          IF (rr(1,n) .LE. xref + rr(2,n)*sqonethird)THEN
             vv(n) = 8.25*COSA
          ELSE
             vv(n) = 0.d0
          END IF
       END DO
    ELSE
       DO n = 1, SIZE(vv)
          IF (rr(2,n).EQ.0.d0 .AND. rr(1,n).GE.x0) THEN
             vv(n) = 0.d0
          ELSE IF (rr(1,n) .LE. xref + rr(2,n)*sqonethird)THEN
             vv(n) = -8.25*SINA
          ELSE
             vv(n) = 0.d0
          END IF
       END DO
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
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(4,SIZE(rr,2)),   INTENT(OUT):: un
    un(1,:) = rho_anal(rr)
    un(2,:) = mt_anal(1,rr)
    un(3,:) = mt_anal(2,rr)
    un(4,:) = E_anal(rr)
  END SUBROUTINE init

  FUNCTION sol_anal(comp,rr,t)RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                                 INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),            INTENT(IN) :: rr
    REAL(KIND=8),                            INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
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
