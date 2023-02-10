MODULE boundary_conditions
  USE input_data
  USE space_dim
  PUBLIC :: flux, sol_anal, pressure, init, rho_anal, press_anal, mt_anal, E_anal
  REAL(KIND=8), PUBLIC :: gamma
  PRIVATE
  REAL(KIND=8) :: rhol, pl, ul
  REAL(KIND=8) :: rhor, pr, ur
  REAL(KIND=8) :: rhom, pm, um
  REAL(KIND=8) :: long, x0, x1
  REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0)
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
  SUBROUTINE init(un,rr)
    USE def_of_gamma
    USE lambda_module
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),                INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(k_dim+2,SIZE(rr,2)), INTENT(OUT):: un
    REAL(KIND=8) :: cl, cr
    SELECT CASE(inputs%type_test)
    CASE('laxt')
       gamma = 1.4d0
       long=1.d0
       x0 = long*0.5d0
       rhol = 0.445d0
       rhor = 0.5d0
       pl = 3.528d0
       pr = 0.571d0 
       ul = 0.698d0
       ur = 0.d0
    CASE ('sodt')
       gamma = 1.4d0
       long=1.d0
       x0 = long*0.5d0
       rhol=1.d0
       rhor=0.125d0
       pl = 1.d0
       pr = 0.1d0
       ul = 0.d0
       ur = 0.d0
    CASE ('cont')
       gamma = 1.4d0
       long=1.d0
       x0 = long*0.15d0 !long*0.15
       rhol=2.d0
       rhor=1.d0
       pl = 1.d0
       pr = 1.d0
       ul = 0.d0
       ur = 0.d0
    CASE ('lebl')
       gamma = 5.d0/3.d0
       long=1.d0
       x0 = long*0.33d0
       rhol = 1.d0 
       rhor = 0.001d0
       pl = (gamma-1.d0)*rhol*1.d-1
       pr = (gamma-1.d0)*rhor*1.d-7
       ul = 0.d0
       ur = 0.d0
    CASE ('soni')
       gamma = 1.4d0
       long=1.d0
       x0 = long*0.5d0
       rhol=3.857d0
       rhor=1.d0
       pr=1.d0
       pl=10.333d0
       ul=0.92d0
       ur=3.55d0
    CASE('expa')
       gamma = 1.4d0
       long = 1.d0
       x0   = 0.2d0*long
       rhol = 3.0
       rhor = 0.5
       !rhol = 0.003
       !rhor = 0.0005
       pl   = 1.0
       pr   = pl*(rhor/rhol)**gamma
       cl = SQRT(gamma*pl/rhol)
       cr = SQRT(gamma*pr/rhor)
       ul = cl
       ur   = ul + 2*sqrt(gamma*pl/rhol)/(gamma-1.d0)-2*sqrt(gamma*pr/rhor)/(gamma-1.d0)
       inputs%time = 0.2d0/(ur-cr)
    CASE('blas')
       gamma = 1.4d0
       long = 1.d0
       x0   = 0.1d0*long
       x1   = 0.9d0*long
       rhol = 1.d0
       rhor = 1.d0
       rhom = 1.d0
       pl = 1000.d0
       pm = .01d0 !I had a bug here pm = .1d0
       pr = 100.d0
       ul = 0.d0
       ur = 0.d0
       um = 0.d0
    CASE('smth')
       gamma = 1.4d0
       long = 4.d0
       x0=0.1!-1!0.1
       x1=0.3!0!0.3
    CASE('RPE4') !Qartapelle page 78
       gamma = 1.d0 +0.008d0
       long = 1.d0
       x0   = 0.25d0*long
       rhol = 1.d0/0.818
       rhor = 1.d0/1.089
       pl = 1.256
       pr = 1.028
       ul = 0.141
       ur = 1.894
    CASE('sltz') !Saltman problem
       gamma = 5.d0/3.d0
       long = 1.d0
       x0   = 0.05d0*long
       rhol = 1.d0
       rhor = 1.d0
       pl = 4.d-4/3.d0
       pr = 4.d-4/3.d0
       ul = 2.d0
       ur = 0.d0
    CASE DEFAULT
       WRITE(*,*) ' BUG in init'
       STOP
    END SELECT
    un(1,:) = rho_anal(rr)
    un(2,:) = mt_anal(1,rr)
    un(3,:) = E_anal(rr)
    CALL set_gamma_for_riemann_solver(gamma)
  END SUBROUTINE init
  
  FUNCTION rho_anal(rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv
    REAL(KIND=8) :: xi, cl, cr, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    SELECT CASE(inputs%type_test)
    CASE ('laxt','sodt','cont','soni','sltz')
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = rhol  
          ELSE
             vv(n) = rhor
          END IF
       END DO
    CASE('lebl')
       rhostarL=5.40793353493162d-2
       rhostarR=3.99999806043000d-3
       vstar   =0.621838671391735
       pstar   =0.515577927650970d-3
       lambda1 =0.495784895188979
       lambda3 =0.829118362533470
       DO n = 1, SIZE(vv)
          xi=(rr(1,n)-x0)/inputs%time
          IF (xi.LE. -1.d0/3) THEN
             vv(n) = rhol
          ELSE IF (xi.LE. lambda1) THEN
             vv(n) = (0.75d0-0.75d0*xi)**3
          ELSE IF (xi.LE. vstar) THEN
             vv(n) =rhostarL
          ELSE IF (xi.LE. lambda3) THEN
             vv(n) =rhostarR
          ELSE
             vv(n) = rhoR
          END IF
       END DO
    CASE('expa')
       cl = SQRT(gamma*pl/rhol)
       cr = SQRT(gamma*pr/rhor)
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0+inputs%time*(ul-cl)) THEN
             vv(n) = rhol
          ELSE IF (rr(1,n)<x0+inputs%time*(ur-cr)) THEN
             vv(n) = rhol*(2.d0/(gamma+1.d0) &
                  + ((gamma-1.d0)/((gamma+1.d0)*cl))*(ul-(rr(1,n)-x0)/inputs%time))**(2.d0/(gamma-1.d0))
          ELSE
             vv(n) = rhor
          END IF
       END DO
    CASE('blas')
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = rhol
          ELSE IF (rr(1,n)<x1) THEN
             vv(n) = rhom
          ELSE
             vv(n) = rhor
          END IF
       END DO
    CASE('smth')
       DO n = 1, SIZE(vv)
          IF ((rr(1,n)-inputs%time)<x0 .OR. (rr(1,n)-inputs%time)>x1) THEN
             vv(n) = 1.d0
          ELSE
             !vv(n) = 1 + 1*SIN((rr(1,n)-inputs%time-x0)*pi/(x1-x0))**3
             vv(n) = 1 + (2/(x1-x0))**6*(rr(1,n)-inputs%time-x0)**3*(x1-rr(1,n)+inputs%time)**3
          END IF
       END DO
    CASE('RPE4') !Qartapelle page 78
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
    REAL(KIND=8) :: xi, cl, cr, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    SELECT CASE(inputs%type_test)
    CASE ('laxt','sodt','cont','soni','sltz')
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = pl
          ELSE
             vv(n) = pr
          END IF
       END DO
    CASE('lebl')
       rhostarL=5.40793353493162d-2
       rhostarR=3.99999806043000d-3
       vstar   =0.621838671391735
       pstar   =0.515577927650970d-3
       lambda1 =0.495784895188979
       lambda3 =0.829118362533470
       DO n = 1, SIZE(vv)
          xi=(rr(1,n)-x0)/inputs%time
          IF (xi.LE. -1.d0/3) THEN
             vv(n) = pl
          ELSE IF (xi.LE. lambda1) THEN
             vv(n) = (0.75d0-0.75d0*xi)**5/15
          ELSE IF (xi.LE. lambda3) THEN
             vv(n) = pstar
          ELSE
             vv(n) = pr
          END IF
       END DO
    CASE('expa')
       cl = SQRT(gamma*pl/rhol)
       cr = SQRT(gamma*pr/rhor)
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0+inputs%time*(ul-cl)) THEN
             vv(n) = pl
          ELSE IF (rr(1,n)<x0+inputs%time*(ur-cr)) THEN
             vv(n) = pl*(2.d0/(gamma+1.d0) &
                  + ((gamma-1.d0)/((gamma+1.d0)*cl))*(ul-(rr(1,n)-x0)/inputs%time))**(2.d0*gamma/(gamma-1.d0))
          ELSE
             vv(n) = pr
          END IF
       END DO
    CASE('blas')
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = pl
          ELSE IF (rr(1,n)<x1) THEN
             vv(n) = pm
          ELSE
             vv(n) = pr
          END IF
       END DO
    CASE('smth')
       vv = 1.d0
    CASE('RPE4') !Qartapelle page 78
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
    REAL(KIND=8) :: xi, cl, cr, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    SELECT CASE(inputs%type_test)
    CASE ('laxt','sodt','cont','soni','sltz')
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = ul 
          ELSE
             vv(n) = ur
          END IF
       END DO
    CASE('lebl')
       rhostarL=5.40793353493162d-2
       rhostarR=3.99999806043000d-3
       vstar   =0.621838671391735
       pstar   =0.515577927650970d-3
       lambda1 =0.495784895188979
       lambda3 =0.829118362533470
       DO n = 1, SIZE(vv)
          xi=(rr(1,n)-x0)/inputs%time
          IF (xi.LE. -1.d0/3) THEN
             vv(n) = ul
          ELSE IF (xi.LE. lambda1) THEN
             vv(n) = 0.75d0*(1./3.d0 +xi)
          ELSE IF (xi.LE. lambda3) THEN
             vv(n) = vstar
          ELSE
             vv(n) = ur
          END IF
       END DO
    CASE('expa')
       cl = SQRT(gamma*pl/rhol)
       cr = SQRT(gamma*pr/rhor)
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0+inputs%time*(ul-cl)) THEN
             vv(n) = ul
          ELSE IF (rr(1,n)<x0+inputs%time*(ur-cr)) THEN
             vv(n) = (2.d0/(gamma+1))*(cl + ((gamma-1.d0)/2.d0)*ul +(rr(1,n)-x0)/inputs%time)
          ELSE
             vv(n) = ur
          END IF
       END DO
    CASE('blas')
       DO n = 1, SIZE(vv)
          IF (rr(1,n)<x0) THEN
             vv(n) = ul
          ELSE IF (rr(1,n)<x1) THEN
             vv(n) = um
          ELSE
             vv(n) = ur
          END IF
       END DO
    CASE('smth')
       vv = 1.d0
    CASE ('RPE4')
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
    IF (inputs%type_test=='RPE4') THEN
       vv =  (press_anal(rr) + 3*rho_anal(rr)**2)*(1-rho_anal(rr)/3)/(gamma-1.d0) -3*rho_anal(rr)**2 &
            + rho_anal(rr)*(vit_anal(1,rr)**2)/2
       RETURN
    END IF
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
    REAL(KIND=8), PARAMETER :: a=3.d0, b=1.d0/3.d0
    IF (inputs%type_test=='RPE4') THEN
       vv = (gamma-1.d0)*(un(3,:)-0.5d0*(un(2,:)**2)/un(1,:) + a*un(1,:)**2)/(1-b*un(1,:))&
            - a*un(1,:)**2
    ELSE
       vv=(un(3,:)-0.5d0*(un(2,:)**2)/un(1,:))*(gamma-1)
    END IF
  END FUNCTION pressure

!!$  FUNCTION sound_speed(un) RESULT(vv)
!!$    IMPLICIT NONE 
!!$    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
!!$    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: vv
!!$    vv= SQRT(gamma*pressure(un)/un(1,:))
!!$  END FUNCTION sound_speed

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
