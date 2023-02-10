MODULE boundary_conditions
  USE space_dim
  USE vdw
  PUBLIC :: flux, sol_anal, pressure, init, rho_anal, press_anal, mt_anal, E_anal
  INTEGER, PUBLIC      :: VdW_test_case
  REAL(KIND=8), PUBLIC :: gamma = 1.4d0 !===FIXME
  PRIVATE
  REAL(KIND=8) :: rhol, pl, ul
  REAL(KIND=8) :: rhor, pr, ur
  REAL(KIND=8) :: long, x0, x1
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
    USE lambda_module
    USE character_strings
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),                INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(k_dim+2,SIZE(rr,2)), INTENT(OUT):: un
    REAL(KIND=8) :: cl, cr, rho_plus
    REAL(KIND=8) :: in_state(2), in_data(3), out_state(4)
    INTEGER :: in_unit=10
    OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(in_unit, "===VdW test case: (0,2)===") 
    READ (in_unit,*) VdW_test_case
    CLOSE(in_unit)

    SELECT CASE(VdW_test_case)
    CASE(0) !===Convergence test Section 6.1, SISC Vol. 44, No. 1, pp. A444-A470
       rho_plus = 0.35d0
       rhol = 0.1d0
       rhor = 0.39d0
       in_state(1) = rhol
       in_state(2) = rhor
       in_data(1) = 1.d0
       in_data(2) = 1.d0
       in_data(3) = 1.02d0
       CALL initialize_vdw(rho_plus, in_state, in_data, out_state)
       uL = out_state(1)
       uR = out_state(2)
       pL = out_state(3)
       pR = out_state(4)
    CASE(1) !===Stability test 1, Eq. (6.4) section 6.2, SISC Vol. 44, No. 1, pp. A444-A470
       avdw = 1.d0
       bvdw = 1.d0
       gamma_vdw = 1.02d0
       rhol= (0.5d0)*(2-gamma_vdw)/(2*bvdw)  !===rho must be smaller than (2-gamma_vdw)/(2*bvdw)
       rhor= (0.25d0)*(2-gamma_vdw)/(2*bvdw)
       pl = 1.01*avdw*rhol**2*(2-gamma_vdw-2*bvdw*rhol)/(gamma_vdw)
       pr = 1.913*avdw*rhor**2*(2-gamma_vdw-2*bvdw*rhor)/(gamma_vdw)
       ul = 0.
       ur = 0.
    CASE(2) !===Stability test 1, Eq. (6.5) section 6.2, SISC Vol. 44, No. 1, pp. A444-A470
       avdw = 1.d0
       bvdw = 1.d0
       gamma_vdw = 1.02d0
       rhol = 2.5d-1
       rhor = 4.9d-5
       ul = 0
       ur = 0
       pl = 3.d-2
       pr = 5.d-8
    CASE(3) !===Stability test 1, Eq. (6.6) section 6.2, SISC Vol. 44, No. 1, pp. A444-A470
       avdw = 1.d0
       bvdw = 1.d0
       gamma_vdw = 1.02d0
       rhol=0.9932
       rhor = 0.95
       ul = 3
       ur = -3
       pl = 2
       pr = 2
    END SELECT
    WRITE(*,*) avdw, bvdw, gamma_vdw
    WRITE(*,*) 'rhol', rhol, 'rhor', rhor
    write(*,*) 'vl', ul, 'vr', ur
    write(*,*) 'pl', pl, 'pr', pr
    WRITE (*,*) 'cl', SQRT(gamma_vdw*(pl+avdw*rhol**2)/(rhol*(1-bvdw*rhol)) - 2*avdw*rhol)
    WRITE (*,*) 'cR', SQRT(gamma_vdw*(pr+avdw*rhor**2)/(rhor*(1-bvdw*rhor)) - 2*avdw*rhor)
    un(1,:) = rho_anal(rr)
    un(2,:) = mt_anal(1,rr)
    un(3,:) = E_anal(rr)

  END SUBROUTINE init

  FUNCTION rho_anal(rr) RESULT(vv)
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
    vv =  (press_anal(rr) + avdw*rho_anal(rr)**2)*(1-bvdw*rho_anal(rr))/(gamma_vdw-1.d0) &
         - avdw*rho_anal(rr)**2 + rho_anal(rr)*(vit_anal(1,rr)**2)/2
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
    vv = pressure_vdw(un)
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
