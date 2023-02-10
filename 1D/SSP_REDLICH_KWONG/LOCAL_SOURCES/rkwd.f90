MODULE rdkw
  PUBLIC :: pressure_rdkw, internal_energy_in_space_rdkw, initialize_rdkw, Temperature_rdkw
  PRIVATE
  REAL(KIND=8) :: R !=0.4d0
  REAL(KIND=8) :: alpha != 0.5d0
  REAL(KIND=8) :: b  != 0.5d0
  REAL(KIND=8) :: cv != 1.5d0 !===CASE3
  REAL(KIND=8) :: r1 != 0.d0
  REAL(KIND=8) :: r2 !=-1.d0
CONTAINS

    SUBROUTINE initialize_rdkw(in_data)
    IMPLICIT NONE
    REAL(KIND=8) :: rhop
    REAL(KIND=8), DIMENSION(:) :: in_data
    R     = in_data(1)
    alpha = in_data(2)
    b     = in_data(3)
    cv    = in_data(4)
    r1    = in_data(5)
    r2    = in_data(6)
    write(*,*) R, alpha, b, cv, r1, r2
  END SUBROUTINE initialize_rdkw
    
  FUNCTION Temperature_rdkw(un) RESULT(vv)
    USE space_dim
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2)) :: vv, e
    REAL(KIND=8) :: pp, qq, ccos, theta, x
    REAL(KIND=8), PARAMETER:: pi = -ACOS(-1.d0)
    INTEGER :: n
    e = un(k_dim+2,:)/un(1,:) - 0.5*(un(2,:)/un(1,:))**2
    DO n = 1, SIZE(un,2)
       pp = -e(n)/cv
       qq = (3*alpha/(2*b*cv))*(1/(r1-r2))*log((1-b*r1*un(1,n))/(1-b*r2*un(1,n)))
       x = qq**2/4+pp**3/27
       IF (x>0.d0) THEN !===qq always negative
          IF (pp.LE.0.d0) THEN
             vv(n) = (SQRT(x)-qq/2)**(1.d0/3) + (-SQRT(x)-qq/2)**(1.d0/3)
          ELSE
             vv(n) = (SQRT(x)-qq/2)**(1.d0/3) - (SQRT(x)+qq/2)**(1.d0/3)
          END IF
          vv(n) = vv(n)**2
       ELSE
          theta = ATAN(SQRT(-x),-qq/2)
          ccos = COS(theta/3)
          ccos = MAX(ccos,COS((theta+2*pi)/3))
          ccos = MAX(ccos,COS((theta+4*pi)/3))
          vv(n) = 2*SQRT(-pp/3)*ccos
          vv(n) = vv(n)**2
       END IF
    END DO
  END FUNCTION Temperature_rdkw

  FUNCTION pressure_rdkw(un)  RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2)) :: vv, temp, v
    REAL(KIND=8) :: pp, qq
    INTEGER :: n
    temp =  Temperature_rdkw(un)
    DO n = 1, SIZE(temp)
       WRITE(10,*) n, temp(n)
    END DO
    v = 1/un(1,:)
    vv = R*temp/(v-b) - alpha/(SQRT(temp)*(v-b*r1)*(v-b*r2)) 
  END FUNCTION pressure_rdkw
  
  FUNCTION internal_energy_rdkw(un)  RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2)) :: vv, temp, v
    REAL(KIND=8) :: R, alpha, b, cv, r1, r2, pp, qq
    temp =  Temperature_rdkw(un)
    v  = 1/un(1,:)
    vv = cv*temp + 3*alpha/(2*b*SQRT(Temp))*(1/(r1-r2))*log((v-b*r1)/(v-b*r2))
  END FUNCTION internal_energy_rdkw

  FUNCTION internal_energy_in_space_rdkw(press,rho) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: press, rho
    REAL(KIND=8), DIMENSION(SIZE(press)) :: vv
    REAL(KIND=8) :: pp, qq, ccos, theta, x, v, temp
    REAL(KIND=8), PARAMETER:: pi = -ACOS(-1.d0)
    INTEGER :: n
    DO n = 1, SIZE(press)
       v = 1/rho(n)
       pp = -press(n)*(v-b)/R
       qq = -alpha*(v-b)/(R*(v-b*r1)*(v-b*r2))
       x = qq**2/4+pp**3/27
       IF (x>0.d0) THEN !===qq is negative here
          IF (pp.LE.0.d0) THEN
             temp = (SQRT(x)-qq/2)**(1.d0/3) + (-SQRT(x)-qq/2)**(1.d0/3)
          ELSE
             temp = (SQRT(x)-qq/2)**(1.d0/3) - (SQRT(x)+qq/2)**(1.d0/3)
          END IF
          temp = temp**2
       ELSE
          theta = ATAN(SQRT(-x),-qq/2)
          ccos = COS(theta/3)
          ccos = MAX(ccos,COS((theta+2*pi)/3))
          ccos = MAX(ccos,COS((theta+4*pi)/3))
          temp = 2*SQRT(-pp/3)*ccos
          temp = temp**2
       END IF
       !write(*,*) press(n), R*temp/(v-b) -alpha/(SQRT(Temp)*(v-b*r1)*(v-b*r2))
       vv(n) = cv*temp + (3*alpha/(2*b*SQRT(Temp)))*(1/(r1-r2))*log((v-b*r1)/(v-b*r2))
       !write(*,*) vv(n), temp, press(n)
    END DO
  END FUNCTION internal_energy_in_space_rdkw
END MODULE rdkw
