MODULE shock_tub

CONTAINS

!  |                           i                                           |
!  o----o--------0-----O================O--------0----o--------------------o
!  |             n3    n1               n2       n4                        |


!  alternatively and equivalently


!  |                           i                                           |
!  o----o--------0-----O================O--------0----o--------------------o
!  |             n4    n2               n1       n3                        |




!  |  i                                                                    |
!  O====O--------0-----o----------------o--------o----o--------------------o
!n1=n3  n2       n3                                                        |



!  |                                                             i         |
!  o---o--------o-----o----------------o---------0----O====================O
!  |                                             n3   n1                 n2=n4

  SUBROUTINE sol_exa_stub(a,b,np,itest,t_final,x_e,v_e,r_e,p_e)

    USE  grid_1d  ! xx, cell_size, node_pair
    USE  ideal_polytropic_gas     ! Definition of gamma
    USE  exact_Riemann_solver_i   ! Solution of Riemann problem
    USE  shock_tube_profiles_i    ! Profile of solution in space

    !USE  Roe_linearization

    IMPLICIT NONE

    REAL(KIND=8) , INTENT(IN) :: a,b
    INTEGER,       INTENT(IN) :: itest
    REAL (KIND=8), INTENT(IN) :: t_final
    REAL(KIND=8), DIMENSION(:), POINTER :: x_e, v_e, r_e, p_e

    INTEGER, INTENT(IN) :: np 

    REAL (KIND=8), DIMENSION(:, :), ALLOCATABLE :: FF
    REAL (KIND=8), DIMENSION(:, :), ALLOCATABLE :: uu
    REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: rr_e,  vv_e,  PP_e,  ee_e,  PP

    REAL (KIND=8) :: c_l, c_r, r_l, r_r,  v_c, P_c, c_l_c, c_r_c, r_l_c, r_r_c

    REAL (KIND=8) :: rho_l, v_l, P_l, m_l,  Et_l, & ! Sod
         rho_r, v_r, P_r, m_r,  Et_r    ! Sod

!!!REAL (KIND=8) :: rho_l = 1,   v_l = 3,   P_l = 1,   m_l,   Et_l,  &  ! 2sw
    !   REAL (KIND=8) :: rho_l = 1,   v_l = -0.5,P_l = 1,   m_l,   Et_l,  &  ! 2sw
    !                    rho_r = 2,   v_r = 0.25, P_r = 1,   m_r,   Et_r      ! 2sw

    INTEGER, DIMENSION(3) :: limiters = 1

    CHARACTER, DIMENSION(2) :: l_r_w

    LOGICAL :: uniform = .TRUE.
    INTEGER :: n, n1, n2,  i,  k
    INTEGER :: k1 = 40,  k2 = 55 ! TEST unordered xx
    REAL(KIND=8) :: x            ! TEST unordered xx
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    CALL compute_gamma_stuff

    IF(ALLOCATED(rr_e)) THEN
       DEALLOCATE(rr_e,vv_e,PP_e,ee_e,PP,uu,x_e,v_e,r_e,p_e,FF)
    END IF

    ALLOCATE(rr_e(np),  vv_e(np),  PP_e(np),  ee_e(np),  PP(np), uu(3,np))

    SELECT CASE(itest)

    CASE(1)
       rho_l = 1;      v_l = 0;   P_l = 1;    & ! Sod
            rho_r = 0.125;  v_r = 0;   P_r = 0.1;    ! Sod
    CASE(2)
       rho_l= 0.445d0; rho_r= 0.5d0
       v_l = 0.698d0;   v_r = 0
       P_l = 3.528d0;   P_r = 0.571d0   !Lax 
    CASE(3)  ! CT + VD
       rho_l = 120;          rho_r = 1.2
       v_l = 0;              v_r = 0
       p_l = rho_l / gamma_;  p_r = rho_r / gamma_
    CASE(4) !Leblanc
       rho_l = 1.0d0;  rho_r = 0.001d0
       v_l   = 0.d0;   v_r   = 0.d0
       p_l = (gamma_-1.d0)*rho_l*1.d-1
       p_r = (gamma_-1.d0)*rho_r*1.d-7
    CASE(5) !Sonic
       rho_l=1.d0; rho_r=1.d0
       v_l = -1.d0; v_r =1.d0
       p_l =0.01d0; p_r =0.01d0
    CASE(6) !expa
       rho_l = 3.0
       rho_r = 0.5
       p_l   = 1.0
       p_r   = p_l*(rho_r/rho_l)**gamma_
       v_l   = 1.0
       v_r   = v_l + 2.d0*sqrt(gamma_*p_l/rho_l)/(gamma_-1.d0)&
            -2.d0*sqrt(gamma_*p_r/rho_r)/(gamma_-1.d0)   
    END SELECT


    ALLOCATE(x_e(np),v_e(np),r_e(np),p_e(np))

    CALL  gen_grid_1d (a, b, np,  uniform) ! last par is optional

    ALLOCATE ( FF(3, SIZE(node_pair,2)) )

    !  TEST unordered xx
    !  PERMUTATION of NODES

    !  INVERSION
    !   k1 = 47;  k2 = 48
    !   x = xx(k1);  xx(k1) = xx(k2);  xx(k2) = x
    !   node_pair(2,k1-1) = k2
    !   node_pair(1,k1  ) = k2
    !   node_pair(2,k2-1) = k1
    !   node_pair(1,k2  ) = k1
    !  TEST unordered xx


    ! initial conditions:  left  state ---> _l
    ! initial conditions:  right state ---> _r

    m_l = rho_l * v_l;   Et_l = P_l/gamma_m1_ + rho_l * v_l**2/2
    m_r = rho_r * v_r;   Et_r = P_r/gamma_m1_ + rho_r * v_r**2/2

    DO n = 1, np
       IF (xx(n) < (a+b)/2) THEN
          uu(1,n) = rho_l;   uu(2,n) = m_l;   uu(3,n) = Et_l
       ELSE
          uu(1,n) = rho_r;   uu(2,n) = m_r;   uu(3,n) = Et_r
       ENDIF
    ENDDO

    !   dt = 0.00411  ! with k_tot = 35  --->  Sod's test problem
    !   dt = 0.002;  k_tot = 1000 

    !  dt = 0.001;  k_tot = 70        ! --->  RLV test problem

    !   DO k = 1, k_tot   ! time loop
    !
    !      CALL  HR_Roe_Euler_1d (node_pair, xx, uu, dt,  FF,  limiters=limiters)
    !
    !      DO i = 1, SIZE(node_pair, 2)
    !
    !         n1 = node_pair(1,i);   n2 = node_pair(2,i)
    !         uu(:,n1) = uu(:,n1)  -  (dt/cell_size(n1)) * FF(:,i)
    !         uu(:,n2) = uu(:,n2)  +  (dt/cell_size(n2)) * FF(:,i)
    !
    !      ENDDO
    !
    !      ! NO INTERACTION OF WAVES WITH INTERVAL EXTREMES IS ACCOUNTED FOR
    !      uu(1, 1) = rho_l;  uu(2, 1) = m_l;  uu(3, 1) = Et_l  ! left  B.C.
    !      uu(1,np) = rho_r;  uu(2,np) = m_r;  uu(3,np) = Et_r  ! right B.C.
    !
    !   ENDDO

    ! exact solution for comparison

    c_l = SQRT(gamma_*P_l/rho_l)
    c_r = SQRT(gamma_*P_r/rho_r)

    CALL  exact_Riemann_solver (v_l, P_l, c_l,  v_r, P_r, c_r,   &
         v_c, P_c, c_l_c, c_r_c,  l_r_w = l_r_w)
    !Character variable with two 2 comp.
    !for left and right

    !WRITE (*,*) l_r_w

    r_l = rho_l;  r_r = rho_r

    r_l_c = gamma_ * P_c / c_l_c**2
    r_r_c = gamma_ * P_c / c_r_c**2

    CALL  Riemann_tube_profiles (r_l, v_l, P_l,   &
         r_l_c, r_r_c, v_c, P_c,  r_r, v_r, P_r,   &
         l_r_w,  (a+b)/2, t_final,  xx,           &
         rr_e,  vv_e,  PP_e,  ee_e)

    !   CALL plot_profile (xx, uu(1,:), rr_e, "density")
    !   CALL plot_profile (xx, uu(2,:)/uu(1,:), vv_e, "velocity")
    !   CALL plot_profile (xx, uu(3,:), ee_e, "energy_t")

    PP = gamma_m1_ * (uu(3,:) - 0.5*uu(2,:)**2/uu(1,:))
    CALL plot_profile (xx, PP, PP_e, "pressure")

    x_e = xx
    v_e = vv_e
    r_e = rr_e
    p_e = PP_e


  CONTAINS
    !=======

    !___________________________________________________________________________

    SUBROUTINE plot_profile (xx, uu, ue, file_name)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xx, uu, ue
      CHARACTER(*),               INTENT(IN) :: file_name

      INTEGER :: i
      REAL(KIND=8) :: x_span, x_min, x_max,  u_span, u_min, u_max

      OPEN (UNIT = 20, FILE = "HR_"//file_name, FORM = 'formatted', STATUS = 'unknown')

      x_span = MAXVAL(xx) - MINVAL(xx)
      x_min = MINVAL(xx) - x_span/10
      x_max = MAXVAL(xx) + x_span/10

      IF (file_name /= 'velocity') THEN
         u_span = MAXVAL(uu) - 0
         u_min = - u_span/10
      ELSE
         u_span = MAXVAL(uu) - MINVAL(uu)
         u_min  = MINVAL(uu) - u_span/10
      ENDIF
      u_max = MAX(MAXVAL(uu), 1.0d0) + u_span/10

      WRITE (20,*) '$ DATA = CURVE2D'
      WRITE (20,*) '% toplabel = "', file_name, '"'
      WRITE (20,*) '% xlabel = "x"'
      WRITE (20,*) '% ylabel = "', file_name, '"'

      IF (ABS(x_span - MAXVAL(uu)) < 0.1) WRITE (20,*) '% equalscale = true'

      WRITE (20,*) '% xmin = "', x_min, '"'
      WRITE (20,*) '% xmax = "', x_max, '"'
      WRITE (20,*) '% ymin = "', u_min, '"'
      WRITE (20,*) '% ymax = "', u_max, '"'

      WRITE (20,*) '% grid = true'


1000  FORMAT(1x, a18, i4, a11, d10.4, a7, d10.4, a1)

      WRITE (20,*) '% markertype = 2'
      IF ( ANY( xx(1:SIZE(xx)-1) > xx(2:SIZE(xx)) ) )  &
           WRITE (20,*) '% linetype = 0'

      DO i = 1, SIZE(xx)
         write (20,*) xx(i), uu(i)
      ENDDO


      WRITE (20,*)
      IF ( ANY( xx(1:SIZE(xx)-1) > xx(2:SIZE(xx)) ) )  THEN
         WRITE (20,*) '% markertype = 12'
         WRITE (20,*) '% linetype = 0'
      ELSE
         WRITE (20,*) '% markertype = 1'
      ENDIF

      DO i = 1, SIZE(xx)
         write (20,*) xx(i), ue(i)
      ENDDO

      CLOSE (20)

    END SUBROUTINE plot_profile

  END SUBROUTINE sol_exa_stub

END MODULE shock_tub
