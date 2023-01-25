PROGRAM  high_resolution_solver_1d


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

   USE  grid_1d  ! xx, cell_size, node_pair

   USE  ideal_polytropic_gas     ! Definition of gamma
   USE  exact_Riemann_solver_i   ! Solution of Riemann problem
   USE  shock_tube_profiles_i    ! Profile of solution in space

   USE  Roe_linearization

   IMPLICIT NONE

   INTEGER, PARAMETER :: np = 100

   REAL (KIND=8), DIMENSION(:, :), ALLOCATABLE :: FF
   REAL (KIND=8), DIMENSION(3, np) :: uu
   REAL (KIND=8), DIMENSION(np)    :: rr_e,  vv_e,  PP_e,  ee_e,  PP

   REAL (KIND=8) :: a = 0,  b = 1,  dt,  &
                    c_l, c_r, r_l, r_r,  v_c, P_c, c_l_c, c_r_c, r_l_c, r_r_c

!  REAL (KIND=8) :: rho_l = 1,      v_l = 0,   P_l = 1,     m_l,   Et_l, & ! Sod
!                   rho_r = 0.125,  v_r = 0,   P_r = 0.1,   m_r,   Et_r    ! Sod
                    !                                              Et_l = 2.5
                    !                                              Et_r = 0.25

!!!REAL (KIND=8) :: rho_l = 1,   v_l = 3,   P_l = 1,   m_l,   Et_l,  &  ! 2sw
   REAL (KIND=8) :: rho_l = 1,   v_l = -0.5,P_l = 1,   m_l,   Et_l,  &  ! 2sw
                    rho_r = 2,   v_r = 0.25, P_r = 1,   m_r,   Et_r      ! 2sw

   INTEGER, DIMENSION(3) :: limiters = 1

   CHARACTER, DIMENSION(2) :: l_r_w

   LOGICAL :: uniform = .TRUE.
   INTEGER :: n, n1, n2,  i,  k, k_tot = 35
   INTEGER :: k1 = 40,  k2 = 55 ! TEST unordered xx
   REAL(KIND=8) :: x            ! TEST unordered xx
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

   dt = 0.00411  ! with k_tot = 35  --->  Sod's test problem
   dt = 0.002;  k_tot = 70

!  dt = 0.001;  k_tot = 70        ! --->  RLV test problem

   DO k = 1, k_tot   ! time loop

      CALL  HR_Roe_Euler_1d (node_pair, xx, uu, dt,  FF,  limiters=limiters)

      DO i = 1, SIZE(node_pair, 2)

         n1 = node_pair(1,i);   n2 = node_pair(2,i)
         uu(:,n1) = uu(:,n1)  -  (dt/cell_size(n1)) * FF(:,i)
         uu(:,n2) = uu(:,n2)  +  (dt/cell_size(n2)) * FF(:,i)

      ENDDO

      ! NO INTERACTION OF WAVES WITH INTERVAL EXTREMES IS ACCOUNTED FOR
      uu(1, 1) = rho_l;  uu(2, 1) = m_l;  uu(3, 1) = Et_l  ! left  B.C.
      uu(1,np) = rho_r;  uu(2,np) = m_r;  uu(3,np) = Et_r  ! right B.C.

   ENDDO

   ! exact solution for comparison

   c_l = SQRT(gamma_*P_l/rho_l)
   c_r = SQRT(gamma_*P_r/rho_r)

   CALL  exact_Riemann_solver (v_l, P_l, c_l,  v_r, P_r, c_r,   &
                               v_c, P_c, c_l_c, c_r_c,  l_r_w = l_r_w)
                                                      !Character variable with two 2 comp.
                                                      !for left and right

   WRITE (*,*) l_r_w

   r_l = rho_l;  r_r = rho_r

   r_l_c = gamma_ * P_c / c_l_c**2
   r_r_c = gamma_ * P_c / c_r_c**2

   CALL  Riemann_tube_profiles (r_l, v_l, P_l,   &
       r_l_c, r_r_c, v_c, P_c,  r_r, v_r, P_r,   &
       l_r_w,  (a+b)/2, dt*k_tot,  xx,           &
           rr_e,  vv_e,  PP_e,  ee_e)

   CALL plot_profile (xx, uu(1,:), rr_e, "density")
   CALL plot_profile (xx, uu(2,:)/uu(1,:), vv_e, "velocity")
   CALL plot_profile (xx, uu(3,:), ee_e, "energy_t")

   PP = gamma_m1_ * (uu(3,:) - 0.5*uu(2,:)**2/uu(1,:))
   CALL plot_profile (xx, PP, PP_e, "pressure")


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

   WRITE (20,1000) '% comment = " N = ', np,   &
                   '   dt/dx = ', dt/(x_span/(SIZE(xx)-1)),   &
                   '   T = ', dt*k_tot, '"'
   1000 FORMAT(1x, a18, i4, a11, d10.4, a7, d10.4, a1)

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


END PROGRAM  high_resolution_solver_1d


!//////////////////////////////////////////////////////////////////////////////


MODULE  Roe_linearization


CONTAINS
!=======

SUBROUTINE  HR_Roe_Euler_1d (node_pair, xx, uu, dt,  FF,  limiters)

!  Calculation of the Numerical Flux  FF  at the interface of each
!  node-pair from the two known node-pair states  uu(:,n1)  and  uu(:,n2).
!  Solution to the approximate Riemman problem by Roe linearization
!  for the compressible Euler equations for an ideal gas in 1d.

!  METHOD:    Ideal polytropic gas
!             Roe linearization
!             Harten-Hyman entropy fixing, a la LeVeque
!             limiter expressed in Rebay's form:
!             psi(du,du') = du * phi(du'/du),  where
!             du = centred  and  du' = upwinded variations.

!  NOTATION:  ht ---> TOTAL specific (per unit mass) enthalpy

!  node_pair(1:2, i),  i = 1, 2, ..., number_of_node-pairs,
!  are the node number of the first and second point of the node-pair

!  xx(n),  n = 1, 2, ..., number_of_nodes,  are the abscissas
!  of the nodal points

!  uu(1:3, n),  n = 1, 2, ..., number_of_nodes,  is (in 1d, 2d, 3d)
!  the one-dimensional array of the conservation variables:
!
!             uu(1)  --->          mass density (rho)
!             uu(2)  --->      momentum density ( m )
!             uu(3)  --->  TOTAL energy density (E^t)  per unit volume

!  dt  time step

!  FF(1:3, i),  i = 1, 2, ..., number_of_node-pairs,
!  is the numerical flux at the interface  i  of the node-pair

!  limiters(1:3) (OPTIONAL) defines the limiter used for each
!  scalar equation of the hyperbolic system, as follows:
!
!    DEFAULT      --->  van Leer
!
!  limiters == 0  --->  no limiter, first-order upwind
!
!  limiters == 1  --->  van Leer
!  limiters == 2  --->  minmod
!  limiters == 3  --->  superbee
!  limiters == 4  --->  Monotonized Central
!  limiters == 5  --->  van Albada
!


   USE  ideal_polytropic_gas


   IMPLICIT NONE

   INTEGER,       DIMENSION(:,:), INTENT(IN)  :: node_pair
   REAL (KIND=8), DIMENSION(:),   INTENT(IN)  :: xx
   REAL (KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu
   REAL (KIND=8),                 INTENT(IN)  :: dt
   REAL (KIND=8), DIMENSION(:,:), INTENT(OUT) :: FF
   INTEGER,       DIMENSION(:),   INTENT(IN), OPTIONAL :: limiters

   REAL (KIND=8), DIMENSION(3,3) :: R, L
   REAL (KIND=8), DIMENSION(3)   :: ul, ur, ucl, ucr, lambda, lambda_ef_N,   &
                                    alpha,  alpha_l,  alpha_r,  u_sl

   REAL (KIND=8) :: sl,  sr,               v,   ht,        c,     &
                    vl,  htl,   Pl,  cl,   vr,  htr,  Pr,  cr,    &
                    vcl, htcl,  Pcl, ccl,  vcr, htcr, Pcr, ccr,   &
                    el_N,  ecl_P,  ecr_N,  er_P,                  &
                    dx,  f_o_upw,  LW_corr,  psi,  a, b,          &
                    zero = 0,  half = 0.5,  one = 1,  two = 2,    &
                    Co,  Co_max,   c_g, cc_g, cv_g

   INTEGER :: nl, nr, nle, nre, i, p


   ! evaluation of variations of uu on each node-pair:
   ! a sequential ordering of the nodes is assumed, i.e.,
   !
   !   i = 1    --->    node_pair(1,1) = 1,   node_pair(2,1) = 2
   !   i = 2    --->    node_pair(1,2) = 2,   node_pair(2,1) = 3
   !   i = 3    --->    node_pair(1,3) = 3,   node_pair(2,3) = 4
   !   i = 4    --->    node_pair(1,4) = 4,   node_pair(2,4) = 5
   !   .....    --->    ...................   ..................
   !
   ! PROGRAM VALID ONLY IN 1D
   !=========================


   Co_max = 0

   ! evaluation of the numerical fluxes

   DO i = 1, SIZE(node_pair, 2)  !  i  --->  Interface

      nl  = node_pair(1,i);   nr  = node_pair(2,i)  ! pair nodes
      nle = node_pair(3,i);   nre = node_pair(4,i)  ! external nodes

      dx = xx(nr) - xx(nl)

      ul = uu(:, nl)
      ur = uu(:, nr)


      vl = ul(2)/ul(1)
      vr = ur(2)/ur(1)
     htl = gamma_*ul(3)/ul(1) - alpha_*vl**2         ! TOTAL specific (per
     htr = gamma_*ur(3)/ur(1) - alpha_*vr**2         ! unit mass) enthalpy
      Pl = gamma_m1_ * (ul(3) - 0.5*ul(2)**2/ul(1))
!!!   Pr = gamma_m1_ * (ur(3) - 0.5*ur(2)**2/ur(1))  !! not used
      cl = SQRT(gamma_m1_ * (htl - vl*vl/2))
      cr = SQRT(gamma_m1_ * (htr - vr*vr/2))

      ! centred contribution to flux:  left-biased asymmetric expression

      FF(1,i) = ul(2)                        ! <<<<
      FF(2,i) = ul(2)**2/ul(1) + Pl          ! <<<<
      FF(3,i) = (ul(3) + Pl) * ul(2)/ul(1)   ! <<<<


      ! Roe's averaging

      sl = SQRT(ul(1));    sr = SQRT(ur(1))

      v = ( vl * sl  +   vr * sr) / (sl + sr)
     ht = (htl * sl  +  htr * sr) / (sl + sr)  ! TOTAL enthalpy per unit mass

      c = SQRT(gamma_m1_ * (ht - v*v/2))

      Co_max = MAX(Co_max, (dt/ABS(dx))*(ABS(v) + c))
      IF (Co_max > 1) THEN;  WRITE (*,*) 'Co_max > 1 for i = ', i; STOP;  ENDIF

      ! eigenvalues and right/left eigenvectors,

      lambda(1) = v - c;   lambda(2) = v;    lambda(3) = v + c

      R(1,1) = 1;          R(1,2) = 1;       R(1,3) = 1
      R(2,1) =  v - c;     R(2,2) = v;       R(2,3) =  v + c
      R(3,1) = ht - c*v;   R(3,2) = v*v/2;   R(3,3) = ht + c*v

!!!   R(3,1) = v*v/2 - c*v + c*c/gamma_m1_  ! equivalent quadratic expression
!!!   R(3,3) = v*v/2 + c*v + c*c/gamma_m1_  ! equivalent quadratic expression

      R = R / (2*c)

       c_g =   c/gamma_m1_
      cc_g = c*c/gamma_m1_
      cv_g = c*v/gamma_m1_

      L(1,1) =  v*v/2 + cv_g;   L(1,2) = - v - c_g;   L(1,3) =  1
      L(2,1) = -v*v + 2*cc_g;   L(2,2) = 2*v;         L(2,3) = -2
      L(3,1) =  v*v/2 - cv_g;   L(3,2) = - v + c_g;   L(3,3) =  1

      L = (gamma_m1_/c) * L

      alpha = MATMUL (L, ur - ul)

      ! ENTROPY FIXING (Harten and Hyman a la LeVeque)
      ! states and TRUE characteristic speeds
      ! no fixing for the second wave (contact discontinuity)

      ! notation:    l  |   cl    contact discontinuity     cr  |   r
      !             ul  |  ucl    contact discontinuity    ucr  |  ur

      ucl = ul  +  R(:,1) * alpha(1)
      ucr = ucl +  R(:,2) * alpha(2)

      vcl = ucl(2)/ucl(1)
      vcr = ucr(2)/ucr(1)
     htcl = gamma_*ucl(3)/ucl(1) - alpha_*vcl**2
     htcr = gamma_*ucr(3)/ucr(1) - alpha_*vcr**2
      ccl = SQRT(gamma_m1_ * (htcl - vcl**2/2))
      ccr = SQRT(gamma_m1_ * (htcr - vcr**2/2))

      el_N = MIN(vl - cl, zero);   ecl_P = MAX(vcl - ccl, zero)
      er_P = MAX(vr + cr, zero);   ecr_N = MIN(vcr + ccr, zero)

      lambda_ef_N  =  MIN(lambda, zero)

      IF (ecl_P /= el_N)  &
      lambda_ef_N(1)  =  el_N * (ecl_P - lambda(1)) / (ecl_P -  el_N)

      IF (er_P /= ecr_N)  &
      lambda_ef_N(3)  = ecr_N * ( er_P - lambda(3)) / ( er_P - ecr_N)

      ! the first-order upwind contribution to the flux is accumulated at
      ! the end, together with the limited second-order contribution


      alpha_l = MATMUL(L, (uu(:,nle) - ul) * (dx/(xx(nle) - xx(nl) + 1.d-8)) )
      alpha_r = MATMUL(L, (uu(:,nre) - ur) * (dx/(xx(nre) - xx(nr) + 1.d-8)) )
      ! the small value 1.d-8 is included to avoid division by zero
      ! (of a zero quantity) at the nodes belonging to the boundary

      DO p = 1, 3

         f_o_upw  =  lambda_ef_N(p) * alpha(p)   ! <<<<

         Co  =  (dt/dx) * lambda(p)

         ! second-order contribution with limiter

         a = alpha(p)                                    ! centred variation
         b = (alpha_l(p) + alpha_r(p)) / 2     &         ! upwind  variation
           + (alpha_l(p) - alpha_r(p)) * SIGN(half, Co)


         IF (PRESENT(limiters)) THEN

            SELECT CASE(limiters(p))

               CASE(0);  ! no limiter, first-order upwind
                   psi = 0

               CASE(1);  ! van Leer
                   psi = (a*ABS(b) + ABS(a)*b)/(ABS(a) + ABS(b) + 1.0d-8)

               CASE(2);  ! minmod
                   psi = (SIGN(half,a) + SIGN(half,b)) * MIN(ABS(a), ABS(b))

               CASE(3);  ! superbee
                   psi = (SIGN(half,a) + SIGN(half,b))   &
                       *  MAX( MIN(ABS(a), 2*ABS(b)),  MIN(2*ABS(a), ABS(b)) )

               CASE(4);  ! Monotonized Central
                   psi = MAX( zero,  MIN((a+b)/2, 2*a, 2*b) )    &
                       + MIN( zero,  MAX((a+b)/2, 2*a, 2*b) )

               CASE(5);  ! van Albada
                   psi = a * b * (a + b)/(a**2 + b**2 + 1.0d-8)
                   ! ATTENTION: this limiter does NOT vanishes if a*b < 0
                   ! never tested

            END SELECT

         ELSE

            ! default limiter: van Leer
            psi = (a*ABS(b) + ABS(a)*b)/(ABS(a) + ABS(b) + 1.0d-8)

         ENDIF

         LW_corr  =  lambda(p) * psi * (SIGN(one,Co) - Co) / 2   ! <<<<

         FF(:,i)  =  FF(:,i)  +  R(:,p) * (f_o_upw + LW_corr)    ! <<<<

      ENDDO


   ENDDO

   WRITE (*,*) 'Co_max = ', Co_max


!CONTAINS
!=======


!-------------------------------------------------------------------------------



END SUBROUTINE  HR_Roe_Euler_1d


END MODULE  Roe_linearization
