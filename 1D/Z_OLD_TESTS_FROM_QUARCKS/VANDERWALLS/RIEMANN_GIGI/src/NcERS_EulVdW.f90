
! =========================================================================
!   Description: Solves Riemann Problem for 1D Euler System with van der
!                Waals gas
!                
!                 W_t + F_x = 0
!
!                where  W = (v, u, p)  and 'v' is specific volume, 'u' is 
!                velocity and 'p' is pressure. F is the vector of fluxes 
!                of mass, momentum and energy.
!                Exact  solution  is  provided even in the case of loss of 
!                genuine nonlinearity. Reduced variables are assumed:
! 
!                 reduced specific volume =  v = V / 3b
!                 reduced velocity        =  u = U / SQRT(a / 9b)
!                 reduced pressure        =  p = P * 27*b^2 / a
!
!                V = dimensional specific volume. U = dimensional velocity.
!                P = dimensional pressure.  a, b = coefficients in Van der 
!                Waals equation of state. 
! 
!
!      Author: Marco Fossati
!      Department of Aerospace Engineering
!      Politecnico di Milano
!      Via La Masa 34, 20156 Milano, ITALY
!      e-mail: fossati@aero.polimi.it
!
!          Year:  2006, February
! =========================================================================

PROGRAM EulVdW_NCRsolver

   USE eigenvalues   ! lambda_Init, lambda, fan_Init, fan
                     ! SUBROUTINE lambda_Init (sign_, w_i_, t_)
                     ! FUNCTION lambda_Vector(v, p, u)   RESULT(l)
                     ! FUNCTION lambda(sign, v, p, u)   RESULT(l)
                     ! SUBROUTINE fan_Init(sign_, sigma_, vi_, ui_, pi_)
                     ! FUNCTION fan(v)   RESULT(l)	
  
   USE nonlinear_wave 
       ! FUNCTION Wave_Structure(sign, v, w_i, vI_i, vE_i)   RESULT(wave)
       ! FUNCTION Inflection_Points(sigma)   RESULT(vI)
       ! FUNCTION Absolute_Envelope(sigma)   RESULT(vE)
       ! FUNCTION Tangency_Point(p0, v0, v_L, v_R)   RESULT(v_tan)
       ! FUNCTION IsoT_Tangent_Points(v0, t0)  RESULT(a_r_v)
       ! FUNCTION Internal_Intersection_Points(Pi_i, nu_i, Pi, nu, v_L, v_R)   RESULT(N_int)
       ! FUNCTION IsoT_Internal_Intersection(v_L, v_R, t, vv)  RESULT(N_int)
       ! FUNCTION Shock_Speed(vv, uu)   RESULT(speed)
       ! SUBROUTINE Invert_Wave (wave, inverted_wave)
  
   USE plot_utils
  
       ! SUBROUTINE Plot_Isentropes (w_L, w_R, sigma_L, sigma_R)
       ! SUBROUTINE Plot_vp_Waves (wave_1, wave_3)
       ! FUNCTION Linear_Interpolation(gg, xx, x) RESULT(gx)
       ! SUBROUTINE Plot_Char_Field (wave_1, wave_3, tt,   vacuum)  
       ! SUBROUTINE Plot_Profile (xx, uu, dt, name)  
       ! SUBROUTINE Animate_Profiles (x, ww, t)
 
   USE problem_params ! SUBROUTINE Load_Param
                      ! SUBROUTINE Check_Data (sigma_L, sigma_R, sigma_star_delta)
                      ! FUNCTION  Inside_Coexistence(P, v)  RESULT(i_c)	
  
   USE vacuum_check
       ! FUNCTION Check_Vacuum(sigma_L, vI_L, vE_L, &
       !                       sigma_R, vI_R, vE_R, &
       !                       wave_1, wave_3)  RESULT(vacuum)
       !
       ! FUNCTION Vacuum_Velocity_Term(sigma_ith, vI_i, vE_i, &
       !                               w_i, sign, wave_i)  RESULT(u_vac)
  
  
   USE var_types ! TYPE nonlin_wave
                 ! TYPE key_points
  
  
   USE vdw_gas ! FUNCTION    p_vdW__vt(v, t)  RESULT(p)
               ! FUNCTION  D_p_vdW__vt(v, t)  RESULT(D_p)
               ! FUNCTION DD_p_vdW__vt(v, t)  RESULT(DD_p)
               ! 
               ! FUNCTION sigma__Pinu(Pi, nu)  RESULT(sigma)
               !
               ! FUNCTION t_vdW__pv(p, v)  RESULT(t)
               ! FUNCTION sigma_star(delta)  RESULT(sigma_star_delta)
               !
               ! FUNCTION   Infl(nu)   RESULT(i)
               ! FUNCTION D_Infl(nu)   RESULT(D_i)
               !
               ! SUBROUTINE Init_Pi_vdW__sigmanu (sigma_i_)
               ! FUNCTION   Pi_vdW__sigmanu(nu)  RESULT(Pi)
               ! FUNCTION  DPi_vdW__sigmanu(nu)  RESULT(DPi)
               ! FUNCTION DDPi_vdW__sigmanu(nu)  RESULT(DDPi)
               !
               ! SUBROUTINE Init_psi (nu0_, sigma00_)
               ! FUNCTION  psi(nu)  RESULT(ps)
               ! FUNCTION Dpsi(nu)  RESULT(Dps)
               !
               ! SUBROUTINE Init_phi (Pi_i_, nu_i_, Pi_hat_, nu_hat_)
               ! FUNCTION  phi(nu)  RESULT(p)
               ! FUNCTION Dphi(nu)  RESULT(Dp)
               !
               ! FUNCTION  h(r)  RESULT(h_)
               ! FUNCTION Dh(r)  RESULT(Dh_)
               !
               ! FUNCTION  Q(v, v0,p0)  RESULT(P)
               ! FUNCTION DQ(v, v0,p0)  RESULT(dP)
               !
               ! SUBROUTINE  Init_IC (pressure, specific_volume)
               !
               ! FUNCTION SQRT_DQ(v, dum)  RESULT(f)
               !
               ! FUNCTION Gzero_Locus(v)  RESULT(p_G0)
               





   !-------------------------------------------------------------------------  
   IMPLICIT NONE
 
   TYPE (nonlin_wave) :: wave_1, wave_3, inverted_wave_3

   REAL(KIND=8) :: v_L, u_L, p_L, sigma_L,     &
                   v_R, u_R, p_R, sigma_R,     &
                   Du_iL, Du_iR, Dp_iL, Dp_iR, &
                   v, w, dv, dw, u, p,         &
                   a11, a12, a21, a22, det,    &
                   u_iL, u_iR, p_iL, p_iR,     &
                   sigma_star_delta, wv, wp,             &
                   v_Ml, u_Ml, p_Ml,           &
                   v_Mr, u_Mr, p_Mr,           &
                   v_a, v_b, xl, xr, dx,       &
                   shock_sigma, shock_posit,   &
                   v0, u0, p0, left_bound,     &
                   right_bound, contact_posit, c_L, c_R, sigma_P

   REAL(KIND=8), DIMENSION(2) :: vI_L, vE_L, &
                                 vI_R, vE_R, &
                                 vE_I


   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ww

   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: xx
  
   CHARACTER(LEN=1) :: wt

   INTEGER :: i, j, k, ws, p_i, np, n_shock, n_fan

   LOGICAL :: vacuum = .FALSE.
   !-------------------------------------------------------------------------

   PRINT*, ''
   PRINT*, ' ------ Riemann Solver for Euler System -----------------------------------'
   PRINT*, '        (van der Waals thermodynamics)'
   PRINT*, ''  

   OPEN(UNIT=12, FILE='verbose.out', FORM='formatted', STATUS='unknown')
  
   CALL Load_Param


  ! Local modification to accept as input the following data:----------
  !  v (specific volume) -> w_L.R(1)
  !  Mach                -> w_L.R(2)
  !  Pressure            -> w_L.R(3)

!   v_L = w_L(1);    v_R = w_R(1)
!   p_L = w_L(3);    p_R = w_R(3)
! 
!   sigma_L = sigma__Pinu(p_L/27, v_L*3)
!   sigma_R = sigma__Pinu(p_R/27, v_R*3)
! 
!   CALL Init_Pi_vdW__sigmanu(sigma_L)
!   c_L = SQRT(-v_L**2*27*DPi_vdW__sigmanu(v_L*3))
! 
!   CALL Init_Pi_vdW__sigmanu(sigma_R)
!   c_R = SQRT(-v_R**2*27*DPi_vdW__sigmanu(v_R*3))
! 
!   u_L = c_L
!   u_R = c_R
  ! END of modification -----------------------------------------------
  


  ! Initializing Variables --------------------------------------------
  ! uncomment following lines to resume from modfication above

   v_L = w_L(1);    v_R = w_R(1)
   u_L = w_L(2);    u_R = w_R(2)
   p_L = w_L(3);    p_R = w_R(3)

   sigma_L = sigma__Pinu(p_L/27, 3*v_L)
   sigma_R = sigma__Pinu(p_R/27, 3*v_R)   

   sigma_star_delta = sigma_star(delta)

   ! Compute the values of  v  for the Inflections Points 
   ! and for the Absolute Envelope when the isentrope 
   ! through the Left thermodynamic state is not convex
   IF (sigma_L < sigma_star_delta) THEN
      vI_L = Inflection_Points(sigma_L)
      vE_L = Absolute_Envelope(sigma_L)
   ENDIF

   ! Compute the values of  v  for the Inflections Points 
   ! and for the Absolute Envelope when the isentrope 
   ! through the Right thermodynamic state is not convex
   IF (sigma_R < sigma_star_delta) THEN
      vI_R = Inflection_Points(sigma_R)
      vE_R = Absolute_Envelope(sigma_R)
   ENDIF

!   DO i = 1, 55
!     sigma_P = 0.3151 - (i-1)*0.0003036
!     CALL Init_Pi_vdW__sigmanu(sigma_P)
!     vE_I = Absolute_Envelope(sigma_P)
!     print*, vE_I(1), 27*Pi_vdW__sigmanu(3*vE_I(1))
!   ENDDO
! 
!   DO i = 1, 55
!      sigma_P = 0.3151 - (i-1)*0.0003036
!      CALL Init_Pi_vdW__sigmanu(sigma_P)
!      vE_I = Absolute_Envelope(sigma_P)
!      print*, vE_I(2), 27*Pi_vdW__sigmanu(3*vE_I(2))
!   ENDDO


   CALL Check_Data (sigma_L, sigma_R, sigma_star_delta)
   CALL Plot_Isentropes (w_L, w_R, sigma_L, sigma_R)


   ! Check vacuum formation         
   PRINT*, ''          
   PRINT*, '    - Check Vacuum'        
   PRINT*, ''          
   vacuum = Check_Vacuum(sigma_L, vI_L, vE_L, &
                         sigma_R, vI_R, vE_R, &
                         wave_1, wave_3) 
    
   IF (vacuum) THEN

      ! Reversing order of component waves succession    
      CALL Invert_Wave     (wave_3, inverted_wave_3)       
      CALL Plot_Char_Field (wave_1, inverted_wave_3, time,   vacuum)

   ELSE

      PRINT*, ''
      PRINT*, '  - Newton method for Intermediate States...'
 
      v = (v_L + v_R) / 2
      w = (v_L + v_R) / 2

      DO i = 1, Max_Iter

         WRITE(12,*) ''                    
         WRITE(12,*) 'ITERATION', i    

         ! 1-Wave ------------------------------------------------------------
         WRITE(12,*) ''
         WRITE(12,*) 'Wave 1 ------------------------------------------------'

         wave_1 = Wave_Structure(-1, v, (/v_L, u_L, p_L/), vI_L, vE_L)

         ! Informations about the component of the composite wave
         ! reaching intermediate state
         ws = wave_1 % sign
         wt = wave_1 % type(SIZE(wave_1 % type))

         wv = wave_1 % v(SIZE(wave_1 % type), 1)
         wp = wave_1 % p(SIZE(wave_1 % type), 1)
 
         Du_iL = Du(ws, wt, v, wv, wp)
         Dp_iL = Dp(    wt, v, wv, wp)
         ! -------------------------------------------------------------------


         ! 3-Wave ------------------------------------------------------------
         WRITE(12,*) ''
         WRITE(12,*) 'Wave 3 ------------------------------------------------'

         wave_3 = Wave_Structure(+1, w, (/v_R, u_R, p_R/), vI_R, vE_R)

         ! Informations about the component of the composite wave
         ! reaching intermediate state
         ws = wave_3 % sign
         wt = wave_3 % type(SIZE(wave_3 % type))

         wv = wave_3 % v(SIZE(wave_3 % type), 1)
         wp = wave_3 % p(SIZE(wave_3 % type), 1)

         Du_iR = Du(ws, wt, w, wv, wp)
         Dp_iR = Dp(    wt, w, wv, wp)
         ! -------------------------------------------------------------------


         u_iL = wave_1 % u(SIZE(wave_1 % type), 2);   u_iR = wave_3 % u(SIZE(wave_3 % type), 2)
         p_iL = wave_1 % p(SIZE(wave_1 % type), 2);   p_iR = wave_3 % p(SIZE(wave_3 % type), 2)


         IF (Inside_Coexistence(p_iL, v)) THEN
            PRINT*, ''
            PRINT*, 'Error. The state (p_iL, v) falls under the coexistence curve'
            PRINT*, ''
            OPEN(UNIT=idf, FILE='stopped')
            WRITE(idf,*) w_L
            WRITE(idf,*) w_R
            WRITE (idf,*) 'The state (p_iL, v) falls under the coexistence curve'
            CLOSE(idf)
            STOP
         ENDIF
 
         IF (Inside_Coexistence(p_iR, w)) THEN
            PRINT*, ''
            PRINT*, 'Error. The state (p_iR, w) falls under the coexistence curve'
            PRINT*, ''
            OPEN(UNIT=idf, FILE='stopped')
            WRITE(idf,*) w_L
            WRITE(idf,*) w_R
            WRITE (idf,*) 'The state (p_iR, w) falls under the coexistence curve'
            CLOSE(idf)
            PRINT*, ''
            STOP
         ENDIF


         u = u_iL - u_iR ! Right hand side of incremental Newton
         p = p_iL - p_iR
 

         a11 = Du_iL;   a12 = - Du_iR

         a21 = Dp_iL;   a22 = - Dp_iR


         det = a11*a22 - a12*a21


         IF (det == 0) THEN

            PRINT*, ''
            PRINT*, 'Error. Null determinant in Newton Method for'
            PRINT*, 'the two by two system of the Riemann problem.'
            PRINT*, ''
            OPEN(UNIT=idf, FILE='stopped')
            WRITE(idf,*) w_L
            WRITE(idf,*) w_R
            WRITE (idf,*) 'Null determinant in Newton Method for the 2X2 Riemann Problem'
            CLOSE(idf)
            STOP

         ENDIF

 
         dv = - (a22*u - a12*p)/det
         dw = - (a11*p - a21*u)/det

         ! Relaxing increment in Newton method for stability -------
         v = v + relax * dv
         w = w + relax * dw
         !----------------------------------------------------------
      

         IF (v < 1.d0/3  .OR.  w < 1.d0/3) THEN

            PRINT*, '' 
            PRINT*, 'Error. Value of v is less than 1/3.'
            PRINT*, ''
            OPEN(UNIT=idf, FILE='stopped')
            WRITE(idf,*) w_L
            WRITE(idf,*) w_R
            WRITE (idf,*) 'Specific Volume less than 0.333.'
            CLOSE(idf)
            STOP

         ENDIF


         IF (ABS(dv) < ABS(v) * eps  .AND.  ABS(dw) < ABS(w) * eps) THEN

            v_Ml = v;  v_Mr = w
            u_Ml = u_iL;  u_Mr = u_iR
            p_Ml = p_iL;  p_Mr = p_iR

            PRINT*, '    ...method has converged in', i, 'iterations.'
            PRINT*, ''
 
            EXIT

         ENDIF

         IF (i == Max_Iter) THEN

            PRINT*, '' 
            PRINT*, 'Error. Maximum number of iterations reached.'
            PRINT*, 'Newton roots', v, w
            PRINT*, ''
            OPEN(UNIT=idf, FILE='stopped')
            WRITE(idf,*) w_L
            WRITE(idf,*) w_R
            WRITE (idf,*) 'Newton method for intermediate state failed'
            CLOSE(idf)
            STOP

         ENDIF

      ENDDO   ! Fine Metodo di Newton


      CALL Init_Wave_Sigma (wave_1)
      CALL Init_Wave_Sigma (wave_3)

      CALL Plot_vp_Waves (wave_1, wave_3)

      ! Reversing order of component waves succession
      ALLOCATE ( inverted_wave_3 % type(SIZE(wave_3 % type)),       &

                 inverted_wave_3 % v(SIZE(wave_3 % type), 2),       &
                 inverted_wave_3 % u(SIZE(wave_3 % type), 2),       &
                 inverted_wave_3 % p(SIZE(wave_3 % type), 2),       &

                 inverted_wave_3 % speeds(SIZE(wave_3 % type), 2),  &
                 inverted_wave_3 % sigma (SIZE(wave_3 % type), 2) )

      CALL Invert_Wave (wave_3, inverted_wave_3)

      CALL Plot_Char_Field (wave_1, inverted_wave_3, time)

 

      ! Identifying classical and nonclassical wave type
      ALLOCATE (type_1(SIZE(wave_1 % type)),  &
                type_3(SIZE(inverted_wave_3 % type)))

      CALL Classical_NonClassical (wave_1, inverted_wave_3,  type_1, type_3)



      PRINT*, '  - Resulting waves   ', type_1, '  ***  ', type_3
      PRINT*, ''
      PRINT*, '           LEFT         INTERMEDIATE-L    INTERMEDIATE-R        RIGHT'
      PRINT*, ''
      WRITE(*, 1000) 'v', v_L, v_Ml, v_Mr, v_R
      WRITE(*, 1000) 'u', u_L, u_Ml, u_Mr, u_R
      WRITE(*, 1000) 'p', p_L, p_Ml, p_Mr, p_R
      WRITE(*, 1000) 's', sigma_L,                                &
                          wave_1 % sigma(SIZE(wave_1 % type), 2), &
                          wave_3 % sigma(SIZE(wave_3 % type), 2), &
                          sigma_R

      1000 FORMAT(3x, a1, 1x, e13.6, 6x, e19.12, 5x, e19.12, 5x, e13.6)


      SELECT CASE (SIZE(wave_1 % type))

         CASE (2)

         IF (type_1(2) == 'F'  .OR.  &
             type_1(2) == 'f')  WRITE(*,*) 'Dsigma LEFT =', sigma_L - wave_1 % sigma(2,2)

         CASE (3)

         IF (type_1(2) == 'F')  WRITE(*,*) 'Dsigma LEFT =', sigma_L - wave_1 % sigma(2,2)
         IF (type_1(3) == 'f')  WRITE(*,*) 'Dsigma LEFT =', sigma_L - wave_1 % sigma(3,2)

      END SELECT


      SELECT CASE (SIZE(wave_3 % type))

         CASE (2)

         IF (type_3(1) == 'f'  .OR.  &
             type_3(1) == 'F')  WRITE(*,*) '                                             ' &
                                             //'Dsigma RIGHT =', sigma_R - wave_3 % sigma(2,2)
         CASE (3)

         IF (type_3(1) == 'f')  WRITE(*,*) '                                             ' &
                                             //'Dsigma RIGHT =', sigma_R - wave_3 % sigma(2,2)
         IF (type_3(2) == 'F')  WRITE(*,*) '                                             ' &
                                             //'Dsigma RIGHT =', sigma_R - wave_3 % sigma(3,2)
  
      END SELECT

      1001 FORMAT(17x, e15.8)
      1002 FORMAT(17x, e15.8)

      PRINT*, ''
      PRINT*, '   Componentwise Dsigma LEFT'
      DO i = 1, SIZE(wave_1 % type)
         WRITE(*,1003) wave_1 % sigma(i,2) - wave_1 % sigma(i,1)
      ENDDO

      PRINT*, ''
      PRINT*, '   Componentwise Dsigma RIGHT '
      DO i = 1, SIZE(wave_3 % type)
         WRITE(*,1003) wave_3 % sigma(i,2) - wave_3 % sigma(i,1)
      ENDDO

      1003 FORMAT(4x, e15.8)

! print*, 'wave_3'
! do i = 1, SIZE(wave_3%type)
! print*, wave_3%type(i)
! print*, wave_3%v(i,:)
! print*, wave_3%u(i,:)
! print*, wave_3%p(i,:)
! print*, wave_3%speeds(i,:)
! print*, wave_3%sigma(i,:)
! enddo
! print*, ''
! print*, 'inverted_wave_3'
! do i = 1, SIZE(wave_3%type)
! print*, inverted_wave_3%type(i)
! print*, inverted_wave_3%v(i,:)
! print*, inverted_wave_3%u(i,:)
! print*, inverted_wave_3%p(i,:)
! print*, inverted_wave_3%speeds(i,:)
! print*, inverted_wave_3%sigma(i,:)
! enddo


      ! The relevant computation stops here. The following instructions are
      ! performed to produce animation of the profiles.

      ! Animate Solution Profiles
      n_shock = COUNT(wave_1 % type == 'S')
      n_shock = n_shock + COUNT(inverted_wave_3 % type == 'S')

      n_fan = COUNT(wave_1 % type == 'F')
      n_fan = n_fan + COUNT(inverted_wave_3 % type == 'F')
      
      ! Number of grid points = 2 for left and right states, 2 for contact discontinuity,
      ! 2*n_shock for shocks and fan_points*n_fan for fans
      np = 4 + 2*n_shock + fan_points*n_fan


      ALLOCATE (xx(np))

      xl = wave_1 % speeds(1, 1) * time  &
         - ABS(0.3d0 * wave_1 % speeds(1, 1) * time)

      IF (xl > -0.5) xl = -0.5d0

      xr = inverted_wave_3 % speeds(SIZE(type_3), 2) * time  &
         + ABS(0.3 * inverted_wave_3 % speeds(SIZE(type_3), 2) * time)

      IF (xr < 0.5) xr = 0.5d0

      xx(1) = xl;   xx(np) = xr

      ALLOCATE (ww(3, np))

      ww(1,1) = w_L(1);   ww(1,np) = w_R(1)
      ww(2,1) = w_L(2);   ww(2,np) = w_R(2)
      ww(3,1) = w_L(3);   ww(3,np) = w_R(3)

      p_i = 1
   
      ! Account for for 1-Wave --------------------------------------------------------------------------
      DO i = 1, SIZE(wave_1 % type)

         ! Shock Wave
         IF (wave_1 % type(i) == 'S') THEN

            shock_sigma = wave_1 % speeds(i,1)
            shock_posit = shock_sigma * time

            xx(p_i + 1) = shock_posit
            xx(p_i + 2) = shock_posit

            ww(1, p_i + 1) = wave_1 % v(i,1);     ww(1, p_i + 2) = wave_1 % v(i,2)
            ww(2, p_i + 1) = wave_1 % u(i,1);     ww(2, p_i + 2) = wave_1 % u(i,2)
            ww(3, p_i + 1) = wave_1 % p(i,1);     ww(3, p_i + 2) = wave_1 % p(i,2)

            p_i = p_i + 2

         ENDIF


        ! Fan Wave
        IF (wave_1 % type(i) == 'F') THEN

           left_bound  = wave_1 % speeds(i,1) * time
           right_bound = wave_1 % speeds(i,2) * time

           xx(p_i + 1)   = left_bound
           xx(p_i + fan_points) = right_bound

           dx = (right_bound - left_bound) / (fan_points - 1)

           DO j = 1, fan_points - 2
              xx(p_i + 1+j) = xx(p_i + 1) + j*dx
           ENDDO

           v0 = wave_1 % v(i,1)
           u0 = wave_1 % u(i,1)
           p0 = wave_1 % p(i,1)

           v_a = wave_1 % v(i,1)
           v_b = wave_1 % v(i,2)

           ww(1, p_i + 1) = wave_1 % v(i,1)
           ww(2, p_i + 1) = u0 + Integral_Curve(wave_1 % sign, v0, ww(1, p_i + 1), v0, p0)
           ww(3, p_i + 1) = Q(ww(1, p_i + 1), v0, p0)


           DO k = p_i + 2, p_i + fan_points-1

              CALL fan_Init(wave_1 % sign, wave_1 % sigma(i,1), v0, u0, p0)

              ww(1, k) = Bisection(fan, v_a, v_b, eps, xx(k)/time)

              ww(2, k) = u0 + Integral_Curve(wave_1 % sign, v0, ww(1, k), v0, p0)

              ww(3, k) = Q(ww(1, k), v0, p0)

           ENDDO

           ww(1, p_i + fan_points) = wave_1 % v(i,2)
           ww(2, p_i + fan_points) = u0 + Integral_Curve(wave_1 % sign, v0, ww(1, p_i + fan_points), v0, p0)
           ww(3, p_i + fan_points) = Q(ww(1, p_i + fan_points), v0, p0)

           p_i = p_i + fan_points

        ENDIF

    ENDDO
    !-------------------------------------------------------------------------------------------------

  

    ! Account for 2-Wave -----------------------------------------------------------------------------
    contact_posit = wave_1 % u(SIZE(type_1),2) * time

    xx(p_i + 1) = contact_posit
    xx(p_i + 2) = contact_posit

    ww(1, p_i + 1) = wave_1 % v(SIZE(type_1),2);   ww(1, p_i + 2) = inverted_wave_3 % v(1,1)
    ww(2, p_i + 1) = wave_1 % u(SIZE(type_1),2);   ww(2, p_i + 2) = inverted_wave_3 % u(1,1)
    ww(3, p_i + 1) = wave_1 % p(SIZE(type_1),2);   ww(3, p_i + 2) = inverted_wave_3 % p(1,1)

    p_i = p_i + 2
    !-------------------------------------------------------------------------------------------------



    ! Account for for 3-Wave --------------------------------------------------------------------------
    DO i = 1, SIZE(inverted_wave_3 % type)
 
       ! Shock Wave
       IF (inverted_wave_3 % type(i) == 'S') THEN
 
          shock_sigma = inverted_wave_3 % speeds(i,1)
          shock_posit = shock_sigma * time

          xx(p_i + 1) = shock_posit
          xx(p_i + 2) = shock_posit

          ww(1, p_i + 1) = inverted_wave_3 % v(i,1);    ww(1, p_i + 2) = inverted_wave_3 % v(i,2)
          ww(2, p_i + 1) = inverted_wave_3 % u(i,1);    ww(2, p_i + 2) = inverted_wave_3 % u(i,2)
          ww(3, p_i + 1) = inverted_wave_3 % p(i,1);    ww(3, p_i + 2) = inverted_wave_3 % p(i,2)

          p_i = p_i + 2

       ENDIF
 
 
       ! Fan Wave
       IF (inverted_wave_3 % type(i) == 'F') THEN
 
          left_bound  = inverted_wave_3 % speeds(i,1) * time
          right_bound = inverted_wave_3 % speeds(i,2) * time

          xx(p_i + 1)   = left_bound
          xx(p_i + fan_points) = right_bound


          dx = (right_bound - left_bound) / (fan_points - 1)

          DO j = 1, fan_points - 2
             xx(p_i + 1+j) = xx(p_i + 1) + j*dx
          ENDDO

            v0 = inverted_wave_3 % v(i,1)
            u0 = inverted_wave_3 % u(i,1)
            p0 = inverted_wave_3 % p(i,1)

            v_a = inverted_wave_3 % v(i,1)
            v_b = inverted_wave_3 % v(i,2)

            ww(1, p_i + 1) = inverted_wave_3 % v(i,1)
            ww(2, p_i + 1) = u0 + Integral_Curve(wave_1 % sign, v0, ww(1, p_i + 1), v0, p0)
            ww(3, p_i + 1) = Q(ww(1, p_i + 1), v0, p0)


            DO k = p_i + 2,  p_i + fan_points - 1

               CALL fan_Init (inverted_wave_3 % sign, inverted_wave_3 % sigma(i,1), v0, u0, p0)
 
               ww(1, k) = Bisection(fan, v_a, v_b, eps, xx(k)/time)

               ww(2, k) = u0 + Integral_Curve(inverted_wave_3 % sign, v0, ww(1, k), v0, p0)

               ww(3, k) = Q(ww(1, k), v0, p0)

            ENDDO

            ww(1, p_i + fan_points) = inverted_wave_3 % v(i,2)
            ww(2, p_i + fan_points) = u0 + Integral_Curve(inverted_wave_3 % sign, v0, ww(1, p_i + fan_points), v0, p0)
            ww(3, p_i + fan_points) = Q(ww(1, p_i + fan_points), v0, p0)

            p_i = p_i + fan_points

         ENDIF

      ENDDO
      !-------------------------------------------------------------------------------------------------

      CALL Plot_Profile (xx, ww(1,:),   time, '1_OVER_RHO          ')
      CALL Plot_Profile (xx, 1/ww(1,:), time, 'DENSITY             ')
      CALL Plot_Profile (xx, ww(2,:),   time, 'VELOCITY            ')
      CALL Plot_Profile (xx, ww(3,:),   time, 'PRESSURE            ')

      CALL Animate_Profiles (xx, ww, time)

   ENDIF   ! Not Vacuum 

   CLOSE (12)

   PRINT*, ''
   PRINT*, ' ---------------- End of Computation --------------------------------------'
   PRINT*, ''


   ! Accessories for run with a set of pseudo-random initial data
   IF ( ANY(type_1 == 'F')   .OR.   ANY(type_1 == 'S')  .OR.   &
        ANY(type_3 == 'F')   .OR.   ANY(type_3 == 'S') ) THEN
  
      OPEN(UNIT = idf,  FILE = 'nonconvex')
         WRITE(idf,*) w_L
         WRITE(idf,*) w_R
         WRITE(idf,*) '  ', type_1, ' *** ', type_3
      CLOSE (idf)

   ELSE

      OPEN(UNIT = idf,  FILE = 'convex')        
         WRITE(idf,*) w_L
         WRITE(idf,*) w_R
         WRITE(idf,*) '  ', type_1, ' *** ', type_3
      CLOSE (idf)

   ENDIF 


CONTAINS


FUNCTION Du(wsign, type, v, v0, p0)   RESULT(D_u)

   ! Computes the first derivative for the solution of the p-th composite
   ! wave with respect the indipendent variable 'v'.  The dependence from
   ! 'v' appears only in the solution for the last component of the wave,
   ! i.e., the component that reaches the intermediate state.   

   !------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER,          INTENT(IN) :: wsign
   CHARACTER(LEN=1), INTENT(IN) :: type
   REAL(KIND=8),     INTENT(IN) :: v, v0, p0

   REAL(KIND=8) :: D_u

   REAL(KIND=8) :: dp, dv
   !------------------------------------------------------------------------

   IF (type ==  'S') THEN 

      dp = p_RH(v, v0, p0) - p0              
  
      dv = v - v0               
                        
      ! Il segno è quello del termine per l'onda di sinistra, negativo          
      D_u = wsign * SIGN(1.d0, dv) * (dp + dv * Dp_RH(v, v0,p0)) / (2 * SQRT(-dv * dp))

   ELSE

      CALL Init_IC(p0, v0)

      D_u = - wsign * SQRT_DQ(v, dv)  ! dv is a dummy argument. Not elegant, to fix

   ENDIF

END FUNCTION Du





FUNCTION Dp(type, v, v0, p0)   RESULT(D_p)
  
   !------------------------------------------------------------------------
   IMPLICIT NONE

   CHARACTER(LEN=1), INTENT(IN) :: type
   REAL(KIND=8),     INTENT(IN) :: v, v0, p0

   REAL(KIND=8) :: D_p
   !------------------------------------------------------------------------

   IF (type ==  'S') THEN     
      D_p = Dp_RH(v, v0, p0)
   ELSE
      D_p = DQ(v, v0, p0)
   ENDIF

END FUNCTION Dp





SUBROUTINE Init_Wave_Sigma (wave)
  
   !------------------------------------------------------------------------  
   IMPLICIT NONE

   TYPE(nonlin_wave), INTENT(INOUT) :: wave 

   REAL(KIND=8) :: Pi, nu

   INTEGER :: nr, i, k
   !------------------------------------------------------------------------  

!  wave % sigma  => NULL()  !  LQ 

   nr = SIZE(wave % p, 1)
 
!  IF (ASSOCIATED(wave % sigma)) DEALLOCATE(wave % sigma)  !  LQ 
   IF (ALLOCATED(wave % sigma)) DEALLOCATE(wave % sigma)  !  LQ 

   ALLOCATE (wave % sigma(nr,2))

   DO i = 1, SIZE(wave % type)

      DO k = 1, 2
 
         Pi = wave % p(i,k) / 27
         nu = wave % v(i,k) * 3
 
         wave % sigma(i,k) = sigma__Pinu(Pi, nu)

      ENDDO

   ENDDO

END SUBROUTINE Init_Wave_Sigma





SUBROUTINE Classical_NonClassical (wave_1, wave_3,  type_1, type_3)
  
   ! This subroutine is used to identify the classical or nonclassical
   ! type of the resulting waves for Riemann problem.   Precisely  the
   ! nonclassical waves are printed in capital letters.
  
   !-------------------------------------------------------------------------
   IMPLICIT NONE
  
   TYPE(nonlin_wave),       INTENT(IN)  :: wave_1, wave_3 
   CHARACTER, DIMENSION(:), INTENT(OUT) :: type_1, type_3 
  
   INTEGER :: i
  
   LOGICAL :: class_shock, class_fan
   !-------------------------------------------------------------------------

   type_1 = wave_1 % type

   DO i = 1, SIZE(wave_1 % type)
  
       IF (wave_1 % type(i) == 'S') THEN
      
          class_shock  =  wave_1 % v(i,1)  >  wave_1 % v(i,2)
      
          IF (class_shock)   type_1(i) = 's'        
        
       ELSE
     
          class_fan  =  wave_1 % v(i,1)  <  wave_1 % v(i,2)
         
          IF (class_fan)   type_1(i) = 'f'        
        
       ENDIF  
  
   ENDDO
 

   type_3 = wave_3 % type

   DO i = 1, SIZE(wave_3 % type)
  
       IF (wave_3 % type(i) == 'S') THEN
      
          class_shock  =  wave_3 % v(i,2)  >  wave_3 % v(i,1)
        
          IF (class_shock)   type_3(i) = 's'        
        
       ELSE
     
          class_fan  =  wave_3 % v(i,2)  <  wave_3 % v(i,1)
        
          IF (class_fan)   type_3(i) = 'f'        
        
       ENDIF    
  
   ENDDO

  
END SUBROUTINE Classical_NonClassical


END PROGRAM EulVdW_NCRsolver
