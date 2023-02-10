
! =========================================================================
!   Description: Module for the Computation  of  wave  form  from convex 
!                envelope construction.  The  output  is  variable  wave
!                that contains the type of the resulting wave connecting
!                L_state  with  R_state, the speeds of the components of
!                the composite wave and the values of the states at each
!                side of the component waves.
!                Entropy solution  is  verified  by  mean of the Oleinik
!                condition, that is:
!
!                 v_L > v_R  --->  LOWER Convex Envelope 
!                 v_L < v_R  --->  UPPER Convex Envelope
!
!                Possible wave types considered are the following.
!        
!                Pure wave types: 
!
!                 S  (Shock)
!                 F  (Fan)
!
!
!                Double Mixed wave types:
!
!                 SF  (Shock - Fan) 
!                 FS  (Fan - Shock)
!
!
!                Triple Mixed wave types:
!
!                 SFS  (Shock - Fan - Shock)
!                 FSF  (Fan - Shock - Fan)
!
!
!   	 Author: Marco Fossati
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!
!          Year: 2005, October
! ==========================================================================

MODULE nonlinear_wave

   USE eigenvalues
   USE problem_params
   USE plot_utils
   USE numerics
   USE rhLocus_intCurve
   USE var_types
   USE vdw_gas

CONTAINS



FUNCTION Wave_Structure(sign, v, w_i, vI_i, vE_i)   RESULT(wave)

   ! Computation of the wave structure connecting state  v  with reference 
   ! state w_i = (v_i, u_i), by means of convex  envelope  considerations.
   ! The structure is saved in variable wave, and accounts for:
   ! 
   !   wave type (S,  SF,  SFS,  FSF,  FS,  F)
   !
   !   wave v: the values of specific volume at the sides of each
   !                  component of the wave
   !
   !   wave u: the values of velocity at the sides of each compo-
   !                  nent of the wave
   !
   !   wave speeds: speeds of the waves. Namely the shock speed  or  the
   !                speeds of the two characteristics delimiting the fan

   !-----------------------------------------------------------------------	     
   IMPLICIT NONE 								     
  
   INTEGER,      	             INTENT(IN) :: sign  								    
   REAL(KIND=8), 	             INTENT(IN) :: v  
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: w_i
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: vI_i, vE_i
 
   TYPE (nonlin_wave) :: wave				     
   
   TYPE (key_points) :: v_tan

   REAL(KIND=8) :: v_i, u_i, p_i,     &
                   sigma_i, nu, Pi,   &
                   s, nu_int, Pi_int, &
                   SL_int, v_tg_i,    &
                   vT_i, pT_i, uT_i,  &
                   vT_v, pT_v, uT_v,  &
                   vE_1, pE_1, uE_1,  &
                   vE_2, pE_2, uE_2,  &
                   dnu
                  
   INTEGER :: j, N_int_pts

   LOGICAL, DIMENSION(5) :: rel_pos

   LOGICAL :: above, follow_P_si, &
              A1, A2, A3, A4, B1, &
              exotic_waves
  !-----------------------------------------------------------------------

   WRITE(12,*) ''
   WRITE(12,*) '   NONLINEAR_WAVE'
  

   wave % sign = sign

   v_i = w_i(1)
   u_i = w_i(2)
   p_i = w_i(3)

   ! Scaling for VdW eos: Reference State (_i)
   nu_i = v_i * 3
   Pi_i = p_i / 27

   sigma_i = sigma__Pinu(Pi_i, nu_i)

   exotic_waves = sigma_i < sigma_star(delta)


   ! Inizialization of sigma_i for function Pi_vdW in Module vdw_gas. 
   CALL Init_Pi_vdW__sigmanu(sigma_i)

   ! Scaling for VdW eos: Variable Intermediate State
   nu = v * 3
   Pi = Pi_vdW__sigmanu(nu)


   WRITE(12,*) '     local L', w_i(1)           
   WRITE(12,*) '     local R', v                
   WRITE(12,*) '     sigma_i', sigma_i          
   WRITE(12,*) '     exotic?', exotic_waves     


   ! Comparing pressure values on curve and Straight Line
   ! at in five points included in [v_i, v]. Scaled variables 
   ! are assumed and NOT reduced variables.

   dnu = ABS(nu_i - nu) / 6

   DO j = 1, 5

      nu_int = MIN(nu_i, nu) + j*dnu
  
      Pi_int = Pi_vdW__sigmanu(nu_int)

      SL_int = Pi_i  +  (Pi_i - Pi) * (nu_int - nu_i) / (nu_i - nu)  

      rel_pos(j) = Pi_int >= SL_int

   ENDDO

   above = ALL(rel_pos)

  
   ! Lower Envelope Construction. Compressive Waves Branch
   IF (v < v_i) THEN

      WRITE(12,*) '     LOWER ENVELOPE'

      IF (DDPi_vdW__sigmanu(3*v_i) > 0) THEN
 
         follow_P_si = .FALSE.

      ELSE

         v_tg_i = v
             
         IF (exotic_waves) THEN

            ! Computing possible tangent point from v_i.
            v_tan = Tangency_Point(p_i, v_i,  v_i, v)

            IF (v_tan % N_pts == 1) v_tg_i = v_tan % x(1)

         ENDIF 

         s = (Pi_vdW__sigmanu(3*v_tg_i) - Pi_vdW__sigmanu(3*v_i)) / (3*v_tg_i - 3*v_i)

         follow_P_si = DPi_vdW__sigmanu(3*v_i) < s

      ENDIF

      WRITE(12,*) '     follow_Pi', follow_P_si



      IF (follow_P_si) THEN

         IF (exotic_waves) THEN

            A1 = MAXVAL(vI_i) <= v
            A2 = MINVAL(vI_i) >= v_i
            A3 = MINVAL(vI_i) <= v
            A4 = MAXVAL(vI_i) >= v_i

         ELSE

            A1 = .TRUE. 
            A2 = .TRUE. 
            A3 = .TRUE. 
            A4 = .TRUE. 

         ENDIF

         B1 = A3 .AND. A4

         IF (A1  .OR.  A2  .OR.  B1) THEN   ! No Infection point is included in (v, v_i)

            WRITE(12,*) ''
            WRITE(12,*) '     TEMPORARY Wave type ---> F'

            ALLOCATE (wave % type(1))
            ALLOCATE (wave % v(1,2), wave % u(1,2), wave % p(1,2))
            ALLOCATE (wave % speeds(1,2))

            wave % type = 'F'

            wave % v(1,1) = v_i;     wave % v(1,2) = v
            wave % p(1,1) = p_i;     wave % p(1,2) = Q(v, v_i, p_i)
            wave % u(1,1) = u_i;     wave % u(1,2) = u_i  +  Integral_Curve(sign, v_i, v, v_i, p_i)

            ! Waves' Speeds
            wave % speeds(1,1) = lambda(sign, v_i, p_i, u_i)
            wave % speeds(1,2) = lambda(sign, v,  wave % p(1,2),  wave % u(1,2))

         ELSE

            IF (exotic_waves) THEN

               A1 = MINVAL(vE_i) >= v_i
               A2 = MAXVAL(vE_i) <= v

               B1 = A1 .AND. A2

               IF (.NOT.B1) THEN   ! No Absolute Envelope point is included in (v, v_i)

                  WRITE(12,*) ''
                  WRITE(12,*) '     TEMPORARY Wave type ---> FS'

                  ! Compute Tangent point from guess state v
                  v_tan = Tangency_Point(27*Pi, v, v_i, v)

                  vT_v = v_tan % x(1);     pT_v = Q(vT_v, v_i, p_i)


                  ALLOCATE (wave % type(2))
                  ALLOCATE (wave % v(2,2), wave % u(2,2), wave % p(2,2))
                  ALLOCATE (wave % speeds(2,2))

                  wave % type = (/'F','S'/)

                  ! First Wave Component
                  uT_v = u_i  +  Integral_Curve(sign, v_i, vT_v, v_i, p_i)

                  wave % v(1,1) = v_i;     wave % v(1,2) = vT_v
                  wave % p(1,1) = p_i;     wave % p(1,2) = pT_v
                  wave % u(1,1) = u_i;     wave % u(1,2) = uT_v

                  ! Second Wave Component
                  wave % v(2,1) = vT_v;    wave % v(2,2) = v
                  wave % p(2,1) = pT_v;    wave % p(2,2) = p_RH(v,  vT_v, pT_v)
                  wave % u(2,1) = uT_v;    wave % u(2,2) = uT_v  +  u_RH(sign, v, vT_v, pT_v)

                  ! Waves' Speeds
                  wave % speeds(1,1) = lambda(sign, v_i,  p_i,  u_i)
                  wave % speeds(1,2) = lambda(sign, vT_v, pT_v, uT_v)

                  wave % speeds(2,:) = Shock_Speed(wave % v(2,:), wave % u(2,:))

               ENDIF

            ENDIF   ! End Exotic Waves

         ENDIF 


      ELSE  ! Not follow_P_si


         N_int_pts = Internal_Intersection_Points(Pi_i, nu_i, Pi, nu, v_i, v)


         IF (N_int_pts <= 0  .AND.  .NOT. above) THEN   ! No intersection point included and isentrope below

            WRITE(12,*) ''
            WRITE(12,*) '     TEMPORARY Wave type ---> S'

            ALLOCATE (wave % type(1))
            ALLOCATE (wave % v(1,2), wave % u(1,2), wave % p(1,2))
            ALLOCATE (wave % speeds(1,2))
 
            wave % type = 'S'

            wave % v(1,1) = v_i;     wave % v(1,2) = v
            wave % p(1,1) = p_i;     wave % p(1,2) = p_RH(v,  v_i, p_i)
            wave % u(1,1) = u_i;     wave % u(1,2) = u_i  +  u_RH(sign, v, v_i, p_i)

            ! Waves' Speeds
            wave % speeds(1,:) = Shock_Speed(wave % v(1,:), wave % u(1,:))

         ELSE

            IF (exotic_waves) THEN

               ! Compute Tangent point from I state
               v_tan = Tangency_Point(p_i, v_i, v_i, v)

               IF (v_tan % N_pts == 0) THEN

                  A1 = .FALSE.
                  A2 = .FALSE.
                  B1 = .FALSE.

               ELSE

                  vT_i = v_tan % x(1)

                  A1 = MAXVAL(vI_i) <= v
                  A2 = MINVAL(vI_i) >= vT_i
                  A3 = MINVAL(vI_i) <= v
                  A4 = MAXVAL(vI_i) >= vT_i

                  B1 = A3 .AND. A4

               ENDIF

               IF (A1  .OR.  A2  .OR.  B1) THEN   ! No Infection point in (v, vT_i)

                  WRITE(12,*) ''
                  WRITE(12,*) '     TEMPORARY Wave type ---> SF'

                  pT_i = p_RH(vT_i,  v_i, p_i)

                  ALLOCATE (wave % type(2))
                  ALLOCATE (wave % v(2,2), wave % u(2,2), wave % p(2,2))
                  ALLOCATE (wave % speeds(2,2))
 
                  wave % type = (/'S', 'F'/)

                  ! First Wave Component
                  uT_i = u_i  +  u_RH(sign, vT_i, v_i, p_i)

                  wave % v(1,1) = v_i;     wave % v(1,2) = vT_i
                  wave % p(1,1) = p_i;     wave % p(1,2) = pT_i
                  wave % u(1,1) = u_i;     wave % u(1,2) = uT_i


                  ! Second Wave Component
                  wave % v(2,1) = vT_i;    wave % v(2,2) = v
                  wave % p(2,1) = pT_i;    wave % p(2,2) = Q(v, vT_i, pT_i)
                  wave % u(2,1) = uT_i;    wave % u(2,2) = uT_i  +  Integral_Curve(sign, vT_i, v, vT_i, pT_i)


                  ! Waves' Speeds
                  wave % speeds(1,:) = Shock_Speed(wave % v(1,:), wave % u(1,:))
 
                  wave % speeds(2,1) = lambda(sign, vT_i, pT_i, uT_i)
                  wave % speeds(2,2) = lambda(sign, v,  wave % p(2,2),  wave % u(2,2))

               ELSE

                  WRITE(12,*) ''
                  WRITE(12,*) '     TEMPORARY Wave type ---> SFS'

                  ! Compute Tangent point from I state
                  v_tan = Tangency_Point(p_i, v_i, v_i, v)

                  vT_i = v_tan % x(1);   pT_i = p_RH(vT_i,  v_i, p_i)

                  ! Compute Tangent point from guess state v
                  v_tan = Tangency_Point(27*Pi, v, v_i, v)

                  vT_v = v_tan % x(1);   pT_v = Q(vT_v, vT_i, pT_i)


                  ALLOCATE (wave % type(3))
                  ALLOCATE (wave % v(3,2), wave % u(3,2), wave % p(3,2))
                  ALLOCATE (wave % speeds(3,2))

                  wave % type = (/'S', 'F', 'S'/)

                  ! First Wave Component
                  uT_i = u_i  +  u_RH(sign, vT_i, v_i, p_i)

                  wave % v(1,1) = v_i;     wave % v(1,2) = vT_i
                  wave % p(1,1) = p_i;     wave % p(1,2) = pT_i
                  wave % u(1,1) = u_i;     wave % u(1,2) = uT_i


                  ! Second Wave Component
                  uT_v = uT_i  +  Integral_Curve(sign, vT_i, vT_v, vT_i, pT_i)

                  wave % v(2,1) = vT_i;    wave % v(2,2) = vT_v
                  wave % p(2,1) = pT_i;    wave % p(2,2) = pT_v
                  wave % u(2,1) = uT_i;    wave % u(2,2) = uT_v
 

                  ! Third Wave Component
                  wave % v(3,1) = vT_v;    wave % v(3,2) = v
                  wave % p(3,1) = pT_v;    wave % p(3,2) = p_RH(v,  vT_v, pT_v)
                  wave % u(3,1) = uT_v;    wave % u(3,2) = uT_v  +  u_RH(sign, v, vT_v, pT_v)



                  ! Waves' Speeds
                  wave % speeds(1,:) = Shock_Speed(wave % v(1,:), wave % u(1,:))

                  wave % speeds(2,1) = lambda(sign, vT_i, pT_i, uT_i)
                  wave % speeds(2,2) = lambda(sign, vT_v, pT_v, uT_v)

                  wave % speeds(3,:) = Shock_Speed(wave % v(3,:), wave % u(3,:))


               ENDIF

            ENDIF  ! End Exotic Waves

         ENDIF

      ENDIF

   ENDIF





  ! Upper Envelope Construction. Rarefactive Waves Branch
  IF (v_i < v) THEN

      WRITE(12,*) '     UPPER ENVELOPE'

      IF (DDPi_vdW__sigmanu(3*v_i) < 0) THEN

         follow_P_si = .FALSE.
 
      ELSE
 
         v_tg_i = v
 
         IF (exotic_waves) THEN
 
            ! Computing possible tangent point from v_i.
            v_tan = Tangency_Point(p_i, v_i,  v_i, v)

            IF (v_tan % N_pts == 1) v_tg_i = v_tan % x(1)
 
         ENDIF
 
         s = (Pi_vdW__sigmanu(3*v_tg_i) - Pi_vdW__sigmanu(3*v_i)) / (3*v_tg_i - 3*v_i)
 
         follow_P_si = DPi_vdW__sigmanu(3*v_i) < s
 
      ENDIF

      WRITE(12,*) '     follow_Pi', follow_P_si



      IF (follow_P_si) THEN

         IF (exotic_waves) THEN

            A1 = MAXVAL(vI_i) <= v_i
            A2 = MINVAL(vI_i) >= v
            A3 = MINVAL(vI_i) <= v_i
            A4 = MAXVAL(vI_i) >= v

         ELSE

            A1 = .TRUE. 
            A2 = .TRUE. 
            A3 = .TRUE. 
            A4 = .TRUE. 

         ENDIF

         B1 = A3 .AND. A4

         IF (A1  .OR.  A2  .OR.  B1) THEN   ! No Infection point is included in (v, v_i)

            WRITE(12,*) ''
            WRITE(12,*) '     TEMPORARY Wave type ---> F wave'

            ALLOCATE (wave % type(1))
            ALLOCATE (wave % v(1,2), wave % u(1,2), wave % p(1,2))
            ALLOCATE (wave % speeds(1,2))

            wave % type = 'F'

            wave % v(1,1) = v_i;     wave % v(1,2) = v
            wave % p(1,1) = p_i;     wave % p(1,2) = Q(v, v_i, p_i)
            wave % u(1,1) = u_i;     wave % u(1,2) = u_i  +  Integral_Curve(sign, v_i, v, v_i, p_i)

            ! Waves' Speeds
            wave % speeds(1,1) = lambda(sign, v_i, p_i, u_i)
            wave % speeds(1,2) = lambda(sign, v,  wave % p(1,2),  wave % u(1,2))

         ELSE

            IF (exotic_waves) THEN

               A1 = MINVAL(vE_i) >= v_i
               A2 = MAXVAL(vE_i) <= v

               B1 = A1 .AND. A2


               IF (.NOT.  B1) THEN   ! No both Absolute Envelope points are included in (v_i, v)

                  WRITE(12,*) ''
                  WRITE(12,*) '     TEMPORARY Wave type ---> FS'

                  ! Compute Tangent point from guess state v
                  v_tan = Tangency_Point(27*Pi, v, v_i, v)

                  vT_v = v_tan % x(1);     pT_v = Q(vT_v, v_i, p_i)


                  ALLOCATE (wave % type(2))
                  ALLOCATE (wave % v(2,2), wave % u(2,2), wave % p(2,2))
                  ALLOCATE (wave % speeds(2,2))

                  wave % type = (/'F','S'/)

                  ! First Wave Component
                  uT_v = u_i  +  Integral_Curve(sign, v_i, vT_v, v_i, p_i)

                  wave % v(1,1) = v_i;     wave % v(1,2) = vT_v
                  wave % p(1,1) = p_i;     wave % p(1,2) = pT_v
                  wave % u(1,1) = u_i;     wave % u(1,2) = uT_v

                  ! Second Wave Component
                  wave % v(2,1) = vT_v;    wave % v(2,2) = v
                  wave % p(2,1) = pT_v;    wave % p(2,2) = p_RH(v,  vT_v, pT_v)
                  wave % u(2,1) = uT_v;    wave % u(2,2) = uT_v  +  u_RH(sign, v, vT_v, pT_v)

                  ! Waves' Speeds
                  wave % speeds(1,1) = lambda(sign, v_i,  p_i,  u_i)
                  wave % speeds(1,2) = lambda(sign, vT_v, pT_v, uT_v)

                  wave % speeds(2,:) = Shock_Speed(wave % v(2,:), wave % u(2,:))

               ELSE

                  WRITE(12,*) ''
                  WRITE(12,*) '     TEMPORARY Wave type ---> FSF'

                  vE_1 = MINVAL(vE_i);     pE_1 =    Q(vE_1,  v_i,  p_i)
                  vE_2 = MAXVAL(vE_i);     pE_2 = p_RH(vE_2,  vE_1, pE_1)


                  ALLOCATE (wave % type(3))
                  ALLOCATE (wave % v(3,2), wave % u(3,2), wave % p(3,2))
                  ALLOCATE (wave % speeds(3,2))
 
                  wave % type = (/'F','S','F'/)
 
                  ! First Wave Component
                  uE_1 = u_i  +  Integral_Curve(sign, v_i, vE_1, v_i, p_i)
 
                  wave % v(1,1) = v_i;     wave % v(1,2) = vE_1
                  wave % p(1,1) = p_i;     wave % p(1,2) = pE_1
                  wave % u(1,1) = u_i;     wave % u(1,2) = uE_1
 

                  ! Second Wave Component
                  uE_2 = uE_1  +  u_RH(sign, vE_2, vE_1, pE_1)
 
                  wave % v(2,1) = vE_1;    wave % v(2,2) = vE_2
                  wave % p(2,1) = pE_1;    wave % p(2,2) = pE_2
                  wave % u(2,1) = uE_1;    wave % u(2,2) = uE_2
 

                  ! Third Wave Component
                  wave % v(3,1) = vE_2;    wave % v(3,2) = v
                  wave % p(3,1) = pE_2;    wave % p(3,2) = Q(v, vE_2, pE_2)
                  wave % u(3,1) = uE_2;    wave % u(3,2) = uE_2  +  Integral_Curve(sign, vE_2, v, vE_2, pE_2)


                  ! Waves' Speeds
                  wave % speeds(1,1) = lambda(sign, v_i, p_i, u_i)
                  wave % speeds(1,2) = lambda(sign, vE_1, pE_1, uE_1)
 
                  wave % speeds(2,:) = Shock_Speed(wave % v(2,:), wave % u(2,:))
 
                  wave % speeds(3,1) = lambda(sign, vE_2, pE_2, uE_2)
                  wave % speeds(3,2) = lambda(sign, v,  wave % p(3,2),  wave % u(3,2))

               ENDIF

            ENDIF

         ENDIF


      ELSE  ! Not follow_P_si

         N_int_pts = Internal_Intersection_Points(Pi_i, nu_i, Pi, nu, v_i, v)

         IF (N_int_pts <= 0  .AND.  above) THEN   ! No intersection point included and isentrope below

            WRITE(12,*) ''
            WRITE(12,*) '     TEMPORARY Wave type ---> S'

            ALLOCATE (wave % type(1))
            ALLOCATE (wave % v(1,2), wave % u(1,2), wave % p(1,2))
            ALLOCATE (wave % speeds(1,2))
 
            wave % type = 'S'

            wave % v(1,1) = v_i;     wave % v(1,2) = v
            wave % p(1,1) = p_i;     wave % p(1,2) = p_RH(v,  v_i, p_i)
            wave % u(1,1) = u_i;     wave % u(1,2) = u_i  +  u_RH(sign, v, v_i, p_i)

            ! Waves' Speeds
            wave % speeds(1,:) = Shock_Speed(wave % v(1,:), wave % u(1,:))

         ELSE

            IF (exotic_waves) THEN

               ! Compute Tangent point from I state
               v_tan = Tangency_Point(p_i, v_i, v_i, v)

               IF (v_tan % N_pts == 0) THEN

                  A1 = .FALSE.
                  A2 = .FALSE.
                  B1 = .FALSE.

               ELSE

                  vT_i = v_tan % x(1)

                  A1 = MAXVAL(vI_i) <= vT_i
                  A2 = MINVAL(vI_i) >= v
                  A3 = MINVAL(vI_i) <= vT_i
                  A4 = MAXVAL(vI_i) >= v

                  B1 = A3 .AND. A4

               ENDIF


               IF (A1  .OR.  A2  .OR.  B1) THEN   ! No Inflection point in (v, vT_i)

                  WRITE(12,*) ''
                  WRITE(12,*) '     TEMPORARY Wave type ---> SF'

                  pT_i = p_RH(vT_i,  v_i, p_i)

                  ALLOCATE (wave % type(2))
                  ALLOCATE (wave % v(2,2), wave % u(2,2), wave % p(2,2))
                  ALLOCATE (wave % speeds(2,2))
 
                  wave % type = (/'S', 'F'/)

                  ! First Wave Component
                  uT_i = u_i  +  u_RH(sign, vT_i, v_i, p_i)

                  wave % v(1,1) = v_i;     wave % v(1,2) = vT_i
                  wave % p(1,1) = p_i;     wave % p(1,2) = pT_i
                  wave % u(1,1) = u_i;     wave % u(1,2) = uT_i


                  ! Second Wave Component
                  wave % v(2,1) = vT_i;    wave % v(2,2) = v
                  wave % p(2,1) = pT_i;    wave % p(2,2) = Q(v, vT_i, pT_i)
                  wave % u(2,1) = uT_i;    wave % u(2,2) = uT_i  +  Integral_Curve(sign, vT_i, v, vT_i, pT_i)


                  ! Waves' Speeds
                  wave % speeds(1,:) = Shock_Speed(wave % v(1,:), wave % u(1,:))
 
                  wave % speeds(2,1) = lambda(sign, vT_i, pT_i, uT_i)
                  wave % speeds(2,2) = lambda(sign, v,  wave % p(2,2),  wave % u(2,2))

               ENDIF

            ENDIF

         ENDIF

      ENDIF

   ENDIF


END FUNCTION Wave_Structure
 




!***************************************************************************

FUNCTION Inflection_Points(sigma)   RESULT(vI)
   
   USE numerics
 
   IMPLICIT NONE 
    
   !-----------------------------------------------------------------------
 
   REAL(KIND=8), INTENT(IN) :: sigma

   REAL(KIND=8), DIMENSION(2) :: vI

   TYPE (key_points) :: v_inf

   REAL(KIND=8), DIMENSION(2) :: nuI
 
   REAL(KIND=8), PARAMETER :: t = 1.03d0
   
   INTEGER :: iter
   !-----------------------------------------------------------------------

   !INTEGER :: j
   ! 
   !REAL(KIND=8), DIMENSION(1000) :: vv, Newt_fun 

   WRITE(12,*) ''
   WRITE(12,*) '   INFLECTION_POINTS'


   !!-------------------------------------------------------------
   !DO j = 1, 1000
   !
   ! vv(j) = 0.34 + j*0.0566
   ! Newt_fun(j) = Infl(vv(j)*3) - (delta+1)*(delta+2)*sigma / 6
   !
   !ENDDO
   !
   !CALL Plot_Profile(vv, Newt_fun, vv(1), 'NewtonInflections   ')
   !!-------------------------------------------------------------



   ! Initial Guesses for Newton are assumed equal to inflection points
   ! for the isothermal case (delta=0) with t=1.03
   
   v_inf = IsoT_Inflection_Points(t)
   
   iter = Max_Iter

   nuI(1) = Newton(Infl, D_Infl, 3 * v_inf % x(1), eps, iter, &
                   (delta+1)*(delta+2)*sigma / 6 )

   iter = Max_Iter

   nuI(2) = Newton(Infl, D_Infl, 3 * v_inf % x(2), eps, iter, &
                   (delta+1)*(delta+2)*sigma / 6 )

	     

   !IF (ABS(nuI(2) - nuI(1)) < 1.d-6) THEN
   !
   !  WRITE(12,*) '     Assuming perturbed inflection points as'
   !  WRITE(12,*) '     initial guess for Absolute Envelope...'
   !
   !  ! Using perturbed Inflection Points as guesses for Newton Method
   !  v_inf % x = (v_inf % x + (/-1.d-2, 1.d-2/))
   !
   !  iter = Max_Iter
   !
   !  nuI(1) = Newton(Infl, D_Infl, 3 * v_inf % x(1), eps, iter, &
   !                  (delta+1)*(delta+2)*sigma / 6 )
   !
   !  iter = Max_Iter
   !
   !  nuI(2) = Newton(Infl, D_Infl, 3 * v_inf % x(2), eps, iter, &
   !                  (delta+1)*(delta+2)*sigma / 6 )
   !
   !ENDIF


   ! Check if Newton finds two distinct roots
   IF (ABS(nuI(2) - nuI(1)) < 1.d-6) THEN

      PRINT*, ''                   
      PRINT*, 'Error. INFLECTION_POINTS:'
      PRINT*, 'Newton Roots:', nuI/3     
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')       
      WRITE(idf,*) w_L          
      WRITE(idf,*) w_R          
      WRITE (idf,*) 'Wrong Inflection Points'   
      CLOSE(idf)           
      STOP             
                   
   ENDIF



   IF (iter < Max_Iter) THEN

      vI = nuI / 3

      WRITE(12,*) '     Inflection points', vI
 
   ELSE

      WRITE(12,*) '     Maximum number of iterations', Max_Iter, 'reached.'
      WRITE(12,*) '     Newton Roots', nuI/3
      WRITE(12,*) '     Continues...'
     
   ENDIF


CONTAINS

FUNCTION IsoT_Inflection_Points(t)  RESULT(v_inf)

   USE numerics                         

   !-----------------------------------------------------------------------  
   IMPLICIT NONE                                                             
 
   REAL(KIND=8), INTENT(IN) :: t                                             
 
   TYPE(key_points) :: v_inf                                                 
 
   REAL(KIND=8) r_guess, I1, I2 
   
   INTEGER :: iter                                             
   !-----------------------------------------------------------------------  
 
   iter = Max_Iter
 
   r_guess = 0.1                                                             
   I1 = Newton(h, Dh, r_guess, eps, iter, 8*t)                          


   iter = Max_Iter
   
   r_guess = 1.5                                                             
   I2 = Newton(h, Dh, r_guess, eps, iter, 8*t)                          
 
   v_inf % N_pts = 2                                                         
 
   ALLOCATE (v_inf % x(2))                                                   
 
   v_inf % x(1) = 1 / I2                                                     
   v_inf % x(2) = 1 / I1                                                     
 
END FUNCTION IsoT_Inflection_Points                                      
     

END FUNCTION Inflection_Points





!***************************************************************************

FUNCTION Absolute_Envelope(sigma)   RESULT(vE)

   USE numerics

   !-----------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: sigma

   REAL(KIND=8), DIMENSION(2) :: vE

   REAL(KIND=8), DIMENSION(2) :: v_E, nu_E
 
   REAL(KIND=8), PARAMETER :: t = 1.03d0

   INTEGER :: iter
   !-----------------------------------------------------------------------

   WRITE(12,*) ''
   WRITE(12,*) '   ABSOLUTE_ENVELOPE '


   ! Initial Guesses for Newton are assumed equal to envelope points
   ! for the isothermal case (delta=0) with t=1.03

   v_E = IsoT_Absolute_Envelope(t)
   
   iter = Max_Iter

   CALL Init_Pi_vdW__sigmanu(sigma)

!!! nu_E = Newton_Sys2_Diag_Jac(  Pi_vdW__sigmanu, &
   nu_E = Newton_Sys2_Jac(  Pi_vdW__sigmanu, &
                           DPi_vdW__sigmanu, &
                          DDPi_vdW__sigmanu, &
                           3*v_E, eps, iter)


   IF (ABS(nu_E(2) - nu_E(1)) < 1.d-6) THEN

      WRITE(12,*) '     Assuming perturbed inflection points as'
      WRITE(12,*) '     initial guess for Absolute Envelope...'

      ! Using Inflection Points as guesses for Newton Method
      v_E = Inflection_Points(sigma)

      v_E = (v_E + (/-1.d-2, 1.d-2/))

      iter = Max_Iter

!!!   nu_E = Newton_Sys2_Diag_Jac(  Pi_vdW__sigmanu, &
      nu_E = Newton_Sys2_Jac(  Pi_vdW__sigmanu, &
                              DPi_vdW__sigmanu, &
                             DDPi_vdW__sigmanu, &
                              3*v_E, eps, iter)

   ENDIF


   ! Check if Newton finds two distinct roots
   IF (ABS(nu_E(2) - nu_E(1)) < 1.d-6) THEN

      PRINT*, ''
      PRINT*, 'Error. ABSOLUTE_ENVELOPE:'
      PRINT*, 'Newton Roots:', nu_E / 3
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')       
        WRITE(idf,*) w_L          
        WRITE(idf,*) w_R          
        WRITE (idf,*) 'Wrong Absolute Envelope'   
      CLOSE(idf)           
      STOP             
                   
   ENDIF



   IF (iter < Max_Iter) THEN

      vE = nu_E / 3

      WRITE(12,*) '     Absolute Envelope points', vE
 
   ELSE
   
      WRITE(12,*) '     Maximum number of iterations', Max_Iter, 'reached.'
      WRITE(12,*) '     Newton Roots', nu_E/3
      WRITE(12,*) '     Continues...'
     
   ENDIF


CONTAINS

FUNCTION IsoT_Absolute_Envelope(t)  RESULT(v_E)                           

   ! Computes the two absolute envelope points. These are the two points of  
   ! tangency between the Pressure curve and a straight line in v-P plane.   

   !-----------------------------------------------------------------------   
   IMPLICIT NONE                                                             

   REAL(KIND=8), INTENT(IN) :: t                                             
   
   REAL(KIND=8), DIMENSION(2) :: v_E                                         
   !-----------------------------------------------------------------------   

   v_E(1) = 4 / (3*(1 + SQRT(9 - 16*SQRT(t*8/27))))         

   v_E(2) = 4 / (3*(1 - SQRT(9 - 16*SQRT(t*8/27))))         

END FUNCTION IsoT_Absolute_Envelope 
 
 
END FUNCTION Absolute_Envelope





!***************************************************************************

FUNCTION Tangency_Point(p0, v0, v_L, v_R)   RESULT(v_tan)
   !-----------------------------------------------------------------------   
   IMPLICIT NONE                    

   REAL(KIND=8), INTENT(IN) :: p0, v0, v_L, v_R

   TYPE(key_points) :: v_tan


   TYPE(key_points) :: a_r_v

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: nu_tan
   
   REAL(KIND=8) :: t0, v          

   INTEGER :: n, N_arv, iter
   !-----------------------------------------------------------------------

   !INTEGER :: j
   !
   !REAL(KIND=8), DIMENSION(1000) :: vv, Newt_fun


   WRITE(12,*) ''
   WRITE(12,*) '   TANGENCY_POINT'    
   WRITE(12,*) '     (v0, p0) = ', v0, p0

   t0 = t_vdW__pv(p0, v0)

   
   a_r_v = IsoT_Tangent_Points(v0, t0)     

   N_arv = a_r_v % N_pts
 
   IF (N_arv == 0) THEN
   
      ! If no tangent point  is  found for the isothermal  case (guess)
      ! then we are assuming that even in the general case no tangent point
      ! exists (???)

      WRITE(12,*) '     No tangent point in the isothermal case.'

      v_tan % N_pts = 0
    
      RETURN 

   ELSE
      
      v_tan % N_pts = 0   ! Initialize to zero 
 
      CALL Init_psi(3*v0, sigma_i)

      IF (ALLOCATED(nu_tan)) DEALLOCATE (nu_tan) ! LUIGI QUARTAPELLE

      ALLOCATE (nu_tan(N_arv))


      !!-------------------------------------------------------------
      !DO j = 1, 1000
      !
      !vv(j) = 0.34 + j*0.0566
      !Newt_fun(j) = psi(vv(j)*3)
      !
      !ENDDO
      !
      !CALL Plot_Profile(vv, Newt_fun, vv(1), 'NewtonTangent       ')
      !!-------------------------------------------------------------


      DO n = 1, N_arv

         iter = Max_Iter

         ! Custom_Newton works in adimensional variables (Pi, nu)
         nu_tan(n) = Custom_Newton(psi, Dpsi, 3*a_r_v % x(n), eps, iter)

      ENDDO


      ! Scegliamo solo quella che soddisfa le condizioni per essere
      ! vT_v0.

      WRITE(12,*) '     Existing tangents =', nu_tan / 3

      DO n = 1, N_arv
 
         ! Transforming to reduced specific volume
         v = nu_tan(n) / 3


         IF (ABS(v - v0) < 1.e-4   .OR.   v < 1.d0/3) CYCLE   

         IF ( (v < MAX(v_L, v_R) .AND. v > MIN(v_L, v_R))  .AND.  &
 
             ((v_L < v_R  .AND.  DDPi_vdW__sigmanu(3*v) > 0)   .OR.  &
              (v_R < v_L  .AND.  DDPi_vdW__sigmanu(3*v) < 0)) )    THEN
 
            v_tan % N_pts = 1
           
            ALLOCATE (v_tan % x(1))
            v_tan % x(1) = v

            WRITE(12,*) '     VALID tangent point abscissa', v_tan % x(1)

            RETURN

         ENDIF
 
      ENDDO

      WRITE(12,*) '     NO VALID tangent point continue...'

   ENDIF

END FUNCTION Tangency_Point



FUNCTION IsoT_Tangent_Points(v0, t0)  RESULT(a_r_v)

   ! Selects the tangency point through v0  included in [v_L, v_R]  such
   ! that the pressure second  derivative  in this point  is  consistent
   ! with convex envelope required.

   ! a_r_v = admissible real v solutions

   !-----------------------------------------------------------------------   
   IMPLICIT NONE                    

   REAL(KIND=8), INTENT(IN) :: v0, t0

   TYPE(key_points) :: a_r_v
	              

   REAL(KIND=8), DIMENSION(3) :: nu
   
   REAL(KIND=8) :: AA, BB, CC, &
                   tp_0, nu_0, den, nu_s            

   INTEGER :: N_sol                   

   LOGICAL :: single_real_root                  
   !-----------------------------------------------------------------------   

   tp_0 =  8 * t0 / 27
   nu_0 =  3 * v0

   den = nu_0 - 1 - tp_0 * nu_0**2

   IF (den == 0) THEN

      PRINT*, ''
      PRINT*, 'Error. TANGENCY_POINT:'
      PRINT*, 'Divide by zero.'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
        WRITE(idf,*) w_L
        WRITE(idf,*) w_R
        WRITE(idf,*) 'Tangency_Point is STOPPED. Divide by zero.'
      CLOSE(idf)
      STOP
 
   ENDIF


   AA = 2*(nu_0 - 1)**2
   BB = (nu_0 - 1) * (1 - 4*nu_0)
   CC = 2*(nu_0 - 1) * nu_0

   AA = AA / den
   BB = BB / den
   CC = CC / den

   nu = Cubic_Equation(AA, BB, CC,  single_real_root)

   
   ! IF (MAXVAL(nu) <= 1) ! E' equivalente dato che cubic equations mette
   ! zero in nu se esistono radici complesse.
   
   IF ((single_real_root .AND. nu(1) <= 1)  .OR. (MAXVAL(nu) <= 1)) THEN

       a_r_v % N_pts = 0
 
   ELSE  ! Esiste almeno una o due o tre' radici reali accettabili (> 1)
 
      ! Quante sono?
      N_sol = COUNT (nu > 1)

      a_r_v % N_pts = N_sol

      ! Ordinamento decrescente per considerare solo le radici > 1
      IF (nu(1) < nu(2)) THEN
         nu_s = nu(1);   nu(1) = nu(2);   nu(2) = nu_s
      ENDIF

      IF (nu(2) < nu(3)) THEN
         nu_s = nu(2);   nu(2) = nu(3);   nu(3) = nu_s
      ENDIF

      IF (nu(1) < nu(2)) THEN
         nu_s = nu(1);   nu(1) = nu(2);   nu(2) = nu_s
      ENDIF

      ALLOCATE (a_r_v % x(N_sol))
         
      a_r_v % x = nu(1:N_sol) / 3

   ENDIF

END FUNCTION IsoT_Tangent_Points





!***************************************************************************

FUNCTION Internal_Intersection_Points(Pi_i, nu_i, Pi, nu, v_L, v_R)   RESULT(N_int)
   !-----------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: Pi_i, nu_i, Pi, nu, v_L, v_R  !Bello schifo

   INTEGER :: N_int, iter

   REAL(KIND=8), DIMENSION(2) :: v_isoT_inter

   REAL(KIND=8) :: t_i, v1, v2
   !-----------------------------------------------------------------------

   ! INTEGER :: j
   !
   ! REAL(KIND=8), DIMENSION(1000) :: vv, Newt_fun


   WRITE(12,*) ''
   WRITE(12,*) '   INTERNAL_INTERSECTION_POINTS'


   ! Initial Guess for intersection points is obtained computing         
   ! solution of the isothermal case for the reduced temperature         
   ! associated with the reference state (_i)                            

   t_i = t_vdW__pv(27*Pi_i, nu_i/3)
                                                       !OUTPUT
   IF (IsoT_Internal_Intersection(nu/3, nu_i/3, t_i,   v_isoT_inter) == 0) THEN
                                                       !OUTPUT
      N_int = 0

      RETURN

   ELSE

      N_int = 2
   
      CALL Init_phi(Pi_i, nu_i, Pi, nu)                  


      !!-----------------------------------------------------------------
      !DO j = 1, 1000
      !
      !  vv(j) = 0.34 + j*0.166
      !  Newt_fun(j) = 27 * Phi(vv(j)*3)
      !
      !ENDDO
      !
      !CALL Plot_Profile(vv, Newt_fun, vv(1), 'NewtonIntersections ')
      !!-----------------------------------------------------------------



      iter = Max_Iter
                                                                                
      v1 = Custom_Newton(phi, Dphi, 3*v_isoT_inter(1), eps, iter) / 3

      iter = Max_Iter
      
      v2 = Custom_Newton(phi, Dphi, 3*v_isoT_inter(2), eps, iter) / 3


      IF (ABS(v1 - v_L) < 1.d-10) N_int = N_int - 1 
      IF (ABS(v1 - v_R) < 1.d-10) N_int = N_int - 1 
 
      IF (ABS(v2 - v_L) < 1.d-10) N_int = N_int - 1 
      IF (ABS(v2 - v_R) < 1.d-10) N_int = N_int - 1 
                       
      IF (v1 > MAX(v_L, v_R)  .OR.  v1 < MIN(v_L, v_R)) N_int = N_int - 1
      IF (v2 > MAX(v_L, v_R)  .OR.  v2 < MIN(v_L, v_R)) N_int = N_int - 1

   ENDIF

END FUNCTION Internal_Intersection_Points



FUNCTION IsoT_Internal_Intersection(v_L, v_R, t, vv)  RESULT(N_int)
       
   !-----------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v_L, v_R, t

   REAL(KIND=8), DIMENSION(2), INTENT(OUT) :: vv

   INTEGER :: N_int

   REAL(KIND=8) :: alpha, beta, discr, delta_p, v_1, v_2
   !-----------------------------------------------------------------------

   delta_p = p_vdW__vt(v_R, t) - p_vdW__vt(v_L, t)

   alpha = -1.d0/3.d0  +  v_L  +  v_R  & 
                       + ( v_R * p_vdW__vt(v_L, t) &
                         - v_L * p_vdW__vt(v_R, t) ) / delta_p

   beta  = - (1/(v_L * v_R)) * (v_R - v_L) / delta_p


   discr = alpha**2 - 4 * beta

   IF (discr >= 0) THEN

      ! Two distinct intersection points
      N_int = 2
 
      v_1 = (-alpha - SQRT(discr)) / 2
      v_2 = (-alpha + SQRT(discr)) / 2

      vv = (/v_1, v_2/)

   ELSE

      ! No intersection points
      N_int = 0
 
   ENDIF

END FUNCTION IsoT_Internal_Intersection





FUNCTION Shock_Speed(vv, uu)   RESULT(speed)
   !-----------------------------------------------------------------------
   
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: vv, uu

   REAL(KIND=8) :: speed
   !-----------------------------------------------------------------------
   
   IF (vv(1) == vv(2)) THEN

      PRINT*, ''
      PRINT*, 'Error. SHOCK_SPEED:'
      PRINT*, 'Equal specific volumes at the sides of the shock'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
        WRITE(idf,*) w_L
        WRITE(idf,*) w_R
        WRITE(idf,*) 'Equal specific volumes at the sides of the shock'
      CLOSE(idf)
      STOP

   ENDIF

   speed = uu(2) - vv(2) * (uu(2) - uu(1)) / (vv(2) - vv(1))


END FUNCTION Shock_Speed




   
SUBROUTINE Invert_Wave (wave, inverted_wave)

   ! Invert the order of component waves.  This is necessary  only  for  the  
   ! solution visualization and not for the identification of the intermediate  
   ! state.  The order of the subcopmonents of the waves is changed as follows:
   !
   !  from  1 2    to  2 1      for the double composite wave
   !  from  1 2 3  to  3 2 1    for the triple composite wave
   !
   !--------------------------------------------------------------------------
   IMPLICIT NONE
                                                                             
   TYPE (nonlin_wave), INTENT(IN) :: wave

   TYPE (nonlin_wave) :: inverted_wave

   INTEGER :: i, j, nw
   !--------------------------------------------------------------------------

   nw = SIZE(wave % type)

   DO i = 1, nw        
                                                                             
      inverted_wave % type(i) = wave % type(nw+1 - i)
                                                                             
      DO j = 1, 2                                         
                                         
         inverted_wave % v(i, j) = wave % v(nw+1 - i, 3 - j)
         inverted_wave % u(i, j) = wave % u(nw+1 - i, 3 - j)
         inverted_wave % p(i, j) = wave % p(nw+1 - i, 3 - j)

         inverted_wave % speeds(i, j) = wave % speeds(nw+1 - i, 3 - j)
         inverted_wave % sigma (i, j) = wave % sigma (nw+1 - i, 3 - j)

      ENDDO                                                                   
                                                                             
   ENDDO                                                                     
                                                                              
END SUBROUTINE Invert_Wave                                                

   

END MODULE nonlinear_wave
