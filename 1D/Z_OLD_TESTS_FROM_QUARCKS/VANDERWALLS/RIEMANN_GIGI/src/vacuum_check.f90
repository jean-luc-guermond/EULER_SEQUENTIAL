
! =========================================================================
!   Description: Vacuum Control Procedures
!
!   	 Author: Marco Fossati
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!
!          Year: 2006  February
! =========================================================================

MODULE vacuum_check

   USE eigenvalues
   USE nonlinear_wave
   USE problem_params
   USE rhLocus_intCurve
   USE var_types
   USE vdw_gas


CONTAINS


FUNCTION Check_Vacuum(sigma_L, vI_L, vE_L, &
                      sigma_R, vI_R, vE_R, &
                      wave_1, wave_3)           RESULT(vacuum) 
   !------------------------------------------------------------------------  
   IMPLICIT NONE
  
   REAL(KIND=8),               INTENT(IN)    :: sigma_L
   REAL(KIND=8), DIMENSION(2), INTENT(IN)    :: vI_L, vE_L
   REAL(KIND=8),               INTENT(IN)    :: sigma_R
   REAL(KIND=8), DIMENSION(2), INTENT(IN)    :: vI_R, vE_R
   TYPE(nonlin_wave),          INTENT(INOUT) :: wave_1, wave_3

   LOGICAL :: vacuum
  
   REAL(KIND=8) :: u_vac_L, u_vac_R, du_vac
    
   INTEGER :: i				    
   !------------------------------------------------------------------------  

   WRITE(12,*) ''
   WRITE(12,*) '   CHECK_VACUUM'

   u_vac_L = Vacuum_Velocity_Term(sigma_L, vI_L, vE_L, w_L,  -1, wave_1)
   WRITE(12,*) '     u_vac_l', u_vac_L

   u_vac_R = Vacuum_Velocity_Term(sigma_R, vI_R, vE_R, w_R,  +1, wave_3)
   WRITE(12,*) '     u_vac_R', u_vac_r


   du_vac = u_vac_L + u_vac_R


   IF (w_R(2) - w_L(2) >= du_vac) THEN
  
      vacuum = .TRUE. 

      ALLOCATE (type_1(SIZE(wave_1 % type)),  &
   	            type_3(SIZE(wave_3 % type)))

      type_1 = wave_1 % type
 
      DO i = 1, SIZE(type_1)
         IF (type_1(i) == 'F') type_1(i) = 'f'	 
      ENDDO

      type_3 = wave_3 % type
   	  
      DO i = 1, SIZE(type_3)
         IF (type_3(i) == 'F') type_3(i) = 'f'	 
      ENDDO

      PRINT*, '	  VACUUM FORMED...'
      PRINT*, '' 
      PRINT*, '	   - Resulting waves   ', type_1, '  ***  ', type_3
      PRINT*, ''	
      PRINT*, '	   - left  bound velocity = ',  w_L(2) + u_vac_L
      PRINT*, '	   - right bound velocity = ',  w_R(2) - u_vac_R
      PRINT*, ''
      PRINT*, '	  Newton Method for Intermediate States is SKIPPED.'
      PRINT*, ''

  ELSE
  
     vacuum = .FALSE.

  ENDIF
    
END FUNCTION Check_Vacuum 



FUNCTION Vacuum_Velocity_Term(sigma_ith, vI_i, vE_i, &
                              w_i, sign, wave_i)    RESULT(u_vac)
   !------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8),               INTENT(IN) :: sigma_ith
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: vI_i, vE_i
   REAL(KIND=8), DIMENSION(3), INTENT(IN) :: w_i
   INTEGER,                    INTENT(IN) :: sign

   TYPE(nonlin_wave),          INTENT(INOUT) :: wave_i

   REAL(KIND=8) :: u_vac

   TYPE(key_points) :: v_tan

   REAL(KIND=8) :: v_i, u_i, p_i, sigma_st, &
                   vT_i, uT_i, pT_i,        &
                   vE_1, uE_1, pE_1,        &
                   vE_2, uE_2, pE_2
   !------------------------------------------------------------------------

   WRITE(12,*) ''
   WRITE(12,*) '   VACUUM_VELOCITY_TERM'


   v_i = w_i(1) 
   u_i = w_i(2)
   p_i = w_i(3) 


   CALL Init_Pi_vdW__sigmanu(sigma_ith)
   CALL Init_psi(3*v_i, sigma_ith)
   CALL Init_phi(p_i/27, 3*v_i, 1.d-5, 1.d+5)

   WRITE(12,*) '     sigma', sigma_i 



   sigma_st = sigma_star(delta)
   

   IF (sigma_ith > sigma_st   .OR.  (sigma_ith < sigma_st  .AND.  v_i > MAXVAL(vI_i))) THEN

      WRITE(12,*) '     Wave type:   f'

      u_vac = ABS(Integral_Curve(1, v_i, 1.d+5, v_i, p_i))

      ! Computing Wave Components

      ALLOCATE (wave_i % type(1))
      ALLOCATE (wave_i % v(1,2), wave_i % u(1,2), wave_i % p(1,2))
      ALLOCATE (wave_i % speeds(1,2))

      wave_i % type = 'F'

      wave_i % v(1,1) = v_i;     wave_i % v(1,2) = 1.d+5
      wave_i % p(1,1) = p_i;     wave_i % p(1,2) = Q(1.d+5, v_i, p_i)
      wave_i % u(1,1) = u_i;     wave_i % u(1,2) = u_i  +  Integral_Curve(sign, v_i, 1.d+5, v_i, p_i)

      ! Waves' Speeds
      wave_i % speeds(1,1) = lambda(sign, v_i, p_i, u_i)
      wave_i % speeds(1,2) = lambda(sign, 1.d+5,  wave_i % p(1,2),  wave_i % u(1,2))

      RETURN

   ENDIF


   IF (sigma_ith < sigma_st  .AND.  (v_i > MINVAL(vE_i)  .AND.  v_i < MAXVAL(vI_i))) THEN

      WRITE(12,*) '     Wave_type:   Sf'      

      v_tan = Tangency_Point(p_i, v_i, v_i, 1.d+5)

      vT_i = v_tan % x(1)
      uT_i = u_RH(-1, vT_i, v_i, p_i)
      pT_i = p_RH(vT_i,  v_i, p_i)

      u_vac = uT_i  +  Integral_Curve(-1, vT_i, 1.d+5, vT_i, pT_i)



      ! Computing Wave Components

      ALLOCATE (wave_i % type(2))
      ALLOCATE (wave_i % v(2,2), wave_i % u(2,2), wave_i % p(2,2))
      ALLOCATE (wave_i % speeds(2,2))
 
      wave_i % type = (/'S','F'/)

      uT_i = u_i  +  u_RH(sign, vT_i, v_i, p_i)

      ! First Wave Component
      wave_i % v(1,1) = v_i;     wave_i % v(1,2) = vT_i
      wave_i % p(1,1) = p_i;     wave_i % p(1,2) = pT_i
      wave_i % u(1,1) = u_i;     wave_i % u(1,2) = uT_i

      ! Second Wave Component
      wave_i % v(2,1) = vT_i;    wave_i % v(2,2) = 1.d+5
      wave_i % p(2,1) = pT_i;    wave_i % p(2,2) = Q(1.d+5, vT_i, pT_i)
      wave_i % u(2,1) = uT_i;    wave_i % u(2,2) = uT_i  +  Integral_Curve(sign, vT_i, 1.d+5, vT_i, pT_i)

      ! Waves' Speeds 
      wave_i % speeds(1,:) = Shock_Speed(wave_i % v(1,:), wave_i % u(1,:))
 
      wave_i % speeds(2,1) = lambda(sign, vT_i, pT_i, uT_i)
      wave_i % speeds(2,2) = lambda(sign, 1.d+5,  wave_i % p(2,2),  wave_i % u(2,2))

      RETURN

   ENDIF


   IF (sigma_ith < sigma_st  .AND.  v_i < MINVAL(vI_i)) THEN

      WRITE(12,*) '     Wave_type:   fSf'

      vE_1 = MINVAL(vE_i)
      pE_1 = Q(vE_1,  v_i,  p_i)
      uE_1 = Integral_Curve(-1, v_i, vE_1, v_i, p_i)

      vE_2 = MAXVAL(vE_i)     
      pE_2 = p_RH(vE_2,  vE_1, pE_1)
      uE_2 = uE_1  +  u_RH(-1, vE_2, vE_1, pE_1)

      u_vac = uE_2  +  Integral_Curve(-1, vE_2, 1.d+5, vE_2, pE_2)



      ! Computing Wave Components

      ALLOCATE (wave_i % type(3))
      ALLOCATE (wave_i % v(3,2), wave_i % u(3,2), wave_i % p(3,2))
      ALLOCATE (wave_i % speeds(3,2))

      wave_i % type = (/'F','S','F'/)

      uE_1 =  u_i  +  Integral_Curve(sign, v_i, vE_1, v_i, p_i)
      uE_2 = uE_1  +  u_RH(sign, vE_2, vE_1, pE_1)

      ! First Wave Component
      wave_i % v(1,1) = v_i;     wave_i % v(1,2) = vE_1
      wave_i % p(1,1) = p_i;     wave_i % p(1,2) = pE_1
      wave_i % u(1,1) = u_i;     wave_i % u(1,2) = uE_1

      ! Second Wave Component
      wave_i % v(2,1) = vE_1;    wave_i % v(2,2) = vE_2
      wave_i % p(2,1) = pE_1;    wave_i % p(2,2) = pE_2
      wave_i % u(2,1) = uE_1;    wave_i % u(2,2) = uE_2

      ! Third Wave Component
      wave_i % v(3,1) = vE_2;    wave_i % v(3,2) = 1.d+5
      wave_i % p(3,1) = pE_2;    wave_i % p(3,2) = Q(1.d+5, vE_2, pE_2)
      wave_i % u(3,1) = uE_2;    wave_i % u(3,2) = uE_2  +  Integral_Curve(sign, vE_2, 1.d+5, vE_2, pE_2)

      ! Waves' Speeds 
      wave_i % speeds(1,1) = lambda(sign, v_i, p_i, u_i)
      wave_i % speeds(1,2) = lambda(sign, vE_1, pE_1, uE_1)

      wave_i % speeds(2,:) = Shock_Speed(wave_i % v(2,:), wave_i % u(2,:))

      wave_i % speeds(3,1) = lambda(sign, vE_2, pE_2, uE_2)
      wave_i % speeds(3,2) = lambda(sign, 1.d+5,  wave_i % p(3,2),  wave_i % u(3,2))

   ENDIF


END FUNCTION Vacuum_Velocity_Term



END MODULE vacuum_check
