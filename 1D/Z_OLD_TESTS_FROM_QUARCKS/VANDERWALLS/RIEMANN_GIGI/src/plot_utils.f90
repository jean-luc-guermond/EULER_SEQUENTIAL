
! ========================================================================
!   Description: Plot procedures 
!
!   	 Author: Marco Fossati
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!
!          Year:  2005, October
! =========================================================================
 
MODULE plot_utils
  
   USE problem_params
   USE rhLocus_intCurve
   USE var_types
   USE vdw_gas


CONTAINS


SUBROUTINE Plot_Isentropes (w_L, w_R, sigma_L, sigma_R)

   !-------------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: w_L, w_R
   REAL(KIND=8),               INTENT(IN) :: sigma_L, sigma_R

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v, v_gas, v_liq, &
                                              p, p_sat, p_curve
   
   REAL(KIND=8) :: dumm, sigma_I, nu, v_

   REAL(KIND=8), PARAMETER :: d_t = 0.01

   INTEGER :: i, nj, pind
   INTEGER, PARAMETER :: idf = 11, n = 1000
   !-------------------------------------------------------------------------

   ! Loading Saturation Curve data from file
   ALLOCATE (v_gas(n), v_liq(n), p_sat(n), p_curve(n))

   OPEN (UNIT = idf, FILE = 'curva_di_saturazione.dat', STATUS = 'unknown')

      READ(idf,*)
     
      DO i = 1, n
         READ(idf,*) dumm, p_sat(i), v_gas(i), v_liq(i)       
      ENDDO

   CLOSE (idf)
   

   ALLOCATE (p(n), v(n))

   OPEN (UNIT = idf, FILE = 'TEC__isentropes.plt', STATUS = 'unknown')
     
   WRITE(idf,*) 'TITLE = "Isentropes"'                   
   WRITE(idf,*) 'VARIABLES = "v", "p"'                  

   WRITE (idf,*)  ''                                        
   WRITE (idf,*)  'ZONE T="Left state", I=', 1
   WRITE (idf,*)  w_L(1), w_L(3)

   WRITE (idf,*)  ''                                        
   WRITE (idf,*)  'ZONE T="Right state", I=', 1
   WRITE (idf,*)  w_R(1), w_R(3)


   ! Writing Saturation Curve
   WRITE(idf,*) 'ZONE T ="saturazione", I =', 2*n-1 

   DO i = 1, n-1
      WRITE(idf,*) v_liq(n-i), p_sat(n-i)
   ENDDO 

   DO i = 1, n
      WRITE(idf,*) v_gas(i), p_sat(i)
   ENDDO 



   ! Writing Locus of G = 0. Polytropic Van der Waals                     
   ! gas 'c_v' is assumed constant.                                       

   DO i = 1, n                                                            

      v(i) = 1.d0 + (i-1)*0.002                                            

      p(i) = Gzero_Locus(v(i))                                      

      p_curve(i) = Linear_Interpolation(p_sat, v_gas, v(i))                

   ENDDO                                                                  

   nj = COUNT(p > p_curve)                                                

   IF (nj > 0) THEN

      WRITE(idf,*) ''                                                      
      WRITE(idf,*) 'ZONE T="d=', delta, '", I=', nj                        

      DO i = 1, n                                                          
         IF (p(i) > p_curve(i))  WRITE(idf,*) v(i), p(i)                    
      ENDDO
                                                              
   ENDIF 


   ! Writing Van der Waals Isentropes for left and right states


   ! Temporary modification to plot envelope points----------------------------
   DO i = 1, 55

      sigma_I = 0.31518 - (i-1)*0.0003036
      CALL Init_Pi_vdW__sigmanu(sigma_I)

      WRITE(idf,*) ''                                                      
      WRITE(idf,*) 'ZONE T="sigma L =', sigma_L, '", I=', n
   
      DO pind = 1, n       
         v(pind) = 0.35 + (pind-1)*0.00465                                                 
         WRITE(idf,*) v(pind), 27*Pi_vdW__sigmanu(3*v(pind))
      ENDDO    

   ENDDO 


   WRITE(idf,*) ''                                                      
   WRITE(idf,*) 'ZONE T="DSL", I=', n

   DO pind = 0, n-1

      v_ = 0.35 + pind*0.00465
      nu = 3*v_

      WRITE(idf,*) v_, 27*0.5 * ( (1-delta)**2              &
                              - (1-delta)*(1-3*delta)/nu    &
                              - 8*delta/nu**2               &
                              - 4/nu**3 )                   &
                              / ((1+delta)*(2+delta) * nu)
   ENDDO
   ! END of modification ------------------------------------------------------



   ! Inizialization of sigma_L.  sigma_L = wave_1 % sigma(1,1)
   CALL Init_Pi_vdW__sigmanu (sigma_L)

   DO i = 1, n                                                            

      v(i) = 0.35 + (i-1)*0.00465                                          

      p(i) = 27*Pi_vdW__sigmanu(3*v(i))                                 

   ENDDO                                                                  

   WRITE(idf,*) ''                                                        
   WRITE(idf,*) 'ZONE T="sigma L =', sigma_L, '", I=', n
   
   DO i = 1, n                                                            
      WRITE(idf,*) v(i), p(i)                                              
   ENDDO                                                                  

   
   ! Inizialization of sigma_R.  sigma_R = wave_3 % sigma(1,1)
   CALL Init_Pi_vdW__sigmanu(sigma_R)

   DO i = 1, n                                                            

      v(i) = 0.35 + (i-1)*0.00465                                          

      p(i) = 27*Pi_vdW__sigmanu(3*v(i))                                 

   ENDDO                                                                  

   WRITE(idf,*) ''                                                        
   WRITE(idf,*) 'ZONE T="sigma R =', sigma_R, '", I=', n
   
   DO i = 1, n                                                            
      WRITE(idf,*) v(i), p(i)                                              
   ENDDO                                                                  

   CLOSE (idf)


END SUBROUTINE Plot_Isentropes





SUBROUTINE Plot_vp_Waves (wave_1, wave_3)

   !-------------------------------------------------------------------------
   IMPLICIT NONE 

   TYPE(nonlin_wave), INTENT(IN) :: wave_1, wave_3

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v, v_gas, v_liq,   &
                                              p, p_sat, p_curve, &
                                              pRH, QIS

   REAL(KIND=8) :: v0, p0, dumm, dv, x_min, x_max, y_min, y_max

   INTEGER :: i, j, nj, nm
   INTEGER, PARAMETER :: idf = 11, n = 1000
   !-------------------------------------------------------------------------    

   ! Loading Saturation Curve data from file
   ALLOCATE (v_gas(n), v_liq(n), p_sat(n), p_curve(n))

   OPEN (UNIT = idf, FILE = 'curva_di_saturazione.dat', STATUS = 'unknown')

   READ(idf,*)
   
   DO i = 1, n
      READ(idf,*) dumm, p_sat(i), v_gas(i), v_liq(i)       
   ENDDO

   CLOSE (idf)


   ALLOCATE (p(n), v(n))

   OPEN (UNIT=idf, FILE='TEC__vp_waves.plt', STATUS='unknown')
   
   WRITE(idf,*) 'TITLE = "vp Waves"'                   
   WRITE(idf,*) 'VARIABLES = "v", "p"'                  

   WRITE (idf,*)  ''                                        
   WRITE (idf,*)  'ZONE T="Left state", I=', 1
   WRITE (idf,*)  wave_1 % v(1,1), wave_1 % p(1,1)

   WRITE (idf,*)  ''                                        
   WRITE (idf,*)  'ZONE T="Intermediate state Left", I=', 1
   WRITE (idf,*)  wave_1 % v(SIZE(wave_1 % type), 2), wave_1 % p(SIZE(wave_1 % type), 2)

   WRITE (idf,*)  ''                                        
   WRITE (idf,*)  'ZONE T="Intermediate state Right", I=', 1
   WRITE (idf,*)  wave_3 % v(SIZE(wave_3 % type), 2), wave_3 % p(SIZE(wave_3 % type), 2)
   
   WRITE (idf,*)  ''                                        
   WRITE (idf,*)  'ZONE T="Right state", I=', 1
   WRITE (idf,*)  wave_3 % v(1,1), wave_3 % p(1,1)


   ! Writing Saturation Curve
   WRITE(idf,*) 'ZONE T ="saturazione", I =', 2*n-1 

   DO i = 1, n-1
      WRITE(idf,*) v_liq(n-i), p_sat(n-i)
   ENDDO 

   DO i = 1, n
      WRITE(idf,*) v_gas(i), p_sat(i)
   ENDDO 



   ! Writing Locus of G = 0. Polytropic Van der Waals                     
   ! gas 'c_v' is assumed constant.                                       

   DO i = 1, n                                                            

      v(i) = 1.d0 + (i-1)*0.002                                            
 
      p(i) = Gzero_Locus(v(i))                                      

      p_curve(i) = Linear_Interpolation(p_sat, v_gas, v(i))                

   ENDDO                                                                  

   nj = COUNT(p > p_curve)                                                

   WRITE(idf,*) ''                                                        
   WRITE(idf,*) 'ZONE T="d=', delta, '", I=', nj                          

   DO i = 1, n                                                            
      IF (p(i) > p_curve(i))  WRITE(idf,*) v(i), p(i)                      
   ENDDO                                                                  


   ! Writing waves in vp plane

   ALLOCATE (pRH(n), QIS(n))

   DO i = 1, SIZE(wave_1 % type)

      v0 = wave_1 % v(i, 1)
      p0 = wave_1 % p(i, 1)

      dv = (wave_1 % v(i,2) - wave_1 % v(i,1)) / (n - 1)
   
      v(1) = wave_1 % v(i,1)
                
      DO j = 1, n
     
         v(j) = v(1) + (j-1) * dv      
   
         IF (wave_1 % type(i) == 'S') THEN
            pRH(j) = p_RH(v(j), v0, p0)
         ELSE      
            QIS(j) = Q(v(j), v0, p0)
         ENDIF
     
      ENDDO
         
   
      IF (wave_1 % type(i) == 'S') THEN
     
         WRITE (idf,*)  ''
         WRITE (idf,*)  'ZONE T="RH W1", I=', n       
         DO j = 1, n
            WRITE(idf,*) v(j), pRH(j)
         ENDDO  

      ELSE 
     
         WRITE (idf,*)  ''           
         WRITE (idf,*)  'ZONE T="Q W1", I=', n
         DO j = 1, n
            WRITE(idf,*) v(j), QIS(j)
         ENDDO 

      ENDIF

   ENDDO


   WRITE (idf,*) ''
   WRITE (idf,*) 'ZONE T="Contact", I=2'
   WRITE (idf,*) wave_1 % v(SIZE(wave_1 % type), 2), wave_1 % p(SIZE(wave_1 % type), 2)
   WRITE (idf,*) wave_3 % v(SIZE(wave_3 % type), 2), wave_3 % p(SIZE(wave_3 % type), 2)

   DO i = 1, SIZE(wave_3 % type)

      v0 = wave_3 % v(i, 1)
      p0 = wave_3 % p(i, 1)

      dv = (wave_3 % v(i,2) - wave_3 % v(i,1)) / (n - 1)
   
      v(1) = wave_3 % v(i,1)
                
      DO j = 1, n
     
         v(j) = v(1) + (j-1) * dv      
   
         IF (wave_3 % type(i) == 'S') THEN
            pRH(j) = p_RH(v(j), v0, p0)
         ELSE      
            QIS(j) = Q(v(j), v0, p0)
         ENDIF
     
      ENDDO
         
   
      IF (wave_3 % type(i) == 'S') THEN
     
         WRITE (idf,*)  ''
         WRITE (idf,*)  'ZONE T="RH W3", I=', n
         DO j = 1, n
            WRITE(idf,*) v(j), pRH(j)
         ENDDO  

      ELSE 
     
         WRITE (idf,*)  ''           
         WRITE (idf,*)  'ZONE T="Q W3", I=', n
         DO j = 1, n
            WRITE(idf,*) v(j), QIS(j)
         ENDDO 
 
      ENDIF

   ENDDO


   x_min = MIN( MINVAL(wave_1 % v), MINVAL(wave_3 % v) )
   x_max = MAX( MAXVAL(wave_1 % v), MAXVAL(wave_3 % v) )

   y_min = MIN( MINVAL(wave_1 % p), MINVAL(wave_3 % p) )
   y_max = MAX( MAXVAL(wave_1 % p), MAXVAL(wave_3 % p) )


   OPEN(UNIT=idf, FILE='TEC__vp_waves.lay', STATUS='unknown')                   

   WRITE(idf,1000)
   1000 FORMAT (t1,'#!MC 1000')

   WRITE(idf,*) '$!FRAMELAYOUT SHOWHEADER = NO'
   WRITE(idf,*) '$!REDRAWALL'
   WRITE(idf,*) '$!READDATASET  ''"TEC__vp_waves.plt" '''
   WRITE(idf,*) '$!XYLINEAXIS YDETAIL 1 {RANGEMIN = ',  y_min, '}'
   WRITE(idf,*) '$!XYLINEAXIS YDETAIL 1 {RANGEMAX = ',  y_max, '}'
   WRITE(idf,*) '$!XYLINEAXIS XDETAIL 1 {RANGEMIN = ',  x_min, '}'
   WRITE(idf,*) '$!XYLINEAXIS XDETAIL 1 {RANGEMAX = ',  x_max, '}'
   WRITE(idf,*) '$!LINEPLOTLAYERS SHOWSYMBOLS = YES'
   WRITE(idf,*) '$!ACTIVELINEMAPS += [1-4]'
   WRITE(idf,*) '$!LINEMAP [1-4]  LINES{SHOW = NO}'
   WRITE(idf,*) '$!LINEMAP [1]  SYMBOLS{SYMBOLSHAPE{GEOMSHAPE = LTRI}}'
   WRITE(idf,*) '$!LINEMAP [1]  SYMBOLS{COLOR = BLACK}'
   WRITE(idf,*) '$!LINEMAP [1]  SYMBOLS{SIZE = 2}'
   WRITE(idf,*) '$!LINEMAP [1]  SYMBOLS{LINETHICKNESS = 0.3}'
   WRITE(idf,*) '$!LINEMAP [2]  SYMBOLS{SYMBOLSHAPE{GEOMSHAPE = CIRCLE}}'
   WRITE(idf,*) '$!LINEMAP [2]  SYMBOLS{COLOR = RED}'
   WRITE(idf,*) '$!LINEMAP [2]  SYMBOLS{SIZE = 1}'
   WRITE(idf,*) '$!LINEMAP [2]  SYMBOLS{LINETHICKNESS = 0.3}'
   WRITE(idf,*) '$!LINEMAP [3]  SYMBOLS{SYMBOLSHAPE{GEOMSHAPE = CIRCLE}}'
   WRITE(idf,*) '$!LINEMAP [3]  SYMBOLS{COLOR = RED}'
   WRITE(idf,*) '$!LINEMAP [3]  SYMBOLS{SIZE = 1}'
   WRITE(idf,*) '$!LINEMAP [3]  SYMBOLS{LINETHICKNESS = 0.3}'
   WRITE(idf,*) '$!LINEMAP [4]  SYMBOLS{SYMBOLSHAPE{GEOMSHAPE = RTRI}}'
   WRITE(idf,*) '$!LINEMAP [4]  SYMBOLS{COLOR = BLACK}'
   WRITE(idf,*) '$!LINEMAP [4]  SYMBOLS{SIZE = 2}'
   WRITE(idf,*) '$!LINEMAP [4]  SYMBOLS{LINETHICKNESS = 0.3}'

   WRITE(idf,*) '$!ACTIVELINEMAPS += [5-6]'
   WRITE(idf,*) '$!LINEMAP [5-6] SYMBOLS{SHOW = NO}'
   WRITE(idf,*) '$!LINEMAP [5]  LINES{SHOW = YES}'
   WRITE(idf,*) '$!LINEMAP [5]  LINES{COLOR = BLACK}'
   WRITE(idf,*) '$!LINEMAP [5]  LINES{LINEPATTERN = SOLID}'
   WRITE(idf,*) '$!LINEMAP [5]  LINES{LINETHICKNESS = 0.5}'
   WRITE(idf,*) '$!LINEMAP [6]  LINES{SHOW = YES}'
   WRITE(idf,*) '$!LINEMAP [6]  LINES{COLOR = RED}'
   WRITE(idf,*) '$!LINEMAP [6]  LINES{LINEPATTERN = SOLID}'
   WRITE(idf,*) '$!LINEMAP [6]  LINES{LINETHICKNESS = 0.5}'
   
   nm = 6                                                   

   DO i = 1, SIZE(wave_1 % type)                                                

      IF (wave_1 % type(i) == 'S') THEN                                         
         nm = nm + 1                                                            
         WRITE(idf,*) '$!ACTIVELINEMAPS += [', nm, ']'
         WRITE(idf,*) '$!LINEMAP [', nm,']  SYMBOLS{SHOW = NO}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{SHOW = YES}'                  
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{COLOR = BLACK}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINEPATTERN = DOTTED}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINETHICKNESS = 0.1}'
      ELSE                                                                      
         nm = nm + 1                                                            
         WRITE(idf,*) '$!ACTIVELINEMAPS += [', nm, ']'
         WRITE(idf,*) '$!LINEMAP [', nm,']  SYMBOLS{SHOW = NO}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{SHOW = YES}'                  
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{COLOR = BLACK}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINEPATTERN = SOLID}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINETHICKNESS = 0.1}'
      ENDIF                                                                     
                                                                                
   ENDDO                                                                        

   nm = nm + 1
   WRITE(idf,*) '$!ACTIVELINEMAPS += [', nm, ']'
   WRITE(idf,*) '$!LINEMAP [', nm,']  SYMBOLS{SHOW = NO}'
   WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{SHOW = YES}'                  
   WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{COLOR = BLACK}'
   WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINEPATTERN = DASHED}'
   WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINETHICKNESS = 0.1}'


   DO i = 1, SIZE(wave_3 % type)                                                

      IF (wave_3 % type(i) == 'S') THEN                                         
         nm = nm + 1                                                            
         WRITE(idf,*) '$!ACTIVELINEMAPS += [', nm, ']'
         WRITE(idf,*) '$!LINEMAP [', nm,']  SYMBOLS{SHOW = NO}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{SHOW = YES}'                  
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{COLOR = BLACK}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINEPATTERN = DOTTED}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINETHICKNESS = 0.1}'
      ELSE                                                                      
         nm = nm + 1                                                            
         WRITE(idf,*) '$!ACTIVELINEMAPS += [', nm, ']'
         WRITE(idf,*) '$!LINEMAP [', nm,']  SYMBOLS{SHOW = NO}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{SHOW = YES}'                  
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{COLOR = BLACK}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINEPATTERN = SOLID}'
         WRITE(idf,*) '$!LINEMAP [', nm,']  LINES{LINETHICKNESS = 0.1}'
      ENDIF                                                                     
                                                                                
   ENDDO                                                                        
   
   WRITE(idf,*) '$!REDRAWALL'                                                   
   
   CLOSE (idf)  


END SUBROUTINE Plot_vp_Waves





FUNCTION Linear_Interpolation(gg, xx, x) RESULT(gx)

   ! Compute the value of function gg at point x. gg  is  a vector  of  the
   ! values of the function evaluated at the points specified  in  xx.  The
   ! value g(x) is computed by mean of linear interpolation.
   !
   ! This algorith assumes that the grid values xx(i) are ardered in incre-
   ! asing manner

   !------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: gg, xx
   REAL(KIND=8),               INTENT(IN) :: x

   LOGICAL :: check
   INTEGER :: i

   REAL(KIND=8) :: x_l, x_r,  &
                  g_l, g_r, gx
   !------------------------------------------------------------------------

   check = .FALSE.   ! Manage the possibility that xp is not included
                     ! in the range xx(1) - xx( SIZE(xx) )

   DO i = 1, SIZE(xx) - 1

      IF ( (x >= xx(i)  .AND.  x <  xx(i+1)) .OR. &
           (x <  xx(i)  .AND.  x >= xx(i+1)) ) THEN

          x_l = xx(i);   x_r = xx(i+1)
          g_l = gg(i);   g_r = gg(i+1)

          check = .TRUE.

          EXIT

      ENDIF

   ENDDO

   IF (.NOT. check)  THEN

      PRINT*, ' The point x is external to the interval of the grid'
      PRINT*, ' Linear_Interpolation function is STOPPED.'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
        WRITE(idf,*) ' Linear_Interpolation function is STOPPED.'
      CLOSE (idf)
      STOP

   ENDIF

   ! Linear interpolation for computing gg(x)
   gx =  g_l + (x - x_l) * (g_r - g_l) / (x_r - x_l)

END FUNCTION Linear_Interpolation





SUBROUTINE Plot_Char_Field (wave_1, wave_3, tt,   vacuum)  
   !------------------------------------------------------------------------
   IMPLICIT NONE
 
   TYPE(nonlin_wave), INTENT(IN) :: wave_1, wave_3
   REAL(KIND=8),      INTENT(IN) :: tt

   LOGICAL,           INTENT(IN), OPTIONAL :: vacuum
  
   INTEGER :: i, j, n
  
   REAL(KIND=8) sp, lb, rb, dpts
   !------------------------------------------------------------------------

   OPEN (UNIT = idf, FILE = 'TEC__char_field.plt', STATUS = 'unknown')      

   WRITE (idf,*)  'TITLE = "Characteristic Field"'	         
   WRITE (idf,*)  'VARIABLES = "x", "t"'		         

   DO i = 1, SIZE(wave_1 % type)			         

      IF ( wave_1 % type(i) == 'S' ) THEN		         
 
         sp = wave_1 % speeds(i,1) * tt  		         
 
    	   WRITE(idf,*) '' 				         
    	   WRITE (idf,*) 'ZONE T="Shock", I = ', 2 	         

    	   WRITE(idf,*) '0.0    0.0'			         
    	   WRITE(idf,*) sp, tt				         
    	   WRITE(idf,*) '' 				         
 
      ENDIF						         
 
 
      IF ( wave_1 % type(i) == 'F' ) THEN		         

         lb = wave_1 % speeds(i,1) * tt			         
         rb = wave_1 % speeds(i,2) * tt			         
 
         dpts = (rb - lb)/7.d0				         
 
         DO j = 0, 7					         
 
    	      WRITE(idf,*) ''				         
    	      WRITE (idf,*) 'ZONE T="Fan", I = ', 2	         

    	      WRITE(idf,*) lb + j*dpts, tt 		         
      	    WRITE(idf,*) '0.0   0.0'			         
    	      WRITE(idf,*) ''				         
 
    	   ENDDO						         
 
    	   WRITE(idf,*) '' 				         

      ENDIF						         
 
   ENDDO						         
 
 
   ! Contact Discontinuity				         
   IF (.NOT. PRESENT(vacuum)) THEN			         

       WRITE(idf,*) ''  				         
       WRITE (idf,*) 'ZONE T="Discontinuity", I = ', 2           

       WRITE(idf,*) '0.0    0.0'			         
       WRITE(idf,*) wave_3 % u(1,1) * tt,   tt  	         

   ENDIF						         
    							         
 
   DO i = 1, SIZE(wave_3 % type)			         

      IF (wave_3 % type(i) == 'S') THEN		         
 
         sp = wave_3 % speeds(i,1) * tt  		         
 
    	   WRITE(idf,*) '' 				         
         WRITE (idf,*) 'ZONE T="Shock", I = ', 2 	         

    	   WRITE(idf,*) '0.0    0.0'			         
         WRITE(idf,*) sp, tt				         
    	   WRITE(idf,*) '' 				         
 
      ENDIF						         
 
 
      IF (wave_3 % type(i) == 'F') THEN		         

         lb = wave_3 % speeds(i,1) * tt			         
         rb = wave_3 % speeds(i,2) * tt			         
 
         dpts = (rb - lb)/7.d0				         
 
    	   DO j = 0, 7					         
 
    	      WRITE(idf,*) ''				         
    	      WRITE (idf,*) 'ZONE T="Fan", I = ', 2	         

    	      WRITE(idf,*) lb + j*dpts, tt 		         
    	      WRITE(idf,*) '0.0   0.0'			         
    	      WRITE(idf,*) ''				         
 
    	   ENDDO						         
 
    	   WRITE(idf,*) '' 				         

      ENDIF						         
 
   ENDDO						         
 
   CLOSE (idf)


   OPEN(UNIT=idf, FILE='TEC__char_field.lay', STATUS='unknown')		              

   WRITE(idf,1000)                               
   1000 FORMAT (t1,'#!MC 1000')                            

   WRITE(idf,*) '$!FRAMELAYOUT SHOWHEADER = NO'                        
   WRITE(idf,*) '$!REDRAWALL'                              
   WRITE(idf,*) '$!READDATASET  ''"TEC__char_field.plt" '''        
   WRITE(idf,*) '$!ACTIVELINEMAPS += [1-50]'           
   WRITE(idf,*) '$!LINEMAP [1-50]  LINES{COLOR = BLACK}'       
   WRITE(idf,*) '$!VIEW FIT'               
                         
   n = 0                   
                         
   DO i = 1, SIZE(wave_1 % type)             
                         
      IF (wave_1 % type(i) == 'S') THEN    
         n = n + 1
         WRITE(idf,*) '$!LINEMAP [', n, ']  LINES{LINETHICKNESS = 0.5}'        
      ELSE 
         n = n + 8    
      ENDIF
       
   ENDDO                   

   IF (.NOT. PRESENT(vacuum)) THEN             

      ! Contact Discontinuity               
      n = n + 1                 
      WRITE(idf,*) '$!LINEMAP [', n, ']  LINES{LINETHICKNESS = 0.1}'
      WRITE(idf,*) '$!LINEMAP [', n, ']  LINES{LINEPATTERN = DASHED}'     

   ENDIF                   


                         
   DO i = 1, SIZE(wave_3 % type)             
                         
      IF (wave_3 % type(i) == 'S') THEN    
         n = n + 1
         WRITE(idf,*) '$!LINEMAP [', n, ']  LINES{LINETHICKNESS = 0.5}'        
      ELSE 
         n = n + 8    
      ENDIF
       
   ENDDO                   
                        
                         
   WRITE(idf,*) '$!REDRAWALL'                



   IF (PRESENT(vacuum)) THEN               

      WRITE(idf,*) '$!ATTACHGEOM'             
      WRITE(idf,*) '  ANCHORPOS'              
      WRITE(idf,*) '  {'                
      WRITE(idf,*) '  X = 0'                
      WRITE(idf,*) '  Y = 0'                
      WRITE(idf,*) '  }'                
      WRITE(idf,*) '  COLOR = CUSTOM2'              
      WRITE(idf,*) '  ISFILLED = YES'             
      WRITE(idf,*) '  FILLCOLOR = CUSTOM2'            
      WRITE(idf,*) '  DRAWORDER = BEFOREDATA'           
      WRITE(idf,*) '  RAWDATA'                
      WRITE(idf,*) '1'                  
      WRITE(idf,*) '4'                  
      WRITE(idf,*) '0 0'                
      WRITE(idf,*) wave_1 % speeds(SIZE(wave_1 % type),2) * tt, tt      
      WRITE(idf,*) wave_3 % speeds(1,1) * tt, tt          
      WRITE(idf,*) '0 0'                
                         
   ENDIF

   CLOSE (idf)								          

END SUBROUTINE Plot_Char_Field
  
  
  
  
  
SUBROUTINE Plot_Profile (xx, uu, dt, name)  
 
   !------------------------------------------------------------------------
   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xx
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu 
   REAL(KIND=8),               INTENT(IN) :: dt
   CHARACTER(LEN=20),          INTENT(IN) :: name
 
   REAL(KIND=8) :: du

   INTEGER :: i
   !------------------------------------------------------------------------

   !===JLG
   OPEN(UNIT = idf, FILE = TRIM(name)//'.plt', STATUS = 'unknown') 
   DO i = 1, SIZE(uu)                          
      WRITE(idf,*) (xx(i)-minval(xx))/(maxval(xx)-minval(xx)), uu(i)                       
   ENDDO
   CLOSE(idf)
   !===JLG
   
   OPEN(UNIT = idf, FILE = 'TEC__'//TRIM(name)//'.plt', STATUS = 'unknown') 	        
 
   WRITE(idf,*) 'TITLE = "'//name//'"'                     
   WRITE(idf,*) 'VARIABLES = "x", "'//TRIM(name)//'"'                  
 

   WRITE(idf,*) 'ZONE T ="'//TRIM(name), dt,'", N=', SIZE(uu), &       
        ', E =', SIZE(uu) - 1,       &       
        ', DATAPACKING=POINT, ZONETYPE=FELINESEG'              

   DO i = 1, SIZE(uu)                          
      WRITE(idf,*) xx(i), uu(i)                       
   ENDDO                           
 
   DO i = 1, SIZE(uu) - 1                        
      WRITE(idf,*) i, i + 1                        
   ENDDO                           

   CLOSE(idf)									        


   du = MAXVAL(uu) - MINVAL(uu)							        

   OPEN(UNIT = idf, FILE = 'TEC__'//TRIM(name)//'.lay', STATUS = 'unknown') 	        

   WRITE(idf,1000)                       
   1000 FORMAT (t1,'#!MC 1000')                    

   WRITE(idf,*) '$!FRAMELAYOUT SHOWHEADER = NO'                
   WRITE(idf,*) '$!REDRAWALL'                      
   WRITE(idf,*) '$!READDATASET  ''"TEC__'//TRIM(name)//'.plt" '''            
   WRITE(idf,*) '$!TWODAXIS AXISMODE = INDEPENDENT'                
   WRITE(idf,*) '$!TWODAXIS YDETAIL{RANGEMIN = ', MINVAL(uu) - 0.3 * ABS(du), '}'      
   WRITE(idf,*) '$!TWODAXIS YDETAIL{RANGEMAX = ', MAXVAL(uu) + 0.3 * ABS(du), '}'      
   WRITE(idf,*) '$!TWODAXIS XDETAIL{RANGEMIN = ', xx(1), '}'             
   WRITE(idf,*) '$!TWODAXIS XDETAIL{RANGEMAX = ', xx(SIZE(uu)), '}'            
   WRITE(idf,*) '$!FIELD [1-', SIZE(uu),']  MESH{COLOR = BLACK}'           
   WRITE(idf,*) '$!ACTIVEFIELDZONES = [1]'                 
   WRITE(idf,*) '$!REDRAWALL'                      
  	
   CLOSE (idf)									        
 
END SUBROUTINE Plot_Profile  
 

  
  

SUBROUTINE Animate_Profiles (x, ww, t)
   !------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww 
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: x
   REAL(KIND=8),                 INTENT(IN) :: t

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xx
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: vv
 
   REAL(KIND=8) :: dt

   INTEGER :: i, j, nt_step = 100
   !------------------------------------------------------------------------


   ALLOCATE (vv(SIZE(x)))

   DO i = 1, SIZE(x)
      vv(i) = x(i) / t
   ENDDO


   ALLOCATE (xx(SIZE(x), 0:nt_step))

   dt = t / nt_step

   xx(1,:) = x(1);     xx(SIZE(x),:) = x(SIZE(x))


   DO j = 0, nt_step

      DO i = 2, SIZE(x) - 1
         xx(i,j) = vv(i) * j*dt
      ENDDO

   ENDDO


  OPEN(UNIT=idf, FILE='MOV__profile', STATUS='unknown')                 
 
   WRITE(idf,*) 'TITLE = "Animations"'             
   WRITE(idf,*) 'VARIABLES = "x", "w"'            
 
   DO j = 0, nt_step                

      WRITE(idf,*) 'ZONE T ="Specific Volume', dt,'", N=', SIZE(ww,2), &  
           ', E =', SIZE(ww,2) - 1,        &     
           ', DATAPACKING=POINT, ZONETYPE=FELINESEG'          

      DO i = 1, SIZE(ww,2)              
         WRITE(idf,*) xx(i,j), ww(1,i)                     
      ENDDO                       
 
      DO i = 1, SIZE(ww,2) - 1              
         WRITE(idf,*) i, i + 1              
      ENDDO                   

   ENDDO                  

   WRITE(idf,*) ''                

   DO j = 0, nt_step                

      WRITE(idf,*) 'ZONE T ="Velocity', dt,'", N=', SIZE(ww,2), &     
           ', E =', SIZE(ww,2) - 1,       &       
           ', DATAPACKING=POINT, ZONETYPE=FELINESEG'          

      DO i = 1, SIZE(ww,2)              
         WRITE(idf,*) xx(i,j), ww(2,i)                     
      ENDDO                       
 
      DO i = 1, SIZE(ww,2) - 1              
         WRITE(idf,*) i, i + 1              
      ENDDO                   

   ENDDO                 

   WRITE(idf,*) ''               

   DO j = 0, nt_step               

      WRITE(idf,*) 'ZONE T ="Pressure', dt,'", N=', SIZE(ww,2), &    
             ', E =', SIZE(ww,2) - 1,       &      
             ', DATAPACKING=POINT, ZONETYPE=FELINESEG'         

      DO i = 1, SIZE(ww,2)             
         WRITE(idf,*) xx(i,j), ww(3,i)                    
      ENDDO                      
 
      DO i = 1, SIZE(ww,2) - 1             
         WRITE(idf,*) i, i + 1             
      ENDDO                  

   ENDDO                 

   CLOSE (idf)

END SUBROUTINE Animate_Profiles


END MODULE plot_utils
