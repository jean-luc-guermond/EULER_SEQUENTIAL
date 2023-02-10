MODULE eigenvalues

   USE rhLocus_intCurve
   
   USE vdw_gas

   IMPLICIT NONE

   PRIVATE

   REAL(KIND=8), DIMENSION(2) :: w_i

   REAL(KIND=8) :: t, v_i, u_i, p_i, sigma

   INTEGER :: sign

   PUBLIC :: lambda_Init, lambda, fan_Init, fan


CONTAINS



SUBROUTINE lambda_Init (sign_, w_i_, t_)

   !------------------------------------------------------------------------  
   IMPLICIT NONE								    
 
   INTEGER,                    INTENT(IN) :: sign_
   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: w_i_
   REAL(KIND=8),               INTENT(IN) :: t_
   !------------------------------------------------------------------------   
 
   sign = sign_

   w_i = w_i_  

   t = t_

END SUBROUTINE lambda_Init





FUNCTION lambda_Vector(v, p, u)   RESULT(l)

   ! Computes the eigenvalue. The index is specified by input variable sign
				    
   !------------------------------------------------------------------------  
   IMPLICIT NONE								    
	    
   REAL(KIND=8), INTENT(IN) :: v, p, u				    
  									    
   REAL(KIND=8), DIMENSION(3) :: l

   REAL(KIND=8) :: c2, c
   !------------------------------------------------------------------------  

   c2 = (delta + 1) * (p*v**2 + 3) / (v - 1.d0/3)  -  6/v

   IF (c2 <= 0) THEN

      PRINT*, 'Negative square speed of sound '
      PRINT*, 'Computation is interrupted.'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
        WRITE(idf,*) 'Negative square speed of sound '
      CLOSE(idf)
      STOP

   ENDIF

   c = SQRT(c2)

   l(1) = u - c
							    
   l(2) = u 
 		
   l(3) = u + c
							    
END FUNCTION lambda_Vector
 									   



 									   
FUNCTION lambda(sign, v, p, u)   RESULT(l)

   ! Computes the eigenvalue. The index is specified by input variable sign
				    
   !------------------------------------------------------------------------  
   IMPLICIT NONE								    

   INTEGER,      INTENT(IN) :: sign	    
   REAL(KIND=8), INTENT(IN) :: v, p, u
 									    
   REAL(KIND=8) :: l

   REAL(KIND=8) :: c2, c
   !------------------------------------------------------------------------  

   c2 = (delta+1) * (p*v**2 + 3) / (v - 1.d0/3)  -  6/v
 
   IF (c2 <= 0) THEN

      PRINT*, 'Negative square speed of sound '
      PRINT*, 'Computation is interrupted.'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
        WRITE(idf,*) 'Negative square speed of sound '
      CLOSE(idf)
      STOP

   ENDIF

   c = SQRT(c2)

   l = u + sign * c
							    
END FUNCTION lambda





SUBROUTINE fan_Init(sign_, sigma_, vi_, ui_, pi_)

   !-----------------------------------------------------------------------   
   IMPLICIT NONE                                                             

   INTEGER,      INTENT(IN) :: sign_
   REAL(KIND=8), INTENT(IN) :: sigma_, vi_, ui_, pi_
   !-----------------------------------------------------------------------   
 
   sign = sign_

   sigma = sigma_

   v_i = vi_
   u_i = ui_
   p_i = pi_
 
END SUBROUTINE fan_Init



FUNCTION fan(v)   RESULT(l)			    
 
   !------------------------------------------------------------------------  
   IMPLICIT NONE								    

   REAL(KIND=8), INTENT(IN) :: v
 									    
   REAL(KIND=8) :: l

   REAL(KIND=8) :: c2, c, u
   !------------------------------------------------------------------------  

   c2 = (3*sigma * (delta+1) * v**2 ) / ((v - 1.d0/3)**(delta+2) * (1.d0/3)**(1-delta)) &
       - 6 / v 

   IF (c2 <= 0) THEN

      PRINT*, 'Negative square speed of sound in fan'
      PRINT*, 'Computation is interrupted.'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
        WRITE(idf,*) 'Negative square speed of sound '
      CLOSE(idf)
      STOP

   ENDIF

   c = SQRT(c2)

   u = u_i  +  Integral_Curve(sign, v_i, v, v_i, p_i)

   l = u + sign * c
						    
END FUNCTION fan



END MODULE eigenvalues
