
! ============================================================================
!   Description: Parameters for the problem. 
!
!   	 Author: Marco Fossati
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!
!          Year: 2005, October
! =============================================================================

MODULE problem_params

   !-------------------------------------------------------------------------
   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(3) :: w_L, w_R

   REAL(KIND=8) :: time, eps, delta, relax

   INTEGER, PARAMETER :: idf = 11, Max_Iter = 100000

   INTEGER :: fan_points
  
   CHARACTER, DIMENSION(:), ALLOCATABLE :: type_1, type_3
   !-------------------------------------------------------------------------  


CONTAINS


SUBROUTINE Load_Param

   OPEN ( UNIT = idf, FILE = 'problem.param', STATUS = 'old' )  
 
   READ(idf,*)
   READ(idf,*) delta
   READ(idf,*)
   READ(idf,*)
   READ(idf,*) w_L
   READ(idf,*) w_R
   READ(idf,*)
   READ(idf,*)
   READ(idf,*) fan_points
   READ(idf,*) time
   READ(idf,*)
   READ(idf,*)
   READ(idf,*) relax
   READ(idf,*) eps

 
   CLOSE (idf)


   !OPEN ( UNIT = idf, FILE = 'data.param', STATUS = 'old' )  
 
   !READ(idf,*) w_L
   !READ(idf,*) w_R
 
   CLOSE (idf)

END SUBROUTINE Load_Param





SUBROUTINE Check_Data (sigma_L, sigma_R, s_star)

   !-------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: sigma_L, sigma_R, s_star
  
   REAL(KIND=8) :: v_L, u_L, p_L,  &
                   v_R, u_R, p_R  
   !-------------------------------------------------------------------------

   IF (w_L(1) == w_R(1)) THEN
     
      ! The problem with equal specific volumes is due to the
      ! fact that guesses for  Newton's Method  are computed
      ! referring to the isothermal case ???

      PRINT*, ''
      PRINT*, 'Error. CHECK_DATA:'
      PRINT*, 'Equal specific volumes.'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
      WRITE(idf,*) w_L
      WRITE(idf,*) w_R
      WRITE(idf,*) 'Equal specific volumes.'
      CLOSE(idf)
      PRINT*, ''
      STOP

   ENDIF


   IF (w_L(2) == w_R(2)  .AND.  w_L(3) == w_R(3)) THEN
     
      ! Initial data is contact discontinuity no further processing
      PRINT*, ''
      PRINT*, 'Error. CHECK_DATA:'
      PRINT*, 'Initial data is contact discontinuity'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
      WRITE(idf,*) w_L
      WRITE(idf,*) w_R
      WRITE(idf,*) 'Initial data is contact discontinuity'
      CLOSE(idf)
      PRINT*, ''
      STOP

   ENDIF


   v_L = w_L(1);    v_R = w_R(1)             
   u_L = w_L(2);    u_R = w_R(2)             
   p_L = w_L(3);    p_R = w_R(3)             


   IF (Inside_Coexistence(p_L, v_L)) THEN            
      PRINT*, ''
      PRINT*, 'Error. The state (p_L, v_L) falls under the coexistence curve'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')           
      WRITE(idf,*) w_L             
      WRITE(idf,*) w_R             
      WRITE (idf,*) 'The state (p_L, v_L) falls under the coexistence curve' 
      CLOSE(idf)               
      PRINT*, ''               
      STOP                 
   ENDIF                   
                        
   IF (Inside_Coexistence(p_R, v_R)) THEN            
      PRINT*, ''
      PRINT*, 'Error. The state (p_R, v_R) falls under the coexistence curve'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')           
      WRITE(idf,*) w_L             
      WRITE(idf,*) w_R             
      WRITE (idf,*) 'The state (p_R, v_R) falls under the coexistence curve' 
      CLOSE(idf)               
      PRINT*, ''               
      STOP                 
   ENDIF                   


   WRITE(*,1000) delta                 
   1000 FORMAT ('     Gas Complexity =', e13.5)          
   WRITE(*,1001) s_star                
   1001 FORMAT ('     Starred  Sigma =', e13.5)          

                         
   PRINT*, ''                    
   PRINT*, '              LEFT  STATE       RIGHT STATE'       
   PRINT*, ''                    
   WRITE(*, '(8x, a5, e13.5, 5x, e13.5)') 'v', v_L, v_R        
   WRITE(*, '(8x, a5, e13.5, 5x, e13.5)') 'u', u_L, u_R        
   WRITE(*, '(8x, a5, e13.5, 5x, e13.5)') 'p', p_L, p_R        
   WRITE(*, '(8x, a5, e13.5, 5x, e13.5)') 'nu', 3*v_L, 3*v_R       
   WRITE(*, '(8x, a5, e13.5, 5x, e13.5)') 'Pi', p_L/27, p_R/27       
   WRITE(*, '(8x, a5, e13.5, 5x, e13.5)') 'sigma', sigma_L, sigma_R      
   PRINT*, ''
  
   IF (sigma_L < s_star) THEN
       PRINT*, '    NONCONVEX  Left  Isentrope'
   ELSE
       PRINT*, '    CONVEX     Left  Isentrope'
   ENDIF

   IF (sigma_R < s_star) THEN
       PRINT*, '    NONCONVEX  Right Isentrope'
   ELSE
       PRINT*, '    CONVEX     Right Isentrope'
   ENDIF


END SUBROUTINE Check_Data





FUNCTION  Inside_Coexistence(P, v)  RESULT(i_c)	        	      
 
   !-------------------------------------------------------------------------  
   IMPLICIT NONE 			        			      
  					        			      
   REAL(KIND=8), INTENT(IN) :: P, v	        			      
   LOGICAL :: i_c  	        			      
  		  
   INTEGER, PARAMETER :: JA = 1000
			        			      
   REAL(KIND=8), DIMENSION(JA), SAVE :: vlA, vgA
   			      
   REAL(KIND=8) :: tummy, pummy
  
   INTEGER :: j  			        			      
  
   LOGICAL, SAVE :: first = .TRUE.	        			      
   !-------------------------------------------------------------------------  							      		      
  							      		      
   IF (first) THEN					      		      
  							      		      
      OPEN (UNIT = 35, FILE = 'curva_di_saturazione.dat', & 		      
            FORM = 'formatted', STATUS = 'old') 	    		      
         						    		      
      first = .FALSE.					    		      
        						    		      
      READ (35,*)					    		      
        						    		      
      DO j = 1, JA					    		      
         READ (35,*) tummy, pummy, vgA(j), vlA(j)	    		      
      ENDDO						    		      

      CLOSE (35)				            		      
  							      		      
   ENDIF						      		      
  							      		      

   i_c = .FALSE.                         
                               
   IF (P > 1) RETURN                         
                               
   j = (JA - 1)*(1 - P) + 1                        
                               
   IF (v < vlA(j)  .OR.  v > vgA(j)) RETURN                    
                               
   i_c = .TRUE.                          

END FUNCTION  Inside_Coexistence			        	      

  
END MODULE problem_params
