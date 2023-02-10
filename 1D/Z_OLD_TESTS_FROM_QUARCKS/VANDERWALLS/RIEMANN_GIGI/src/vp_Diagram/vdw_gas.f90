
! =========================================================================
!   Description: Module that contains the definitions of the Van der Waals
!                Pressure equation of state P(v,T). Temperature is assumed  
!                as constant (isothermal gas). 
!
!   	 Author: Marco Fossati, Pierluigi Di Lizia
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!                        dilizia@aero.polimi.it
!
!          Year:  2005, November
! ==========================================================================

MODULE vdw_gas


   CONTAINS


   FUNCTION p_vdW__vt(v, t)   RESULT(p)

   ! Isothermal Van der Waals dimensionless equation of state for Pressure. 
   !
   !   p = p(v,t)

   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v, t

   REAL(KIND=8) :: p
   !----------------------------------------------------------------------- 

   p = -3 / v**2  +  8*t / (3*v - 1)

   END FUNCTION p_vdW__vt




  
   FUNCTION p_vdW__vs(v, sigma, delta)   RESULT(p)

   ! Van der Waals dimensionless equation of state for Pressure. 
   !
   !   p = p(v,s)

   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v, sigma, delta

   REAL(KIND=8) :: p
   !----------------------------------------------------------------------- 

   p = 27.d0*sigma / (3.d0*v - 1)**(delta + 1) - 3.d0 / v**2

   END FUNCTION p_vdW__vs





   FUNCTION Gzero_Locus(v, delta)   RESULT(p_G0)

   ! Locus of G = 0, where G is the fundamental derivative of gasdynamics

   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v, delta

   REAL(KIND=8) :: p_G0
   !----------------------------------------------------------------------- 

   p_G0 = 3.d0/v**2 * (6.d0/((delta+1)*(delta+2)) * (1 - 1/(3*v))**2 - 1) 

   END FUNCTION Gzero_Locus


  
END MODULE vdw_gas
