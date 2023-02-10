
! =========================================================================
!   Description: Module for the computation of the solution for pure waves
!                in the case of isothermal P-System  with  Van  der  Waals
!                equation of state for pressure, p(v;T).
!
!                u_RH(sign, v, v0, t)       -> Hugoniot Locus through v0
!                Integral_Curve(sign, v, v0, t) -> Integral Curve through v0
!
!
!   	 Author: Marco Fossati
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!
!          Year: 2005, November
! ==========================================================================

MODULE rhLocus_intCurve

   USE problem_params
   USE numerics
   USE vdw_gas


CONTAINS


FUNCTION u_RH(wsign, v, v0, p0)   RESULT(u)
   !------------------------------------------------------------------------ 
   IMPLICIT NONE
  
   INTEGER,      INTENT(IN) :: wsign
   REAL(KIND=8), INTENT(IN) :: v, v0, p0

   REAL(KIND=8) :: u
   !------------------------------------------------------------------------

   IF ((v - v0) * (p_RH(v, v0, p0) - p0) > 0) THEN

      PRINT*, 'Error. U_RH:'
      PRINT*, 'Negative argument of square root'
      PRINT*, ''
      OPEN(UNIT=idf, FILE='stopped')
      WRITE(idf,*) w_L
      WRITE(idf,*) w_R
      WRITE(idf,*) 'Negative argument of square root in u_RH'
      CLOSE(idf)
      STOP
  
   ENDIF


   u = - wsign * SIGN(1.d0,(v - v0)) &
               * SQRT(-(v - v0) * (p_RH(v, v0, p0) - p0))

END FUNCTION u_RH





FUNCTION p_RH(v, v0, p0)   RESULT(p)

   ! Rankine-Hugoniot Adiabatic transformation (nonisentropic)

   !------------------------------------------------------------------------ 
   IMPLICIT NONE
  
   REAL(KIND=8), INTENT(IN) :: v, v0, p0

   REAL(KIND=8) :: p

   REAL(KIND=8) :: e0, num, den  
   !------------------------------------------------------------------------

   e0 = (1/delta) * (p0 + 3/v0**2) * (v0 - 1.d0/3) - 3/v0

   num = e0 - p0*(v - v0)/2 + 3*(1 - 1/delta)/v + 1/(delta*v**2)
  
   den = (1.d0/2 + 1/delta)*v - (v0/2 + 1/(3*delta))

   p = num/den

END FUNCTION p_RH





FUNCTION Dp_RH(v, v0, p0)   RESULT(Dp)
   !------------------------------------------------------------------------ 
   IMPLICIT NONE
  
   REAL(KIND=8), INTENT(IN) :: v, v0, p0

   REAL(KIND=8) :: Dp

   REAL(KIND=8) :: A, B, C, e0, num, den, Dnum, Dden 
   !------------------------------------------------------------------------

   A = 3 - 3/delta

   B = 1.d0/2 + 1/delta

   C = v0/2 + 1/(3*delta)


   e0 = (1/delta) * (p0 + 3/v0**2) * (v0 - 1.d0/3) - 3/v
  

    num = e0 - p0*(v-v0)/2 + A/v + 1/(delta*v**2)

   Dnum = -p0/2 - A/v**2 - 2/(delta*v**3) 


    den = B*v - C

   Dden = B


   Dp = (Dnum * den  -  num * Dden) / den**2

END FUNCTION Dp_RH





FUNCTION Integral_Curve(wsign, v1, v2, v0, p0)   RESULT(vel)
   !------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER,      INTENT(IN) :: wsign
   REAL(KIND=8), INTENT(IN) :: v1, v2, v0, p0

   REAL(KIND=8) :: vel
   !------------------------------------------------------------------------ 

   vel = 0.d0

   CALL Init_IC(p0, v0)
            
   CALL ODE_int(v1, v2,  vel,  SQRT_DQ)

   vel = - wsign * vel

END FUNCTION Integral_Curve


END MODULE rhLocus_intCurve
