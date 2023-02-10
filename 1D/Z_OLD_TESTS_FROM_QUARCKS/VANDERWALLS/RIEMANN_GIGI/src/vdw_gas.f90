
! =========================================================================
!   Description: Module that contains the definitions of the Van der Waals
!                Pressure equation of state P(v,T).
!
!   	 Author: Marco Fossati
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!
!          Year:  2005, November
! ==========================================================================

MODULE vdw_gas


   USE problem_params                                                        

   !------------------------------------------------------------------------  
   IMPLICIT NONE

   REAL(KIND=8), PARAMETER :: t_star = 1.0678711d0, &
                              delta_lim = 0.06d0

   REAL(KIND=8) :: sigma_i, Pi_i,   nu_i,   &
                            Pi_hat, nu_hat, &
                   nu0, sigma00, p00, v00
   !------------------------------------------------------------------------  


CONTAINS


FUNCTION p_vdW__vt(v, t)   RESULT(p)  
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v, t

   REAL(KIND=8) :: p
   !----------------------------------------------------------------------- 

   p = -3 / v**2  +  8*t / (3*v - 1)

END FUNCTION p_vdW__vt





FUNCTION D_p_vdW__vt(v, t)   RESULT(D_p) ! partial derivative with respect to v
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v, t

   REAL(KIND=8) :: D_p
   !----------------------------------------------------------------------- 

   D_p = 6 / v**3  -  (24*t) / (3*v - 1)**2

END FUNCTION D_p_vdW__vt





FUNCTION DD_p_vdW__vt(v, t)	RESULT(DD_p) ! second partial derivative with respect to v
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v, t

   REAL(KIND=8) :: DD_p
   !----------------------------------------------------------------------- 

   DD_p = - 18 / v**4  +  (144*t) / (3*v - 1)**3

END FUNCTION DD_p_vdW__vt
  




FUNCTION sigma__Pinu(Pi, nu)   RESULT(sigma)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: Pi, nu

   REAL(KIND=8) :: sigma
   !----------------------------------------------------------------------- 

   sigma = (Pi + 1/nu**2) * (nu - 1)**(delta+1)

END FUNCTION sigma__Pinu





FUNCTION t_vdW__pv(p, v)   RESULT(t)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: p, v

   REAL(KIND=8) :: t
   !----------------------------------------------------------------------- 
   
   t = (3.d0/8) * (p + 3/v**2) * (v - 1.d0/3) 

END FUNCTION t_vdW__pv





FUNCTION sigma_star(delta)   RESULT(s_star)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: delta

   REAL(KIND=8) :: s_star
   !----------------------------------------------------------------------- 

   s_star = (3.d0/2**7) * (3 + delta)**(3 + delta) &
                        * (1 - delta)**(1 - delta) &
                        / ((1 + delta) * (2 + delta))

END FUNCTION sigma_star  




!**************************************************************************
!  SUBPROGRAMS FOR INFLECTION POINTS



FUNCTION Infl(nu)   RESULT(i)
   !-----------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: i
   !-----------------------------------------------------------------------

   i = (nu - 1)**(3 + delta) / nu**4

   END FUNCTION Infl



FUNCTION D_Infl(nu)   RESULT(D_i)
   !-----------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: D_i
   !-----------------------------------------------------------------------

   D_i = (4 - (1 - delta)*nu) * (nu - 1)**(2 + delta) / nu**5

END FUNCTION D_Infl




!**************************************************************************
!  SUBPROGRAMS FOR ABSOLUTE ENVELOPE



SUBROUTINE Init_Pi_vdW__sigmanu (sigma_i_)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: sigma_i_
   !----------------------------------------------------------------------- 

   sigma_i = sigma_i_

END SUBROUTINE Init_Pi_vdW__sigmanu




FUNCTION Pi_vdW__sigmanu(nu)   RESULT(Pi)

   ! Equation of state giving the pressure for isentropic transformations

   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: Pi
   !----------------------------------------------------------------------- 

   Pi = sigma_i / (nu - 1)**(1 + delta)  -  1 / nu**2

END FUNCTION Pi_vdW__sigmanu




FUNCTION DPi_vdW__sigmanu(nu)   RESULT(DPi)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: DPi
   !----------------------------------------------------------------------- 

   DPi = - (1 + delta) * sigma_i / (nu - 1)**(2 + delta)  +  2 / nu**3

END FUNCTION DPi_vdW__sigmanu




FUNCTION DDPi_vdW__sigmanu(nu)   RESULT(DDPi)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: DDPi
   !----------------------------------------------------------------------- 

   DDPi = (1 + delta) * (2 + delta) * sigma_i / (nu - 1)**(3 + delta)  -  6 / nu**4

END FUNCTION DDPi_vdW__sigmanu





!**************************************************************************
!  SUBPROGRAMS FOR TANGENT LINES



SUBROUTINE Init_psi (nu0_, sigma00_)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu0_, sigma00_
   !----------------------------------------------------------------------- 
 
   nu0 = nu0_
   
   sigma00 = sigma00_

END SUBROUTINE Init_psi





FUNCTION psi(nu)   RESULT(ps)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: ps

   REAL(KIND=8) :: nu_nu0, nu_1, nu0_1
   !----------------------------------------------------------------------- 

   nu_nu0 = nu - nu0
   nu_1   = nu - 1
   nu0_1  = nu0 - 1
   
   ps =   (1/nu_nu0**2) * (1/nu0_1**(1+delta)  -  1/nu_1**(1+delta))  &
        - (1+delta) / (nu_nu0 * nu_1**(2+delta))                      &
        - (nu + 2*nu0) / (sigma00 * nu0**2 * nu**3)

END FUNCTION psi





FUNCTION Dpsi(nu)   RESULT(Dps)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: Dps

   REAL(KIND=8) :: nu_nu0, nu_1, nu0_1
   !----------------------------------------------------------------------- 

   nu_nu0 = nu - nu0
   nu_1   = nu - 1
   nu0_1  = nu0 - 1
   
   Dps =   (2/nu_nu0**3) * (1/nu_1**(1+delta)  -  1/nu0_1**(1+delta))  &
        +  2*(1+delta) / (nu_nu0**2 * nu_1**(2+delta))                 &
        +  (1+delta)*(2+delta) / (nu_nu0 * nu_1**(3+delta))            &
        +  (2*nu + 6*nu0) / (sigma00 * nu0**2 * nu**4) 

END FUNCTION Dpsi





!**************************************************************************
!  SUBPROGRAMS FOR INTERSECTIONS



SUBROUTINE Init_phi (Pi_i_, nu_i_, Pi_hat_, nu_hat_)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: Pi_i_, nu_i_, Pi_hat_, nu_hat_
   !----------------------------------------------------------------------- 

   Pi_i = Pi_i_;   Pi_hat = Pi_hat_
   nu_i = nu_i_;   nu_hat = nu_hat_
 
END SUBROUTINE Init_phi




FUNCTION phi(nu)   RESULT(p)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: p

   REAL(KIND=8) :: A, B
   !----------------------------------------------------------------------- 
   
   B = (Pi_hat - Pi_i) / (nu_hat - nu_i)

   A = Pi_i - B * nu_i

   p = sigma_i/(nu - 1)**(1+delta)  -  1/nu**2  -  B*nu  -  A

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   p = p / ((nu - nu_hat) * (nu - nu_i))
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END FUNCTION phi




FUNCTION Dphi(nu)   RESULT(Dp)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: nu

   REAL(KIND=8) :: p, Dp

   REAL(KIND=8) :: A, B, NUM, DEN
   !-----------------------------------------------------------------------    
   
   B = (Pi_hat - Pi_i) / (nu_hat - nu_i)

   A = Pi_i  -  B * nu_i

   p = sigma_i/(nu - 1)**(1+delta)  -  1/nu**2  -  B*nu  -  A


   Dp = - (1+delta)*sigma_i/(nu - 1)**(2+delta)  +  2/nu**3  -  B

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   NUM = Dp*(nu - nu_hat) * (nu - nu_i)  -  p*(2*nu - nu_hat - nu_i)
 
   DEN = ((nu - nu_hat) * (nu - nu_i))**2

   Dp = NUM / DEN
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END FUNCTION Dphi



!**************************************************************************




FUNCTION h(r)   RESULT(h_)

   ! Second derivative for dimensionless Van der Waals Pressure with
   ! reduced density as indipendent variable. Such function will  be
   ! used for inflections points determination

   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: r

   REAL(KIND=8) :: h_
   !----------------------------------------------------------------------- 

   h_ = r * (3 - r)**3

END FUNCTION h





FUNCTION Dh(r)   RESULT(Dh_)

   ! Third derivative for dimensionless Van der Waals Pressure with
   ! reduced density as indipendent variable. Such function will be
   ! used for inflections points determination

   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: r

   REAL(KIND=8) :: Dh_
   !----------------------------------------------------------------------- 

   Dh_ = (3 - r)**2 * (3 - 4*r)

END FUNCTION Dh





FUNCTION Q(v, v0,p0)   RESULT(P)

   ! Equation of state giving the pressure as a function of specific
   ! volume for the isentrope passing through the state (v0, p0)
   
   !-----------------------------------------------------------------------       
   IMPLICIT NONE
   
   REAL(KIND=8), INTENT(IN) :: v, v0, p0
   
   REAL(KIND=8) :: P
   !----------------------------------------------------------------------- 

   P = (P0 + 3/v0**2) * ((v0 - 1.d0/3)/(v - 1.d0/3))**(1+delta)  -  3/v**2

END FUNCTION Q





FUNCTION DQ(v, v0,p0)   RESULT(dP)
   !----------------------------------------------------------------------- 
   IMPLICIT NONE
   
   REAL(KIND=8), INTENT(IN) :: v, v0, p0
   
   REAL(KIND=8) :: dP
   !-----------------------------------------------------------------------       

   dP = (P0 + 3/v0**2) * ((v0 - 1.d0/3)/(v - 1.d0/3))**(1+delta) 

   dP = - (1+delta) * dP / (v - 1.d0/3)  +  2*3/v**3

END FUNCTION DQ





SUBROUTINE  Init_IC (pressure, specific_volume)
   !-----------------------------------------------------------------------      
   IMPLICIT NONE
    
   REAL(KIND=8), INTENT(IN) :: pressure, specific_volume
   !----------------------------------------------------------------------- 
      
   p00 = pressure 
   v00 = specific_volume
       
END SUBROUTINE  Init_IC





FUNCTION SQRT_DQ(v, dum)   RESULT(f)

   ! Integrand function for integral curve.

   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v, dum

   REAL(KIND=8) :: f, dum_
   !----------------------------------------------------------------------- 

   dum_ = dum   !Dummy assignation to support ODE_integration 

   f = SQRT(-DQ(v,   v00, p00))

END FUNCTION SQRT_DQ





FUNCTION Gzero_Locus(v)   RESULT(p_G0)

   ! Locus of G = 0, where G is the fundamental derivative of gasdynamics

   !----------------------------------------------------------------------- 
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: v

   REAL(KIND=8) :: const, p_G0
   !----------------------------------------------------------------------- 
  
   const = 6 / ((1 + delta) * (2 + delta))

   p_G0 = (3/v**2) * (const * (1 - 1/(3*v))**2  -  1) 

END FUNCTION Gzero_Locus


END MODULE vdw_gas
