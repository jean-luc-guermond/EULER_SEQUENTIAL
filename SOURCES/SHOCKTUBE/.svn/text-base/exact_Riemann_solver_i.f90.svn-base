MODULE exact_Riemann_solver_i

CONTAINS

SUBROUTINE  exact_Riemann_solver (ul, pl, cl,  ur, pr, cr,   &
                                  uc, pc, clc, crc,          &
                                  rel_err, iterations,  l_r_w) ! OPTIONAL
!=======================================================================
!     Solution of the Riemann Problem for an Ideal Polytropic Gas
!     Iterative procedure of Gottlieb and Groth
!     Coding by Francesco Bassi
!     Cosmetical modifications by L.Q.  23 November 1996
!-----------------------------------------------------------------------
!
!     INPUT
!
!     ul,  pl,  cl  :  left state  (velocity, pressure, speed of sound)
!     ur,  pr,  cr  :  right state (velocity, pressure, speed of sound)
!
!     optional
!     rel_err       :  relative error (of pressure) to stop the iterations
!                      [default = 1.0d-7]
!     iterations    :  maximum number of iterations
!                      [default = 100]
!
!-----------------------------------------------------------------------
!
!     OUTPUT
!
!     uc,  pc,  clc,  crc     :  states across the contact discontinuity
!
!     optional
!     iterations              :  number of actually computed iterations
!     l_r_w(1) = 'r'  or  's' :  kind of the left wave  (raref or shock)
!     l_r_w(2) = 'r'  or  's' :  kind of the right wave (raref or shock)
!
!=======================================================================

   USE  ideal_polytropic_gas  ! this MODULE contains the value of gamma_

   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN)  :: ul, pl, cl,  ur, pr, cr
   REAL(KIND=8), INTENT(OUT) :: uc, pc, clc, crc

   REAL(KIND=8),            OPTIONAL, INTENT(IN)    :: rel_err
   INTEGER,                 OPTIONAL, INTENT(INOUT) :: iterations
   CHARACTER, DIMENSION(2), OPTIONAL, INTENT(OUT)   :: l_r_w


   REAL(KIND=8) :: gamma,gm,hgm,ghgm,gp,qgp,rocl,rocr,z,ult,urt,  &
                   dul,dur,aa,plc,prc,plcp,prcp, wl,wr,r

   INTEGER                 ::  max_it_Int = 100     ! default
   REAL(KIND=8)            :: rel_err_Int = 1.0d-07 ! default
   CHARACTER, DIMENSION(2) ::   l_r_w_Int
   INTEGER :: it

!-----------------------------------------------------------------------

   IF (PRESENT(iterations))   max_it_Int = iterations

   IF (PRESENT(rel_err))     rel_err_Int = rel_err

   gamma = gamma_  ! this value is taken from MODULE ideal_polytropic_gas

   gm = gamma - 1;  hgm = gm/2;  ghgm = gamma/hgm
   gp = gamma + 1;  qgp = gp/4

   rocl = gamma*pl/cl
   rocr = gamma*pr/cr

   ! initial guess for contact velocity

   z   = (pl/pr)**(1/ghgm) * cr/cl
   ult =  ul + cl/hgm
   urt =  ur - cr/hgm
   uc  = (ult*z + urt)/(1 + z)

   ! solve the Riemann problem

   DO it = 1, max_it_Int

      dul = uc - ul

      IF (dul <= 0) THEN  ! leftward moving shock

         aa   =  qgp*dul/cl
         wl   =  aa - SQRT(1 + aa*aa)
         plc  =  pl + rocl*dul*wl
         plcp =  2*rocl * wl**3/(1 + wl*wl)
         clc  =  cl * SQRT((gp + gm*plc/pl)/(gp + gm*pl/plc))
         l_r_w_Int(1) = 's'

      ELSE  ! left rarefaction wave

         clc  =  cl - hgm*dul
         plc  =  pl * (clc/cl)**ghgm
         plcp = -gamma*plc/clc
         l_r_w_Int(1) = 'r'

      END IF

      dur = uc - ur

      IF (dur >= 0) THEN  ! rightward moving shock

         aa   =  qgp*dur/cr
         wr   =  aa + SQRT(1 + aa*aa)
         prc  =  pr + rocr*dur*wr
         prcp =  2*rocr * wr**3/(1 + wr*wr)
         crc  =  cr * SQRT((gp + gm*prc/pr)/(gp + gm*pr/prc))
         l_r_w_Int(2) = 's'

      ELSE  ! right rarefaction wave

         crc  =  cr + hgm*dur
         prc  =  pr * (crc/cr)**ghgm
         prcp =  gamma*prc/crc
         l_r_w_Int(2) = 'r'

      END IF

      ! convergence check

      IF (ABS(1 - plc/prc) <= rel_err_Int) THEN
         pc = plc
         IF (PRESENT(iterations))  iterations = it
         IF (PRESENT(l_r_w))  l_r_w = l_r_w_Int
         RETURN
      ENDIF

      uc = uc - (plc - prc)/(plcp - prcp)

   ENDDO

   PRINT*,' Riemann solver fails to converge'

END SUBROUTINE  exact_Riemann_solver


END MODULE  exact_Riemann_solver_i
