MODULE  ideal_polytropic_gas
  USE def_of_gamma
  !REAL(KIND=8), PUBLIC, PARAMETER  ::  gamma_ = 1.4
  REAL(KIND=8), PUBLIC  :: gamma_m1_,  gamma_p1_, alpha_, beta_, tau_ 
CONTAINS
  SUBROUTINE compute_gamma_stuff
    IMPLICIT NONE
    gamma_m1_ =  gamma_ - 1              
    gamma_p1_ =  gamma_ + 1             
    alpha_ = (gamma_ - 1)/2          
    beta_ = (gamma_ + 1)/(gamma_ - 1)
    tau_ = (gamma_ - 1)/(2*gamma_)
  END SUBROUTINE compute_gamma_stuff
END MODULE  ideal_polytropic_gas
!!$MODULE  ideal_polytropic_gas
!!$
!!$!   PRIVATE
!!$
!!$   REAL(KIND=8), PUBLIC, PARAMETER  ::  gamma_ = 1.4
!!$
!!$   REAL(KIND=8), PUBLIC, PARAMETER  ::  gamma_m1_ =  gamma_ - 1,              &
!!$                                        gamma_p1_ =  gamma_ + 1,              &
!!$                                           alpha_ = (gamma_ - 1)/2,           &
!!$                                            beta_ = (gamma_ + 1)/(gamma_ - 1),&
!!$                                             tau_ = (gamma_ - 1)/(2*gamma_)
!!$
!!$END MODULE  ideal_polytropic_gas
