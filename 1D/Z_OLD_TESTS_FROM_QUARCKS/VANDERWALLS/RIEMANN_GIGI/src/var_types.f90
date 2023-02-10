
! =========================================================================
!   Description: Contains variables types definitions
!
!   	 Author: Marco Fossati
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!
!          Year:  2005, October
!
!   MODIFICATIONS BY LQ to avoid POINTER 
!
! =========================================================================

MODULE var_types



   TYPE nonlin_wave

      INTEGER :: sign

      CHARACTER(LEN=1), DIMENSION(:),   ALLOCATABLE :: type
    ! CHARACTER(LEN=1), DIMENSION(:),   POINTER :: type

    ! REAL(KIND=8),     DIMENSION(:,:), POINTER :: v, u, p, &
      REAL(KIND=8),     DIMENSION(:,:), ALLOCATABLE :: v, u, p, &
                                                       speeds, sigma

   END TYPE nonlin_wave



   TYPE key_points

      INTEGER :: N_pts

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x
    ! REAL(KIND=8), DIMENSION(:), POINTER :: x
   
   END TYPE key_points



END MODULE var_types
