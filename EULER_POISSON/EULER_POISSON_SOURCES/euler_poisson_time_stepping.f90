MODULE euler_poisson_time_stepping
  USE update_euler_poisson
  USE update_euler
  IMPLICIT NONE
  INTEGER :: syst_size = k_dim+3
CONTAINS
  SUBROUTINE initialization
    IMPLICIT NONE
    CALL read_my_data('data')
    CALL construct_mesh
    CALL construct_matrices
    CALL construct_euler_poisson_matrices
  END SUBROUTINE initialization

  SUBROUTINE time_update(un,unext)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(syst_size,mesh%np),  INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(syst_size,,mesh%np), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(syst_size,mesh%np) :: ui, uo
    REAL(KIND=8) :: current_dt, to

    !===Compute time step
    CALL COMPUTE_DT(un)
    current_dt = inputs%dt

    !===First Strang half step using SSP2
    inputs%dt = current_dt/2
    to = inputs%time
    uo = un
    !===Step 1
    CALL euler(uo,un)
    CALL BC(un,to+inputs%dt) !t+dt
    !===Step 2
    CALL euler(un,ui)
    un = 0.5d0*(uo + ui)
    CALL BC(un,to+inputs%dt) !t+dt

    !===Second Strang step
    inputs%dt = current_dt
    CALL euler_poisson(un,ui)

    !===Second Strang half step using SSP2
    inputs%dt = current_dt/2
    to = inputs%time + inputs%dt !===to+dt/2 
    uo = ui
    !===Step 1
    CALL euler(uo,un)
    CALL BC(un,to+inputs%dt) !t+dt
    !===Step 2
    CALL euler(un,ui)
    unext = uo + ui
    CALL BC(un,to+inputs%dt) !t+dt

    !===Advance 
    inputs%dt = current_dt
    inputs%time = inputs%time + inputs%dt
    
  END SUBROUTINE time_update
  
END MODULE euler_poisson_time_stepping
