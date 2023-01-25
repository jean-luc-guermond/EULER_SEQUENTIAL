MODULE euler_start
  USE euler_bc_arrays
  USE euler_boundary_conditions
  USE boundary_conditions
  USE pardiso_solve
  USE update_euler
  USE space_dim
  USE mesh_handling
  USE input_data
CONTAINS
  SUBROUTINE start_euler
    inputs%syst_size = k_dim+2
    CALL read_my_data('data')
    CALL construct_mesh
    CALL construct_euler_bc
    CALL construct_euler_matrices
    !===Pardiso parameters
    isolve_euler_pardiso = -1 !===
    CALL allocate_pardiso_parameters(1)
    pardiso_param(1)%mtype = 1    !===real and structurally symmetric matrix
    pardiso_param(1)%phase = 33   !===Direct solve
    pardiso_param(1)%parm(4) = 0  !===Direct solve
  END SUBROUTINE start_euler
END MODULE euler_start
