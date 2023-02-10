PROGRAM compressible_euler
  USE mesh_handling
  USE sub_plot
  USE fem_tn
  USE timing_tools
  USE character_strings
  USE regression_test
  USE input_data
  USE IDP_euler_start
  USE IDP_update_euler
  USE IDP_euler_boundary_conditions
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un
  INTEGER                                   :: it, it_max, it_plot, it_time
  REAL(KIND=8)                              :: tps, to, norm, t_plot
  CHARACTER(len=1)                          :: tit
  CALL IDP_start_euler
  ALLOCATE(un(mesh%np,k_dim+2))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(:,1), 'rhoinit.plt')
  CALL IDP_COMPUTE_DT(un)

  !===Plotting
  it_plot = 0
  t_plot = inputs%dt_plot
  !===End plotting

  it_time = 0
  tps = user_time()
  it_max = INT(inputs%Tfinal/inputs%dt)
  inputs%dt = inputs%Tfinal/it_max
  it = 0
  DO WHILE(inputs%time<inputs%Tfinal)
     it = it + 1
     IF (inputs%if_regression_test .AND. it>5) THEN
           EXIT
     END IF
     CALL IDP_COMPUTE_DT(un)
     IF (inputs%time + inputs%dt>inputs%Tfinal) THEN
        inputs%dt=inputs%Tfinal-inputs%time
     END IF
      !write(*,*) 'time ', time, inputs%dt
     CALL full_step_ERK(un) 
     inputs%time = inputs%time + inputs%dt
     WRITE(*,*) 'time, dt ', inputs%time, inputs%dt
     IF ((inputs%time<t_plot .AND. inputs%time+inputs%dt/2>t_plot) .OR. &
          (inputs%time.GE.t_plot .AND. inputs%time-inputs%dt/2<t_plot)) THEN
        it_plot = it_plot + 1
        WRITE(tit,'(i1)') it_plot
        t_plot = t_plot+inputs%dt_plot
        CALL plot_scalar_field(mesh%jj, mesh%rr, un(:,1), 'rho_'//tit//'.plt')
     END IF
     it_time = it_time + 1
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step per grid point', tps/(it_time*mesh%np), 'itmax', it_time

  !===Regression test
  CALL regression(un)
     
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(:,1), 'rho.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(:,2)/un(:,1), 'ux.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(:,3)/un(:,1), 'uy.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, pressure(un), 'press.plt')
  CALL ns_0 (mesh, un(:,1),  norm)
  WRITE(*,*) ' L2 norm of density ', norm
CONTAINS
END PROGRAM compressible_euler
