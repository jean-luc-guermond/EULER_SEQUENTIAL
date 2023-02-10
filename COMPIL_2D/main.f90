PROGRAM compressible_euler
  USE input_data
  USE mesh_handling
  USE sub_plot
  USE fem_tn
  USE timing_tools
  USE euler_start
  USE update_euler
  USE character_strings
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, uo, ui
  INTEGER                                   :: it, it_max, it_plot, it_time
  REAL(KIND=8)                              :: tps, to, norm, t_plot, dt_plot
  CHARACTER(len=1)                          :: tit
  CALL start_euler
  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rhoinit.plt')
  CALL COMPUTE_DT(un)

  !===Plotting
  OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
  CALL read_until(21, '===dt_plot')
  READ (21,*) dt_plot
  it_plot = 0
  t_plot = dt_plot
  CLOSE(21)
  !===End plotting

  it_time = 0
  tps = user_time()
  it_max = INT(inputs%Tfinal/inputs%dt)
  inputs%dt = inputs%Tfinal/it_max
  !DO it = 1, it_max
  DO WHILE(inputs%time<inputs%Tfinal)
     CALL COMPUTE_DT(un)
     IF (inputs%time + inputs%dt>inputs%Tfinal) THEN
        inputs%dt=inputs%Tfinal-inputs%time
     END IF
     to = inputs%time
     uo = un 
     !===Step 1
     CALL euler(uo,un)
     CALL euler_BC(un,to+inputs%dt) !t+dt
     !===Step 2
     inputs%time=to+inputs%dt
     CALL euler(un,ui)
     un = (3*uo+ ui)/4
     CALL euler_BC(un,to+inputs%dt/2) !t+dt/2
     !===Step 3
     inputs%time =  to + inputs%dt/2
     CALL euler(un,ui)
     un = (uo+ 2*ui)/3
     CALL euler_BC(un,to+inputs%dt) !t+dt
     inputs%time = to + inputs%dt
     WRITE(*,*) 'time, dt ', inputs%time, inputs%dt
     IF ((inputs%time<t_plot .AND. inputs%time+inputs%dt/2>t_plot) .OR. &
          (inputs%time.GE.t_plot .AND. inputs%time-inputs%dt/2<t_plot)) THEN
        it_plot = it_plot + 1
        WRITE(tit,'(i1)') it_plot
        t_plot = t_plot+dt_plot
        CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho_'//tit//'.plt')
     END IF
     it_time = it_time + 1
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step per grid point', tps/(it_time*mesh%np), 'itmax', it_time
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uy.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, pressure(un), 'press.plt')
  CALL ns_0 (mesh, un(1,:),  norm)
  WRITE(*,*) ' L2 norm of density ', norm
CONTAINS
END PROGRAM compressible_euler
