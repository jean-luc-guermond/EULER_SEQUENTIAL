PROGRAM compressible_euler
  USE input_data
  USE mesh_handling
  USE sub_plot
  USE fem_tn
  USE timing_tools
  USE regression_test
  USE euler_start
  USE update_euler
  USE character_strings
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, uo, ui
  INTEGER                                   :: i, j, it=0, it_max, it_time=0
  REAL(KIND=8)                              :: tps, to, norm
  CALL start_euler
  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL COMPUTE_DT(un)
  tps = user_time()
  it_max = int(inputs%Tfinal/inputs%dt)
  inputs%dt = inputs%Tfinal/it_max
  tps = user_time() 
  DO WHILE(inputs%time<inputs%Tfinal)
     it = it + 1
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
     CALL euler_BC(un,to+inputs%dt) !t+dt/2
     inputs%time = to + inputs%dt
     it_time = it_time + 1
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step and dof', tps/(it_time*mesh%np)
  
  !===Regression test
  CALL regression(un)
  
  CALL plot_1d(mesh%rr(1,:),un(1,:),'rho.plt')
  CALL plot_1d(mesh%rr(1,:),pressure(un), 'press.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:)/un(1,:),'u.plt')

END PROGRAM compressible_euler
