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
  REAL(KIND=8)                              :: tps, to, norm, t_plot
  REAL(KIND=8)                              :: tinit_checkpoint
  INTEGER                                   :: checkpoint_trigger=1
  CHARACTER(len=2)                          :: tit
  CALL start_euler
  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rhoinit.plt')

  !===Checkpointing
  tinit_checkpoint=inputs%time

  CALL COMPUTE_DT(un)

  !===Plotting
  it_plot = 0
  t_plot = inputs%time+inputs%dt_plot

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
     CALL euler_BC(un,to+inputs%dt) !t+dt/2
     inputs%time = to + inputs%dt
     WRITE(*,*) 'time, dt ', inputs%time, inputs%dt
     !===Check point
     IF (FLOOR((inputs%time-tinit_checkpoint)/inputs%checkpointing_freq)==checkpoint_trigger &
          .OR. inputs%Tfinal.LE.inputs%time + inputs%dt/100) THEN
        WRITE(tit,'(i2)') checkpoint_trigger
        CALL write_restart(un,tit)
        checkpoint_trigger = checkpoint_trigger + 1  
     END IF
     !===Plots
     IF ((inputs%time<t_plot .AND. inputs%time+inputs%dt/2>t_plot) .OR. &
          (inputs%time.GE.t_plot .AND. inputs%time-inputs%dt/2<t_plot)) THEN
        it_plot = it_plot + 1
        WRITE(tit,'(i2)') it_plot
        t_plot = t_plot+inputs%dt_plot
        CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho_'//trim(adjustl(tit))//'.plt')
        CALL schlieren(un(1,:),uo(1,:))
        CALL vtk_2d(mesh, uo(1,:), 10, 'grad_'//trim(adjustl(tit))//'.vtk','scalars')
        CALL plot_scalar_field(mesh%jj, mesh%rr, uo(1,:), 'grad_'//trim(adjustl(tit))//'.plt')
        CALL vtk_2d(mesh, un(1,:), 10, 'density_'//trim(adjustl(tit))//'.vtk','scalars')
     END IF
     it_time = it_time + 1
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step per grid point', tps/(it_time*mesh%np), 'itmax', it_time
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uy.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, pressure(un), 'press.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, pressure(un)*un(1,:)**(-gamma), 'entropy.plt')
  CALL ns_0 (mesh, un(1,:),  norm)
  WRITE(*,*) ' L2 norm of density ', norm
CONTAINS

   SUBROUTINE schlieren(un,grad) 
    implicit none
    REAL(KIND=8), DIMENSION(mesh%np)     :: un
    REAL(KIND=8), DIMENSION(mesh%np)     :: grad
    REAL(KIND=8), DIMENSION(k_dim) :: xx
    REAL(KIND=8) :: mx, mn
    integer :: i, j, k, p    
    !===Gradient
    DO i = 1, mesh%np
       xx = 0.d0
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = mass%ja(p)
          DO k = 1, k_dim
             xx(k) = xx(k) + un(j)*cij(k)%aa(p)
          END DO
       END DO
       grad(i) = SQRT(SUM(xx**2))/lumped(i)
    END DO
    mx=MAXVAL(grad)
    mn=MINVAL(grad)
    grad = EXP(-10*(grad-mn)/(mx-mn))
  END SUBROUTINE schlieren

  SUBROUTINE write_restart(un,tit)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    CHARACTER(len=2) :: tit
    OPEN(unit = 10, file = 'restart_'//trim(adjustl(tit))//'_'//inputs%file_name, form = 'unformatted', status = 'unknown')
    WRITE(10) inputs%time, un
    WRITE(*,*) ' inputs%time at checkpoint', inputs%time 
    CLOSE(10)
  END SUBROUTINE write_restart

END PROGRAM compressible_euler
