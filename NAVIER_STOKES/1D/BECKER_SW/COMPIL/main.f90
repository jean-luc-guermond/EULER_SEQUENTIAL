PROGRAM compressible_ns
  USE input_data
  USE mesh_handling
  USE strang_splitting_compressible_ns
  USE sub_plot
  USE timing_tools
  USE navier_stokes_start
  USE boundary_conditions
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, unext
  REAL(KIND=8) :: tps, tinit_checkpoint, t_plot, dt_plot
  INTEGER :: checkpoint_trigger=1, it_time=0, it_plot, it, j
  CHARACTER(len=1) :: tit
  CALL start_navier_stokes
  ALLOCATE(un(inputs%syst_size,mesh%np),unext(inputs%syst_size,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr,inputs%time)
  CALL plot_1d(mesh%rr(1,:),un(1,:),'rhoinit.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:),'mtinit.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:)/un(1,:),'velinit.plt')
  CALL plot_1d(mesh%rr(1,:),un(3,:),'Einit.plt')
  !===Plots
  it_plot = 0
  t_plot = inputs%time+dt_plot
  !===Checkpointing
  tinit_checkpoint=inputs%time
  !===First dt to get started
  CALL COMPUTE_DT(un)
  tps = user_time()
  !DO it = 1, 2
  DO WHILE(inputs%time.LE.inputs%Tfinal)
     !write(*,*) ' time', inputs%time, inputs%dt, (sum(lumped*un(j,:)), j = 1, inputs%syst_size) 
     CALL time_stepping_compressible_ns(un,unext)
     un = unext
     inputs%time = inputs%time + inputs%dt
     CALL COMPUTE_DT(un)
     IF (inputs%time + inputs%dt>inputs%Tfinal) THEN
        inputs%dt = inputs%Tfinal - inputs%time + 1.d-15
     END IF
     it_time = it_time+1

     !===Check point
     IF (FLOOR((inputs%time-tinit_checkpoint)/inputs%checkpointing_freq)==checkpoint_trigger &
          .OR. inputs%Tfinal.LE.inputs%time + inputs%dt/100) THEN
        CALL write_restart(un)
        checkpoint_trigger = checkpoint_trigger + 1  
     END IF
     !===Plots
     IF ((inputs%time<t_plot .AND. inputs%time+inputs%dt/2>t_plot) .OR. &
          (inputs%time.GE.t_plot .AND. inputs%time-inputs%dt/2<t_plot)) THEN
        it_plot = it_plot + 1
        WRITE(tit,'(i1)') it_plot
        t_plot = t_plot+dt_plot
        CALL plot_1d(mesh%rr(1,:),un(1,:),'rho_'//tit//'.plt')
     END IF
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step per dof', tps/(it_time*mesh%np), 'itmax', it_time
  CALL plot_1d(mesh%rr(1,:),un(1,:),'rho.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:),'mt.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:)/un(1,:),'vel.plt')
  CALL plot_1d(mesh%rr(1,:),un(3,:),'E.plt')
  CALL comput_errors(un,inputs%time)
CONTAINS
  SUBROUTINE write_restart(un)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    OPEN(unit = 10, file = 'restart.'//inputs%file_name, form = 'unformatted', status = 'unknown')
    WRITE(10) inputs%time, un
    WRITE(*,*) ' inputs%time at checkpoint', inputs%time 
    CLOSE(10)
  END SUBROUTINE write_restart

  SUBROUTINE comput_errors(un,time)
    USE mesh_handling
    USE fem_tn
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)     :: un
    REAL(KIND=8),                     INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(mesh%np) :: m_e, E_e, r_e
    REAL(KIND=8) :: error_rho_1, error_rho_2, error_rho_inf, error_mt_1, error_mt_2, error_mt_inf, &
         error_E_1, error_E_2, error_E_inf
    REAL(KIND=8) :: norm_rho_1, norm_rho_2, norm_rho_inf, norm_mt_1, norm_mt_2, norm_mt_inf, &
         norm_E_1, norm_E_2, norm_E_inf
    INTEGER :: itest
    r_e = rho_anal(mesh%rr,time)
    m_e = mt_anal(1,mesh%rr,time)
    E_e = E_anal(mesh%rr,time)

    CALL ns_l1 (mesh, un(1,:) - r_e, error_rho_1)
    CALL ns_l1 (mesh,           r_e, norm_rho_1)
    CALL ns_0  (mesh, un(1,:) - r_e, error_rho_2)
    CALL ns_0  (mesh,           r_e, norm_rho_2)
    CALL ns_l1 (mesh, un(2,:) - m_e, error_mt_1)
    CALL ns_l1 (mesh,           m_e, norm_mt_1)
    CALL ns_0  (mesh, un(2,:) - m_e, error_mt_2)
    CALL ns_0  (mesh,           m_e, norm_mt_2)
    CALL ns_l1 (mesh, un(3,:) - E_e, error_E_1)
    CALL ns_l1 (mesh,           E_e, norm_E_1)
    CALL ns_0  (mesh, un(3,:) - E_e, error_E_2)
    CALL ns_0  (mesh,           E_e, norm_E_2)
    WRITE(*,*) ' Global relative error L1-norm    ', error_rho_1/norm_rho_1+error_mt_1/norm_mt_1+error_E_1/norm_E_1
    WRITE(*,*) ' Global relative error L2-norm    ', error_rho_2/norm_rho_2+error_mt_2/norm_mt_2+error_E_2/norm_E_2
    error_rho_inf = MAXVAL(ABS(un(1,:) - r_e))
    norm_rho_inf  = MAXVAL(ABS(r_e))
    error_mt_inf  = MAXVAL(ABS(un(2,:) - m_e))
    norm_mt_inf   = MAXVAL(ABS(m_e))
    error_E_inf   = MAXVAL(ABS(un(3,:) - E_e))
    norm_E_inf    = MAXVAL(ABS(E_e))
    WRITE(*,*) ' Global relative error Linfty-norm',&
         error_rho_inf/norm_rho_inf+error_mt_inf/norm_mt_inf+error_E_inf/norm_E_inf

    WRITE(*,*) ' Error density relative, L1-norm    ', error_rho_1/norm_rho_1
    WRITE(*,*) ' Error density relative, L2-norm    ', error_rho_2/norm_rho_2
    WRITE(*,*) ' Error density relative, Linfty-norm', error_rho_inf/norm_rho_inf

    CALL plot_1d(mesh%rr(1,:),r_e,'rho_exact')
    CALL plot_1d(mesh%rr(1,:),un(1,:)-r_e,'err_rho')
    CALL plot_1d(mesh%rr(1,:),m_e,'mt_exact')
    CALL plot_1d(mesh%rr(1,:),un(2,:)-m_e,'err_mt')
    CALL plot_1d(mesh%rr(1,:),E_e,'E_exact')
    CALL plot_1d(mesh%rr(1,:),un(3,:)-E_e,'err_E')
 

  END SUBROUTINE comput_errors
  
END PROGRAM compressible_ns
