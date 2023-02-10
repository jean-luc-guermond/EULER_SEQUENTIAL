PROGRAM compressible_euler
  USE input_data
  USE mesh_handling
  USE sub_plot
  USE fem_tn
  USE timing_tools
  USE euler_start
  USE update_euler
  USE character_strings
  USE regression_test
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, uo, ui
  INTEGER                                   :: it, it_reg=0, it_max, it_plot, it_time
  REAL(KIND=8)                              :: tps, to, norm, t_plot, curl_max, dp
  REAL(KIND=8)                              :: tinit_checkpoint
  INTEGER                                   :: checkpoint_trigger=1
  CHARACTER(len=2)                          :: tit
  CALL start_euler
  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rhoinit.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:), 'uxinit.plt')

  !===Curl of initial velocity
  uo(2,:) = un(2,:)/un(1,:)
  uo(3,:) = un(3,:)/un(1,:)
  CALL vorticity(uo(2:3,:),uo(1,:))
  curl_max = MAXVAL(ABS(uo(1,:)))
  uo(1,:) = pressure(un)
  dp = MAXVAL(uo(1,:)) -  MINVAL(uo(1,:))

  !===Checkpointing
  tinit_checkpoint=inputs%time

  CALL COMPUTE_DT(un)

  !===Plotting
  it_plot = 0
  t_plot = inputs%time  !+inputs%dt_plot

  it_time = 0
  tps = user_time()
  it_max = INT(inputs%Tfinal/inputs%dt)
  inputs%dt = inputs%Tfinal/it_max
  !DO it = 1, it_max
  DO WHILE(inputs%time<inputs%Tfinal)
     it_reg = it_reg + 1
     IF (inputs%if_regression_test .AND. it_reg>5) THEN
        EXIT
     END IF
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
     IF ((ABS(inputs%time-t_plot)<inputs%dt/2) .OR. it_plot==0) THEN
        WRITE(*,*) ' Creating vtk files at t=', inputs%time
        it_plot = it_plot + 1
        WRITE(tit,'(i2)') it_plot
        t_plot = t_plot+inputs%dt_plot
        CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho_'//TRIM(ADJUSTL(tit))//'.plt')
        CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux_'//TRIM(ADJUSTL(tit))//'.plt')
        uo(2,:) = un(2,:)/un(1,:)
        uo(3,:) = un(3,:)/un(1,:)
        CALL vorticity(uo(2:3,:),uo(1,:)) 
        CALL vtk_2d(mesh, uo(1,:), 10, 'curl_'//TRIM(ADJUSTL(tit))//'.vtk','scalars')
        CALL plot_scalar_field(mesh%jj, mesh%rr, uo(1,:), 'curl_'//TRIM(ADJUSTL(tit))//'.plt')
        CALL vtk_2d(mesh, un(1,:), 10, 'density_'//TRIM(ADJUSTL(tit))//'.vtk','scalars')
        CALL vtk_2d(mesh, un(2,:), 10, 'ux_'//TRIM(ADJUSTL(tit))//'.vtk','scalars')
        uo(1,:) = (pressure(un)-p_infty)/dp
        CALL plot_scalar_field(mesh%jj, mesh%rr, uo(1,:), 'press_'//TRIM(ADJUSTL(tit))//'.plt')
     END IF
     
     !===Regression test
     CALL regression(un)
     
     uo(2,:) = un(2,:)/un(1,:)
     uo(3,:) = un(3,:)/un(1,:)
     CALL vorticity(uo(2:3,:),uo(1,:))
     WRITE(30,*) inputs%time, MAXVAL(ABS(uo(1,:)))/curl_max
     WRITE(20,*) inputs%time, MAXVAL(SQRT((un(2,:)-u_infty)**2+un(3,:)**2)/un(1,:))/u_infty
     uo(1,:) = (pressure(un)-p_infty)/dp
     WRITE(40,*) inputs%time, MAXVAL(ABS(uo(1,:))), MAXVAL(uo(1,:)), MINVAL(uo(1,:)) 
     it_time = it_time + 1
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step per grid point', tps/(it_time*mesh%np), 'itmax', it_time
  CALL comput_errors(un,inputs%time)
  
CONTAINS

  SUBROUTINE schlieren(un,grad)
    USE lin_alg_tools
    USE matrix_type
    USE space_dim
    implicit none
    REAL(KIND=8), DIMENSION(mesh%np)     :: un
    REAL(KIND=8), DIMENSION(mesh%np)     :: grad
    TYPE(matrice_bloc), DIMENSION(k_dim), SAVE :: cij
    REAL(KIND=8), DIMENSION(k_dim) :: xx
    REAL(KIND=8) :: mx, mn
    integer :: i, j, k, p
    LOGICAL :: once=.true.
    IF (once) THEN
       CALL compute_cij(mesh,cij)
       once=.false.
    END IF
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

    SUBROUTINE vorticity(un,curl)
    USE lin_alg_tools
    USE matrix_type
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim,mesh%np)     :: un
    REAL(KIND=8), DIMENSION(mesh%np)     :: curl
    REAL(KIND=8) :: mx, mn, xx
    INTEGER :: i, j, k, p
    TYPE(matrice_bloc), DIMENSION(k_dim), SAVE :: cij
    LOGICAL :: once=.true.
    IF (once) THEN
       CALL compute_cij(mesh,cij)
       once=.false.
    END IF
    !===Gradient
    DO i = 1, mesh%np
       xx = 0.d0
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = mass%ja(p)
          xx = xx + un(2,j)*cij(1)%aa(p)-un(1,j)*cij(2)%aa(p)
       END DO
       curl(i) = xx/lumped(i)
    END DO
  END SUBROUTINE vorticity

  SUBROUTINE write_restart(un,tit)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    CHARACTER(len=2) :: tit
    OPEN(unit = 10, file = 'restart_'//TRIM(ADJUSTL(tit))//'_'//inputs%file_name, form = 'unformatted', status = 'unknown')
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
    REAL(KIND=8), DIMENSION(mesh%np)  :: mask
    REAL(KIND=8), DIMENSION(:), POINTER :: x_e, v_e, r_e, p_e
    REAL(KIND=8), DIMENSION(mesh%np) :: m_e1, m_e2, E_e, s_e, err_s
    REAL(KIND=8) :: error_rho_1, error_rho_2, error_rho_inf, error_mt_1, error_mt_2, error_mt_inf, &
         error_E_1, error_E_2, error_E_inf
    REAL(KIND=8) :: norm_rho_1, norm_rho_2, norm_rho_inf, norm_mt_1, norm_mt_2, norm_mt_inf, &
         norm_E_1, norm_E_2, norm_E_inf, em1, em2, nm1, nm2, r
    INTEGER :: itest, i


!!$    DO i = 1, mesh%np
!!$       r = SQRT((mesh%rr(1,i)-x0-u_infty*inputs%time)**2+(mesh%rr(2,i)-y0)**2)
!!$       IF  (r<2.5d0) THEN
!!$          mask(i) = 1.d0
!!$       ELSE IF (r<4.5d0) THEN
!!$          r =  (r-2.5d0)/(4.5d0-2.5d0)
!!$          mask(i) = (20*r**3 + 10*r**2 + 4*r + 1)*(r - 1)**4
!!$       ELSE  
!!$          mask(i) = 0.d0
!!$       END IF
!!$    END DO
    mask =1.d0
    !CALL plot_scalar_field(mesh%jj, mesh%rr, mask, 'mask.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uy.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, pressure(un), 'press.plt')
    !CALL ns_0 (mesh, un(1,:),  norm)
    !WRITE(*,*) ' L2 norm of density ', norm

    ALLOCATE(r_e(mesh%np),p_e(mesh%np))
    r_e = rho_anal(mesh%rr)
    m_e1 = mt_anal(1,mesh%rr)
    m_e2 = mt_anal(2,mesh%rr)
    E_e = E_anal(mesh%rr)
    p_e = press_anal(mesh%rr)

    CALL ns_l1 (mesh, mask*(un(1,:) - r_e), error_rho_1)
    CALL ns_l1 (mesh,           mask* r_e, norm_rho_1)
    CALL ns_0  (mesh, mask*(un(1,:) - r_e), error_rho_2)
    CALL ns_0  (mesh,            mask*r_e, norm_rho_2)

    CALL ns_l1 (mesh, mask*(un(2,:) - m_e1), error_mt_1)
    CALL ns_l1 (mesh,           mask*m_e1, norm_mt_1)
    CALL ns_0  (mesh, mask*(un(2,:) - m_e1), error_mt_2)
    CALL ns_0  (mesh,          mask* m_e1, norm_mt_2)
    CALL ns_l1 (mesh, mask*(un(3,:) - m_e2), em1)
    CALL ns_l1 (mesh,          mask* m_e2, nm1)
    CALL ns_0  (mesh, mask*(un(3,:) - m_e2), em2)
    CALL ns_0  (mesh,          mask* m_e2, nm2)
    error_mt_1 =  error_mt_1 + em1
    error_mt_2 =  SQRT(error_mt_2**2 + em2**2)
    norm_mt_1 =  norm_mt_1 + em1
    norm_mt_2 =  SQRT(norm_mt_2**2 + em2**2)
    CALL ns_l1 (mesh, mask*(un(4,:) - E_e), error_E_1)
    CALL ns_l1 (mesh,           mask*E_e, norm_E_1)
    CALL ns_0  (mesh, mask*(un(4,:) - E_e), error_E_2)
    CALL ns_0  (mesh,          mask*E_e, norm_E_2)
    WRITE(*,*) ' Global relative error L1-norm    ', error_rho_1/norm_rho_1+error_mt_1/norm_mt_1+error_E_1/norm_E_1
    WRITE(*,*) ' Global relative error L2-norm    ', error_rho_2/norm_rho_2+error_mt_2/norm_mt_2+error_E_2/norm_E_2

    error_rho_inf = MAXVAL(ABS(mask*(un(1,:) - r_e)))
    norm_rho_inf  = MAXVAL(ABS(r_e))

    error_mt_inf  = MAXVAL(ABS(mask*(un(2,:) - m_e1)))
    norm_mt_inf   = MAXVAL(ABS(m_e1))
    error_mt_inf  = MAX(error_mt_inf,MAXVAL(ABS(mask*(un(3,:) - m_e2))))
    norm_mt_inf   = MAX(norm_mt_inf,MAXVAL(ABS(m_e2)))

    error_E_inf   = MAXVAL(ABS(mask*(un(4,:) - E_e)))
    norm_E_inf    = MAXVAL(ABS(E_e))
    WRITE(*,*) ' Global relative error Linfty-norm',&
         error_rho_inf/norm_rho_inf+error_mt_inf/norm_mt_inf+error_E_inf/norm_E_inf

!!$    WRITE(*,*) ' Error density relative, L1-norm    ', error_rho_1/norm_rho_1
!!$    WRITE(*,*) ' Error density relative, L2-norm    ', error_rho_2/norm_rho_2
    WRITE(*,*) ' Error density relative, Linfty-norm', error_rho_inf/norm_rho_inf
!!$    WRITE(*,*) ' Err momentum relative, L1-norm     ', error_mt_1/norm_mt_1
!!$    WRITE(*,*) ' Err momentum relative, L2-norm     ', error_mt_2/norm_mt_2
    WRITE(*,*) ' Err momentum relative, Linfty-norm ', error_mt_inf/norm_mt_inf
!!$    WRITE(*,*) ' Error energy relative, L1-norm     ', error_E_1/norm_E_1
!!$    WRITE(*,*) ' Error energy relative, L2-norm     ', error_E_2/norm_E_2
    WRITE(*,*) ' Error energy relative, Linfty-norm ', error_E_inf/norm_E_inf

    CALL plot_scalar_field(mesh%jj, mesh%rr, (un(1,:)-r_e)/norm_rho_inf, 'err_rho.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, (un(2,:)-m_e1)/norm_mt_inf, 'err_ux.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, (un(3,:)-m_e2)/norm_mt_inf, 'err_uy.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, (un(4,:)-E_e)/norm_E_inf, 'err_E.plt')
    s_e = (E_e - (m_e1**2+m_e2**2)/(2*r_e))/r_e**gamma
    err_s = ABS(s_e -(un(k_dim+2,:)-un(2,:)**2/(2.d0*un(1,:)))/(un(1,:)**gamma))
    !WRITE(*,*) ' Error log(s), Max-norm', MAXVAL(err_s)/MAXVAL(s_e)
    CALL plot_scalar_field(mesh%jj, mesh%rr, r_e, 'rho_exact.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, m_e1/r_e, 'ux_exact.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, m_e2/r_e, 'uy_exact.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, p_e, 'p_exact.plt')
  END SUBROUTINE comput_errors


END PROGRAM compressible_euler
