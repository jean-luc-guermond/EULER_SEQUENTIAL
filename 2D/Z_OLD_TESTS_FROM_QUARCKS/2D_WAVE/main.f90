PROGRAM scal_cons

  USE input_data
  USE mesh_handling
  USE update
  USE boundary_conditions
  USE sub_plot
  USE lambda_module
  USE fem_tn
  USE character_strings
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, uo, ui
  INTEGER                                   :: it, it_max, it_plot
  REAL(KIND=8)                              :: tps, to, norm, t_plot, dt_plot
  CHARACTER(len=1)                          :: tit
  CALL read_my_data('data')
  CALL construct_mesh
  CALL construct_matrices
  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rhoinit.plt')
  !CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux.plt')
  !CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uy.plt')
  !CALL plot_scalar_field(mesh%jj, mesh%rr, un(4,:), 'E.plt')
  CALL COMPUTE_DT(un)

  !===Plotting
  OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
  CALL read_until(21, '===dt_plot')
  READ (21,*) dt_plot
  it_plot = 0
  t_plot = dt_plot
  CLOSE(21)
  !===End plotting
  
  tps = user_time()
  it_max = int(inputs%Tfinal/inputs%dt)
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
     CALL BC(un,to+inputs%dt) !t+dt
     !===Step 2
     inputs%time=to+inputs%dt
     CALL euler(un,ui)
     un = (3*uo+ ui)/4
     CALL BC(un,to+inputs%dt/2) !t+dt/2
     !===Step 3
     inputs%time =  to + inputs%dt/2
     CALL euler(un,ui)
     un = (uo+ 2*ui)/3
     CALL BC(un,to+inputs%dt) !t+dt/2
     inputs%time = to + inputs%dt
     !write(*,*) 'time, dt ', inputs%time, inputs%dt
     
     IF ((inputs%time<t_plot .AND. inputs%time+inputs%dt/2>t_plot) .OR. &
          (inputs%time.GE.t_plot .AND. inputs%time-inputs%dt/2<t_plot)) THEN
        it_plot = it_plot + 1
        WRITE(tit,'(i1)') it_plot
        t_plot = t_plot+dt_plot
        CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho_'//tit//'.plt')
     END IF
     
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step', tps/it, 'itmax', it
  SELECT CASE(k_dim)
  CASE(1)
     CALL plot_1d(mesh%rr(1,:),un(1,:),'rho.plt')
  CASE(2)
     CALL comput_errors(un,inputs%time)
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uy.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, pressure(un), 'E.plt')
     CALL ns_0 (mesh, un(1,:),  norm)
     WRITE(*,*) ' L2 norm of density ', norm
  END SELECT
CONTAINS

  SUBROUTINE BC(unext,time)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: unext
    REAL(KIND=8)                             :: time
    unext(1,rho_js_D) = sol_anal(1,mesh%rr(:,rho_js_D))
    unext(2,ux_js_D)  = sol_anal(2,mesh%rr(:,ux_js_D))
    unext(3,uy_js_D)  = sol_anal(3,mesh%rr(:,uy_js_D))
    unext(4,rho_js_D)  = sol_anal(4,mesh%rr(:,rho_js_D))
  END SUBROUTINE BC

  SUBROUTINE COMPUTE_DT(u0)
    USE input_data
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: u0
    CALL compute_dij(u0)
    inputs%dt = inputs%CFL*1/MAXVAL(ABS(dij%aa(diag))/lumped)
  END SUBROUTINE COMPUTE_DT

  FUNCTION user_time() RESULT(time)
    IMPLICIT NONE
    REAL(KIND=8) :: time 
    INTEGER :: count, count_rate, count_max
    CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
    time = (1.d0*count)/count_rate
  END FUNCTION user_time

  SUBROUTINE plot_1d(rr,un,file)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rr, un
    INTEGER :: n, unit=10
    CHARACTER(*) :: file
    OPEN(unit,FILE=TRIM(ADJUSTL(file)),FORM='formatted')
    WRITE(unit,*) '%toplabel='' '''
    WRITE(unit,*) '%xlabel='' '''
    WRITE(unit,*) '%ylabel='' '''
    WRITE(unit,*) '%ymax=', MAXVAL(un)
    WRITE(unit,*) '%ymin=', MINVAL(un)
    WRITE(unit,*) '%xyratio=1'
    DO n = 1, SIZE(rr)
       WRITE(unit,*) rr(n), un(n)
    END DO
    CLOSE(unit)
  END SUBROUTINE plot_1d

 SUBROUTINE comput_errors(un,time)
    USE mesh_handling
    USE fem_tn
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)     :: un
    REAL(KIND=8),                     INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(:), POINTER :: x_e, v_e, r_e, p_e
    REAL(KIND=8), DIMENSION(mesh%np) :: m_e1, m_e2, E_e, s_e, err_s
    REAL(KIND=8) :: error_rho_1, error_rho_2, error_rho_inf, error_mt_1, error_mt_2, error_mt_inf, &
         error_E_1, error_E_2, error_E_inf
    REAL(KIND=8) :: norm_rho_1, norm_rho_2, norm_rho_inf, norm_mt_1, norm_mt_2, norm_mt_inf, &
         norm_E_1, norm_E_2, norm_E_inf, em1, em2, nm1, nm2
    INTEGER :: itest

 
    ALLOCATE(r_e(mesh%np),p_e(mesh%np))
    r_e = rho_anal(mesh%rr)
    m_e1 = mt_anal(1,mesh%rr)
    m_e2 = mt_anal(2,mesh%rr)
    E_e = E_anal(mesh%rr)
    p_e = press_anal(mesh%rr)
    
    CALL ns_l1 (mesh, un(1,:) - r_e, error_rho_1)
    CALL ns_l1 (mesh,           r_e, norm_rho_1)
    CALL ns_0  (mesh, un(1,:) - r_e, error_rho_2)
    CALL ns_0  (mesh,           r_e, norm_rho_2)
    
    CALL ns_l1 (mesh, un(2,:) - m_e1, error_mt_1)
    CALL ns_l1 (mesh,           m_e1, norm_mt_1)
    CALL ns_0  (mesh, un(2,:) - m_e1, error_mt_2)
    CALL ns_0  (mesh,           m_e1, norm_mt_2)
    CALL ns_l1 (mesh, un(3,:) - m_e2, em1)
    CALL ns_l1 (mesh,           m_e2, nm1)
    CALL ns_0  (mesh, un(3,:) - m_e2, em2)
    CALL ns_0  (mesh,           m_e2, nm2)
    error_mt_1 =  error_mt_1 + em1
    error_mt_2 =  SQRT(error_mt_2**2 + em2**2)
    norm_mt_1 =  norm_mt_1 + em1
    norm_mt_2 =  SQRT(norm_mt_2**2 + em2**2)
    
    CALL ns_l1 (mesh, un(4,:) - E_e, error_E_1)
    CALL ns_l1 (mesh,           E_e, norm_E_1)
    CALL ns_0  (mesh, un(4,:) - E_e, error_E_2)
    CALL ns_0  (mesh,           E_e, norm_E_2)
    WRITE(*,*) ' Global relative error L1-norm    ', error_rho_1/norm_rho_1+error_mt_1/norm_mt_1+error_E_1/norm_E_1
    WRITE(*,*) ' Global relative error L2-norm    ', error_rho_2/norm_rho_2+error_mt_2/norm_mt_2+error_E_2/norm_E_2
    error_rho_inf = MAXVAL(ABS(un(1,:) - r_e))
    norm_rho_inf  = MAXVAL(ABS(r_e))
    
    error_mt_inf  = MAXVAL(ABS(un(2,:) - m_e1))
    norm_mt_inf   = MAXVAL(ABS(m_e1))
    error_mt_inf  = MAX(error_mt_inf,MAXVAL(ABS(un(3,:) - m_e2)))
    norm_mt_inf   = MAX(norm_mt_inf,MAXVAL(ABS(m_e2)))

    
    error_E_inf   = MAXVAL(ABS(un(4,:) - E_e))
    norm_E_inf    = MAXVAL(ABS(E_e))
    WRITE(*,*) ' Global relative error Linfty-norm',&
         error_rho_inf/norm_rho_inf+error_mt_inf/norm_mt_inf+error_E_inf/norm_E_inf
    
    WRITE(*,*) ' Error density relative, L1-norm    ', error_rho_1/norm_rho_1
    WRITE(*,*) ' Error density relative, L2-norm    ', error_rho_2/norm_rho_2
    WRITE(*,*) ' Error density relative, Linfty-norm', error_rho_inf/norm_rho_inf


    CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:)-r_e, 'err_rho.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, (un(4,:)-E_e)/E_e, 'err_E.plt')
    s_e = (E_e - (m_e1**2+m_e2**2)/(2*r_e))/r_e**gamma
    err_s = ABS(s_e -(un(k_dim+2,:)-un(2,:)**2/(2.d0*un(1,:)))/(un(1,:)**gamma))
    WRITE(*,*) ' Error log(s), Max-norm', MAXVAL(err_s)/MAXVAL(s_e)
  END SUBROUTINE comput_errors


  
END PROGRAM scal_cons
