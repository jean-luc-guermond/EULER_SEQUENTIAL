PROGRAM euler_poisson
  USE input_data
  USE mesh_handling
  USE euler_poisson_time_stepping
  USE boundary_conditions
  USE sub_plot
  USE lambda_module
  USE fem_tn
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, uo, ui
  INTEGER                                   :: i, it=0, it_max
  REAL(KIND=8)                              :: tps, to, norm

  CALL read_my_data('data')
  CALL construct_mesh
  CALL construct_matrices
  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_1d(mesh%rr(1,:),un(1,:),'rhoinit.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:)/un(1,:),'uinit.plt')
  CALL plot_1d(mesh%rr(1,:),un(3,:),'Einit.plt')
  CALL COMPUTE_DT(un)
  !TEST
  !inputs%dt=1.25d-3
  !TEST
  tps = user_time()
  it_max = int(inputs%Tfinal/inputs%dt)
  inputs%dt = inputs%Tfinal/it_max
  tps = user_time() 
  !DO i = 1, it_max
  DO WHILE(inputs%time<inputs%Tfinal)
     it = it + 1 
     CALL COMPUTE_DT(un)
     IF (inputs%time + inputs%dt>inputs%Tfinal) THEN
        inputs%dt=inputs%Tfinal-inputs%time
     END IF
     !write(*,*) 'time, dt ', inputs%time, inputs%dt
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
     !it = it + 1
     !WRITE(*,*) 'mass',sum(un(1,:)*lumped),sum(un(2,:)*lumped),sum(un(3,:)*lumped)
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step and dof', tps/(it*mesh%np)
  SELECT CASE(k_dim)
  CASE(1)
     CALL plot_1d(mesh%rr(1,:),un(1,:),'rho.plt')
     CALL plot_1d(mesh%rr(1,:),1.d0/un(1,:),'v.plt')
     CALL plot_1d(mesh%rr(1,:),un(2,:),'mt.plt')
     CALL plot_1d(mesh%rr(1,:),un(3,:),'E.plt')
     CALL plot_1d(mesh%rr(1,:),pressure(un), 'press.plt')
     CALL plot_1d(mesh%rr(1,:),un(2,:)/un(1,:),'u.plt')
     ui(1,:) = (un(3,:) - un(2,:)**2/(2.d0*un(1,:)))/un(1,:)**gamma
     CALL plot_1d(mesh%rr(1,:),ui(1,:),'es.plt')
     IF (inputs%type_test=='blas' .OR. inputs%type_test=='RPE4') STOP
     CALL comput_errors(un,inputs%time)
  CASE(2)
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uy.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, pressure(un), 'press.plt')
     CALL ns_0 (mesh, un(1,:),  norm)
     WRITE(*,*) ' L2 norm of density ', norm
  END SELECT
CONTAINS

  SUBROUTINE BC(unext,time)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: unext
    REAL(KIND=8)                             :: time
    INTEGER :: k
    unext(1,rho_js_D) = sol_anal(1,mesh%rr(:,rho_js_D))
    DO k = 1, k_dim !k_dim=1 here
       unext(1+k,ux_js_D)  = sol_anal(1+k,mesh%rr(:,ux_js_D))
    END DO
    unext(k_dim+2,rho_js_D) = sol_anal(k_dim+2,mesh%rr(:,rho_js_D))
  END SUBROUTINE BC



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
    WRITE(unit,*) '%mt=4'
    WRITE(unit,*) '%mc=2'
    WRITE(unit,*) '%lc=2'
    DO n = 1, SIZE(rr)
       WRITE(unit,*) rr(n), un(n)
    END DO
    CLOSE(unit)
  END SUBROUTINE plot_1d

  SUBROUTINE comput_errors(un,time)
    USE mesh_handling
    USE fem_tn
    USE shock_tub
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)     :: un
    REAL(KIND=8),                     INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(:), POINTER :: x_e, v_e, r_e, p_e
    REAL(KIND=8), DIMENSION(mesh%np) :: m_e, E_e, s_e, err_s
    REAL(KIND=8) :: error_rho_1, error_rho_2, error_rho_inf, error_mt_1, error_mt_2, error_mt_inf, &
         error_E_1, error_E_2, error_E_inf
    REAL(KIND=8) :: norm_rho_1, norm_rho_2, norm_rho_inf, norm_mt_1, norm_mt_2, norm_mt_inf, &
         norm_E_1, norm_E_2, norm_E_inf
    INTEGER :: itest

    SELECT CASE(inputs%type_test)
    CASE ('laxt') ! Lax
       itest = 2
    CASE ('sodt') ! Sod
       itest = 1
    CASE ('lebl') ! Leblanc
       itest = 4
    CASE ('soni') ! Sonic
       itest = 5
    CASE ('expa') ! Expansion
       itest = 6
    CASE ('smth') ! smooth
       itest = 7  
    CASE DEFAULT
       WRITE(*,*) ' BUG '
       STOP
    END SELECT
    SELECT CASE(inputs%type_test)
    CASE ('expa','lebl', 'smth') 
       ALLOCATE(r_e(mesh%np),p_e(mesh%np))
       r_e = rho_anal(mesh%rr)
       m_e = mt_anal(1,mesh%rr)
       E_e = E_anal(mesh%rr)
       p_e = press_anal(mesh%rr)
    CASE DEFAULT
       CALL sol_exa_stub(MINVAL(mesh%rr(1,:)),MAXVAL(mesh%rr(1,:)),mesh%np,itest,time,x_e,v_e,r_e,p_e)
       m_e = r_e*v_e
       E_e = p_e/(gamma-1) + m_e**2/(2*r_e)
    END SELECT

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
    s_e = (E_e - m_e**2/(2*r_e))/r_e**gamma
    err_s = ABS(s_e -(un(k_dim+2,:)-un(2,:)**2/(2.d0*un(1,:)))/(un(1,:)**gamma))
    WRITE(*,*) ' Error log(s), Max-norm', MAXVAL(err_s)/MAXVAL(s_e)
  END SUBROUTINE comput_errors

END PROGRAM
