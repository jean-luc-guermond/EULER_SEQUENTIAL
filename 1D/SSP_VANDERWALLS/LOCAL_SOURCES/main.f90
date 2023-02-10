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
  USE vdw
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, uo, ui
  INTEGER                                   :: i, j, it=0, it_max, it_time=0
  REAL(KIND=8)                              :: tps, to, norm
  CALL start_euler
  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_1d(mesh%rr(1,:),pressure(un),'pinit.plt')
  CALL plot_1d(mesh%rr(1,:),un(1,:),'rhoinit.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:)/un(1,:),'uinit.plt')
  CALL plot_1d(mesh%rr(1,:),un(3,:),'Einit.plt')
  
  CALL COMPUTE_DT(un)
 
  tps = user_time()
  it_max = int(inputs%Tfinal/inputs%dt)
  inputs%dt = inputs%Tfinal/it_max
  tps = user_time() 
  !DO i = 1, 1
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
     !WRITE(*,*) 'mass',sum(un(1,:)*lumped),sum(un(2,:)*lumped),sum(un(3,:)*lumped)
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step and dof', tps/(it_time*mesh%np)

  !===Regression test
  CALL regression(un)

  CALL plot_1d(mesh%rr(1,:),un(1,:),'rho.plt')
  CALL plot_1d(mesh%rr(1,:),1.d0/un(1,:),'v.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:),'mt.plt')
  CALL plot_1d(mesh%rr(1,:),un(3,:),'E.plt')
  CALL plot_1d(mesh%rr(1,:),pressure(un), 'press.plt')
  CALL plot_1d(mesh%rr(1,:),un(2,:)/un(1,:),'u.plt')
  ui(1,:) = SQRT(gamma_vdw*(pressure(un)+un(1,:)**2)/(un(1,:)*(1-bvdw*un(1,:))) - 2*avdw*un(1,:))
  CALL plot_1d(mesh%rr(1,:),ui(1,:),'c.plt')
  
  IF (VdW_test_case==0) THEN
     CALL comput_errors(un,inputs%time)
  END IF
CONTAINS

  SUBROUTINE comput_errors(un,time)
    USE mesh_handling
    USE fem_tn
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
    REAL(KIND=8) :: x0, long

    long = MAXVAL(mesh%rr(1,:))-MINVAL(mesh%rr(1,:))
    x0   = 0.d0 !MINVAL(mesh%rr(1,:))+ long/2
    ALLOCATE(r_e(mesh%np),p_e(mesh%np))
    CALL rho_v_p_vdw(x0,mesh%rr(1,:)/inputs%time,r_e,m_e,p_e)
    m_e = r_e*m_e
    E_e = (p_e+avdw*r_e**2)*(1.d0-bvdw*r_e)/(gamma_vdw-1)-avdw*r_e**2+0.5d0*m_e**2/r_e
    write(*,*) mesh%np
    s_e = SQRT(gamma_vdw*(pressure(un)+un(1,:)**2)/(un(1,:)*(1-bvdw*un(1,:))) - 2*avdw*un(1,:))
    CALL plot_1d(mesh%rr(1,:),s_e,'c.plt')

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

END PROGRAM compressible_euler
