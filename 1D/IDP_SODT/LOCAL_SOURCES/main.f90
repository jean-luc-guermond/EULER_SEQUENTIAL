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
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, ui
  INTEGER                                   :: i, j, it_max, it_time=0
  REAL(KIND=8)                              :: tps, to, norm
  CALL IDP_start_euler
  ALLOCATE(un(mesh%np,inputs%syst_size))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_1d(mesh%rr(1,:),un(:,1),'rhoinit.plt')
  CALL IDP_COMPUTE_DT(un)
  
  it_max = int(inputs%Tfinal/inputs%dt)
  inputs%dt = inputs%Tfinal/it_max
 
  tps = user_time() 
  DO WHILE(inputs%time<inputs%Tfinal)
     CALL IDP_COMPUTE_DT(un)
     IF (inputs%time + inputs%dt>inputs%Tfinal) THEN
        inputs%dt=inputs%Tfinal-inputs%time
     END IF
      !write(*,*) 'time ', time, inputs%dt
     CALL full_step_ERK(un) 
     inputs%time = inputs%time + inputs%dt
     it_time = it_time + 1
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step and dof', tps/(it_time*mesh%np)

   !===Regression test
  CALL regression(un)

  CALL plot_1d(mesh%rr(1,:),un(:,1),'rho.plt')
  CALL plot_1d(mesh%rr(1,:),1.d0/un(:,1),'v.plt')
  CALL plot_1d(mesh%rr(1,:),un(:,2),'mt.plt')
  CALL plot_1d(mesh%rr(1,:),un(:,3),'E.plt')
  CALL plot_1d(mesh%rr(1,:),pressure(un), 'press.plt')
  CALL plot_1d(mesh%rr(1,:),un(:,2)/un(:,1),'u.plt')
  
  CALL comput_errors(un,inputs%time)

CONTAINS

  SUBROUTINE comput_errors(un,time)
    USE mesh_handling
    USE fem_tn
    USE shock_tub
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2), INTENT(IN)     :: un
    REAL(KIND=8),                     INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(:), POINTER :: x_e, v_e, r_e, p_e
    REAL(KIND=8), DIMENSION(mesh%np) :: m_e, E_e, s_e, err_s
    REAL(KIND=8) :: error_rho_1, error_rho_2, error_rho_inf, error_mt_1, error_mt_2, error_mt_inf, &
         error_E_1, error_E_2, error_E_inf
    REAL(KIND=8) :: norm_rho_1, norm_rho_2, norm_rho_inf, norm_mt_1, norm_mt_2, norm_mt_inf, &
         norm_E_1, norm_E_2, norm_E_inf
    INTEGER :: itest
    
    !===Sodt
    itest = 1

    CALL sol_exa_stub(MINVAL(mesh%rr(1,:)),MAXVAL(mesh%rr(1,:)),mesh%np,itest,time,x_e,v_e,r_e,p_e)
    m_e = r_e*v_e
    E_e = p_e/(gamma-1) + m_e**2/(2*r_e)
    
    CALL ns_l1 (mesh, un(:,1) - r_e, error_rho_1)
    CALL ns_l1 (mesh,           r_e, norm_rho_1)
    CALL ns_0  (mesh, un(:,1) - r_e, error_rho_2)
    CALL ns_0  (mesh,           r_e, norm_rho_2)
    CALL ns_l1 (mesh, un(:,2) - m_e, error_mt_1)
    CALL ns_l1 (mesh,           m_e, norm_mt_1)
    CALL ns_0  (mesh, un(:,2) - m_e, error_mt_2)
    CALL ns_0  (mesh,           m_e, norm_mt_2)
    CALL ns_l1 (mesh, un(:,3) - E_e, error_E_1)
    CALL ns_l1 (mesh,           E_e, norm_E_1)
    CALL ns_0  (mesh, un(:,3) - E_e, error_E_2)
    CALL ns_0  (mesh,           E_e, norm_E_2)
    WRITE(*,*) ' Global relative error L1-norm    ', error_rho_1/norm_rho_1+error_mt_1/norm_mt_1+error_E_1/norm_E_1
    WRITE(*,*) ' Global relative error L2-norm    ', error_rho_2/norm_rho_2+error_mt_2/norm_mt_2+error_E_2/norm_E_2
    error_rho_inf = MAXVAL(ABS(un(:,1) - r_e))
    norm_rho_inf  = MAXVAL(ABS(r_e))
    error_mt_inf  = MAXVAL(ABS(un(:,2) - m_e))
    norm_mt_inf   = MAXVAL(ABS(m_e))
    error_E_inf   = MAXVAL(ABS(un(:,3) - E_e))
    norm_E_inf    = MAXVAL(ABS(E_e))
    WRITE(*,*) ' Global relative error Linfty-norm',&
         error_rho_inf/norm_rho_inf+error_mt_inf/norm_mt_inf+error_E_inf/norm_E_inf

    WRITE(*,*) ' Error density relative, L1-norm    ', error_rho_1/norm_rho_1
    WRITE(*,*) ' Error density relative, L2-norm    ', error_rho_2/norm_rho_2
    WRITE(*,*) ' Error density relative, Linfty-norm', error_rho_inf/norm_rho_inf

    CALL plot_1d(mesh%rr(1,:),r_e,'rho_exact')
    CALL plot_1d(mesh%rr(1,:),un(:,1)-r_e,'err_rho')
    CALL plot_1d(mesh%rr(1,:),m_e,'mt_exact')
    CALL plot_1d(mesh%rr(1,:),un(:,2)-m_e,'err_mt')
    CALL plot_1d(mesh%rr(1,:),E_e,'E_exact')
    CALL plot_1d(mesh%rr(1,:),un(:,3)-E_e,'err_E')
    s_e = (E_e - m_e**2/(2*r_e))/r_e**gamma
    err_s = ABS(s_e -(un(:,k_dim+2)-un(:,2)**2/(2.d0*un(:,1)))/(un(:,1)**gamma))
    WRITE(*,*) ' Error log(s), Max-norm', MAXVAL(err_s)/MAXVAL(s_e)
  END SUBROUTINE comput_errors

END PROGRAM compressible_euler
