PROGRAM scal_cons

  USE input_data
  USE mesh_handling
  USE update
  USE boundary_conditions
  USE sub_plot
  USE lambda_module
  USE fem_tn
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, uo, ui
  INTEGER                                   :: it, it_max, i, j, k, l, p, count = 0
  REAL(KIND=8)                              :: tps, lambda_max, pstar, norm

  CALL read_my_data('data')
  CALL construct_mesh

  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  time =0.d0
  CALL init(un,mesh%rr)

  SELECT CASE(k_dim)
  CASE(1)
     CALL plot_1d(mesh%rr(1,:),un(1,:),'rho_init.plt')
  CASE(2)
     !CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho.plt')
     !CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux.plt')
     !CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uy.plt')
     !CALL plot_scalar_field(mesh%jj, mesh%rr, un(4,:), 'E.plt')
  END SELECT

  CALL COMPUTE_DT(un)
  
  CALL construct_matrices

  tps = user_time()
  it_max = int(inputs%Tfinal/inputs%dt)

  inputs%dt = inputs%Tfinal/it_max
  DO it = 1, it_max
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
     write(*,*) 'time, dt ', inputs%time, inputs%dt
  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step', tps/it_max, 'itmax', it_max, count
  SELECT CASE(k_dim)
  CASE(1)
     CALL plot_1d(mesh%rr(1,:),un(1,:),'rho.plt')
     CALL plot_1d(mesh%rr(1,:),un(2,:),'mt.plt')
     CALL plot_1d(mesh%rr(1,:),un(3,:),'E.plt')
     CALL plot_1d(mesh%rr(1,:),(un(3,:)-0.5d0*un(2,:)**2/un(1,:))/un(1,:), 'internal_energy.plt')
     CALL plot_1d(mesh%rr(1,:),pressure(un), 'Pressure.plt') 
     CALL comput_errors(un,time)
  CASE(2)
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'ux.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uy.plt')
     CALL plot_scalar_field(mesh%jj, mesh%rr, pressure(un), 'E.plt')
     CALL ns_0 (mesh, un(1,:),  norm)
     WRITE(*,*) ' L2 norm of density ', norm
  END SELECT
CONTAINS

  SUBROUTINE COMPUTE_DT(un)
    USE input_data
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN) :: un
    REAL(KIND=8)                                   :: v2_over_h2, speed 
    INTEGER                                        :: m
    REAL(KIND=8), DIMENSION(mesh%np)               :: cn
    cn = sound_speed(un) + ABS(un(2,:))/un(1,:) 
    v2_over_h2 = -1.d0
    DO m = 1, mesh%me
       speed=MAXVAL(cn(mesh%jj(:,m)))
       v2_over_h2 = MAX(v2_over_h2,speed/SUM(mesh%gauss%rj(:,m)))
    END DO
    inputs%dt = inputs%CFL/v2_over_h2
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
    USE shock_tub
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)     :: un
    REAL(KIND=8),                     INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(:), POINTER :: x_e, v_e, r_e, p_e
    REAL(KIND=8) :: error1, error2
    INTEGER :: itest, i

    SELECT CASE(inputs%type_test)
    CASE ('sodt') ! Sod
       itest = 1
    CASE ('lebl') ! Leblanc
       itest = 4
    CASE ('soni') ! Sonic
       itest = 5
    CASE DEFAULT
       WRITE(*,*) ' BUG '
       STOP
    END SELECT
    CALL sol_exa_stub(MINVAL(mesh%rr(1,:)),MAXVAL(mesh%rr(1,:)),mesh%np,itest,time,x_e,v_e,r_e,p_e)
    CALL ns_l1 (mesh, un(1,:) - r_e, error1)
    CALL ns_0  (mesh, un(1,:) - r_e, error2)
    WRITE(*,*) ' Error density, L1-norm', error1 
    WRITE(*,*) ' Error density, L2-norm', error2 
    CALL plot_1d(mesh%rr(1,:),r_e,'rho_exact')
    CALL plot_1d(mesh%rr(1,:),p_e/((gamma-1)*r_e), 'internal_energy_exact')
  END SUBROUTINE comput_errors

END PROGRAM scal_cons
