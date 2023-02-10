PROGRAM scal_cons

  USE input_data
  USE mesh_handling
  USE update
  USE boundary_conditions
  USE sub_plot
  USE lambda_module
  USE fem_tn
  USE character_strings
  USE mesh_interpolation
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, uo, ui
  INTEGER                                   :: i, it=0, it_max, it_plot
  REAL(KIND=8)                              :: tps, to, norm, t_plot, dt_plot
  CHARACTER(len=1)                          :: tit

  CALL read_my_data('data')
  CALL construct_mesh
  CALL construct_matrices
  ALLOCATE(un(k_dim+2,mesh%np),uo(k_dim+2,mesh%np),ui(k_dim+2,mesh%np))
  inputs%time =0.d0
  CALL init(un,mesh%rr)
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rhoinit.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:)/un(1,:), 'uxinit.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:)/un(1,:), 'uyinit.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(4,:), 'Einit.plt')
  CALL COMPUTE_DT(un)

  !===Plotting
  !OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
  !CALL read_until(21, '===dt_plot')
  !READ (21,*) dt_plot
  !it_plot = 0
  !t_plot = dt_plot
  !CLOSE(21)
  t_plot=-1
  !===End plotting
  
  tps = user_time()
  it_max = int(inputs%Tfinal/inputs%dt)
  inputs%dt = inputs%Tfinal/it_max
  tps = user_time() 
  DO i = 1, it_max
     !IF (i==it_max) inputs%time =inputs%Tfinal-inputs%dt
     !DO WHILE(inputs%time<inputs%Tfinal)
    
     !CALL COMPUTE_DT(un)
     !IF (inputs%time + inputs%dt>inputs%Tfinal) THEN
     !   inputs%dt=inputs%Tfinal-inputs%time
     !END IF
     it = it + 1
     
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
  WRITE(*,*) 'total time', tps, 'Time/dt/dof', tps/it/(mesh%np), 'itmax', it
  CALL comput_errors(un,inputs%time)
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'rho.plt')
CONTAINS

  SUBROUTINE BC(unext,time)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: unext
    REAL(KIND=8)                             :: time
    unext(1,rho_js_D) = sol_anal(1,mesh%rr(:,rho_js_D))
    unext(2,ux_js_D)  = sol_anal(2,mesh%rr(:,ux_js_D))
    unext(3,uy_js_D)  = sol_anal(3,mesh%rr(:,uy_js_D))
    unext(4,rho_js_D) = sol_anal(4,mesh%rr(:,rho_js_D))
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
    USE shock_tub
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     :: un
    REAL(KIND=8),                        INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(:), POINTER :: x_e, v_e, r_e, p_e
    REAL(KIND=8), DIMENSION(2)          :: r0_prof, r1_prof
    REAL(KIND=8), DIMENSION(:), POINTER :: rho1d, r1d
    INTEGER :: n1d, itest
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
    
    r0_prof(1)=MINVAL(mesh%rr(1,:))
    r0_prof(2)=0.5d0*(MINVAL(mesh%rr(1,:))+MAXVAL(mesh%rr(2,:)))
    r1_prof(1)=MAXVAL(mesh%rr(1,:))
    r1_prof(2)=0.5d0*(MINVAL(mesh%rr(1,:))+MAXVAL(mesh%rr(2,:)))

    n1d = mesh%np/3
    ALLOCATE(rho1d(n1d),r1d(n1d))
    CALL coupe(mesh,r0_prof,r1_prof,n1d,mesh%rr(1,:),'r_1d.plt',r1d(:))
    CALL coupe(mesh,r0_prof,r1_prof,n1d,un(1,:),'rho_1d.plt',rho1d(:))
    !CALL sol_exa_stub(0.d0,1.d0,101,1,.0225d0,x_e,v_e,r_e,p_e)
    CALL sol_exa_stub(r0_prof(1),r1_prof(1),n1d,itest,time,x_e,v_e,r_e,p_e)

    CALL plot_1d(r1d,r_e,'rho_exact')
    WRITE(*,*) ' Error density relative, L1-norm    ', SUM(ABS(r_e-rho1d))/SUM(ABS(r_e))

  END SUBROUTINE comput_errors


  
END PROGRAM scal_cons
