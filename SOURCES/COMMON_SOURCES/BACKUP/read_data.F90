MODULE input_data
  IMPLICIT NONE
  PUBLIC :: read_my_data
  TYPE my_data
     CHARACTER(len=200)             :: directory
     CHARACTER(len=200)             :: file_name
     INTEGER                        :: nb_dom
     INTEGER, DIMENSION(:), POINTER :: list_dom
     INTEGER                        :: type_fe
     REAL(KIND=8)                   :: Tfinal
     REAL(KIND=8)                   :: CFL
     REAL(KIND=8)                   :: ce
     LOGICAL                        :: if_entropy_viscosity
     LOGICAL                        :: if_fct_limit
     CHARACTER(LEN=10)              :: viscosity_type
     CHARACTER(LEN=4)               :: type_test
     REAL(KIND=8)                   :: dt, time
     INTEGER                        :: rho_nb_Dir_bdy, ux_nb_Dir_bdy, uy_nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: rho_Dir_list, ux_Dir_list, uy_Dir_list
     INTEGER                        :: syst_size
  END type my_data
  TYPE(my_data), PUBLIC  :: inputs
  PRIVATE 
CONTAINS
  SUBROUTINE read_my_data(data_fichier)
    USE character_strings
    USE space_dim
    IMPLICIT NONE
    INTEGER, PARAMETER           :: in_unit=21
    CHARACTER(len=*), INTENT(IN) :: data_fichier
    OPEN(UNIT = in_unit, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(in_unit, "===Name of directory for mesh file===") 
    READ (in_unit,*) inputs%directory
    CALL read_until(in_unit, "===Name of mesh file===")
    READ (in_unit,*) inputs%file_name
    CALL read_until(in_unit, '===Number of subdomains in the mesh===')
    READ(21,*) inputs%nb_dom
    ALLOCATE(inputs%list_dom(inputs%nb_dom))
    CALL read_until(21, '===List of subdomain in the mesh===')
    READ(21,*) inputs%list_dom
    CALL read_until(21, '===Type of finite element===')
    READ(21,*) inputs%type_fe
    CALL read_until(in_unit, "===Final time===") 
    READ (in_unit,*) inputs%Tfinal
    CALL read_until(in_unit, "===CFL number===") 
    READ (in_unit,*) inputs%CFL
    CALL read_until(in_unit, "===Viscosity type (galerkin, Roe, viscous) ===")
    READ (in_unit,*) inputs%viscosity_type
    CALL read_until(in_unit, "===Entropy viscosity or not===")
    READ (in_unit,*) inputs%if_entropy_viscosity
    IF (inputs%if_entropy_viscosity) THEN
       CALL read_until(in_unit, "===ce coefficient===") 
       READ (in_unit,*) inputs%ce
    ELSE
       inputs%ce=0.d0
    END IF
    CALL read_until(in_unit, "======FCT limiter (true/false)===")
    READ (in_unit,*) inputs%if_fct_limit
    CALL read_until(in_unit, "===Test case name===")
    READ (in_unit,*) inputs%type_test
    CALL read_until(in_unit, "===How many Dirichlet boundaries for rho?===")
    READ (in_unit,*)  inputs%rho_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for rho?===")
    ALLOCATE(inputs%rho_Dir_list(inputs%rho_nb_Dir_bdy))
    READ (in_unit,*) inputs%rho_Dir_list
    CALL read_until(in_unit, "===How many Dirichlet boundaries for ux?===")
    READ (in_unit,*)  inputs%ux_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for ux?===")
    ALLOCATE(inputs%ux_Dir_list(inputs%ux_nb_Dir_bdy))
    READ (in_unit,*) inputs%ux_Dir_list
    IF (k_dim==2) THEN
       CALL read_until(in_unit, "===How many Dirichlet boundaries for uy?===")
       READ (in_unit,*)  inputs%uy_nb_Dir_bdy
       CALL read_until(in_unit, "===List of Dirichlet boundaries for uy?===")
       ALLOCATE(inputs%uy_Dir_list(inputs%uy_nb_Dir_bdy))
       READ (in_unit,*) inputs%uy_Dir_list
    ELSE
       inputs%uy_nb_Dir_bdy=0
       ALLOCATE(inputs%uy_Dir_list(inputs%uy_nb_Dir_bdy))
    END IF
    CLOSE(in_unit)

    inputs%syst_size=k_dim+2 !Euler
  END SUBROUTINE read_my_data
END MODULE input_data


