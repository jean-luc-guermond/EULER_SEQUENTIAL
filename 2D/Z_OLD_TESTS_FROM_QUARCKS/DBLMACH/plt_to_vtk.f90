program plt_to_vtk
  USE input_data
  USE mesh_handling
  USE sub_plot
  IMPLICIT NONE
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: un, grad
  
  CALL read_my_data('data')
  CALL construct_mesh
  ALLOCATE(un(mesh%np),grad(mesh%np))
  CALL read_scalar_field(mesh%jj, un, 'rho_in.plt')
  CALL schlieren(un,grad)
  CALL vtk_2d(mesh, un, 10, 'rho_in.vtk','density')
  CALL plot_scalar_field(mesh%jj, mesh%rr, grad, 'grad.plt')
   CALL vtk_2d(mesh, grad, 10, 'grad.vtk','scalars')
CONTAINS
  SUBROUTINE read_scalar_field(jj, uu, file_name)
    IMPLICIT NONE
    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: uu 
    CHARACTER(*) :: file_name
    CHARACTER(1) :: stuff
    REAL(KIND=8) :: xx
    INTEGER :: m, n1, n2, n3, unit_w=47

100 FORMAT(3(f15.8,3x))
    OPEN (UNIT=unit_w, FILE=file_name, FORM='formatted', STATUS='OLD')

    READ (unit_w, '(A)') stuff
    READ (unit_w, '(A)') stuff
    READ (unit_w, '(A)') stuff
    READ (unit_w, '(A)') stuff
    READ (unit_w, '(A)')

    IF (SIZE(jj,1)==3) THEN
       DO m = 1, SIZE(jj, 2)
          n1 = jj(1, m)
          n2 = jj(2, m)
          n3 = jj(3, m)
          READ (unit_w, *) xx, xx, uu(n1)
          READ (unit_w, *) xx, xx, uu(n2)
          READ (unit_w, *) xx, xx, uu(n3)
          READ (unit_w, 100)
       ENDDO
    ELSE IF (SIZE(jj,1)==6) THEN
       DO m = 1, SIZE(jj, 2)        
          n1 = jj(1, m)
          n2 = jj(6, m)
          n3 = jj(5, m)
          READ (unit_w, *) xx, xx, uu(n1)
          READ (unit_w, *) xx, xx, uu(n2)
          READ (unit_w, *) xx, xx, uu(n3)
          READ (unit_w, 100)

          n1 = jj(2, m)
          n2 = jj(4, m)
          n3 = jj(6, m)
          READ (unit_w, *) xx, xx, uu(n1)
          READ (unit_w, *) xx, xx, uu(n2)
          READ (unit_w, *) xx, xx, uu(n3)
          READ (unit_w, 100)

          n1 = jj(4, m)
          n2 = jj(3, m)
          n3 = jj(5, m)
          READ (unit_w, *) xx, xx, uu(n1)
          READ (unit_w, *) xx, xx, uu(n2)
          READ (unit_w, *) xx, xx, uu(n3)
          READ (unit_w, 100)

          n1 = jj(5, m)
          n2 = jj(6, m)
          n3 = jj(4, m)
          READ (unit_w, *) xx, xx, uu(n1)
          READ (unit_w, *) xx, xx, uu(n2)
          READ (unit_w, *) xx, xx, uu(n3)
          READ (unit_w, 100)
       ENDDO

    ELSE
       WRITE(*,*) ' Problem in plot_scalar_field '
       STOP
    END IF

    CLOSE(UNIT=unit_w)

  END SUBROUTINE read_scalar_field

  SUBROUTINE schlieren(un,grad)
    USE space_dim
    USE matrix_type
    USE mesh_handling
    USE st_matrix
    USE fem_s_M
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np)     :: un
    REAL(KIND=8), DIMENSION(mesh%np)     :: grad
    TYPE(matrice_bloc), DIMENSION(k_dim) :: cij
    REAL(KIND=8), DIMENSION(:), POINTER  :: lumped
    TYPE(matrice_bloc)                   :: mass
    REAL(KIND=8), DIMENSION(k_dim) :: xx
    INTEGER :: k, m, ni, i, nj, j, p, mn, mx

    !===cij(1), cij(2)
    DO k = 1, k_dim
       CALL st_csr(mesh%jj, cij(k)%ia, cij(k)%ja)
       ALLOCATE(cij(k)%aa(SIZE(cij(k)%ja)))
       cij(k)%aa = 0.d0
    END DO

    DO m = 1, mesh%me
       DO ni = 1, mesh%gauss%n_w  
          i = mesh%jj(ni, m)
          DO nj = 1, mesh%gauss%n_w  
             j = mesh%jj(nj, m)
             DO k = 1, k_dim
                xx(k) = SUM(mesh%gauss%dw(k,nj,:,m) * mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
             END DO
             DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
                IF (cij(1)%ja(p) == j) THEN 
                   DO k = 1, k_dim
                      cij(k)%aa(p) = cij(k)%aa(p) + xx(k)
                   END DO
                   EXIT
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !===mass
    CALL st_csr(mesh%jj, mass%ia, mass%ja)
    ALLOCATE(mass%aa(SIZE(mass%ja)))
    mass%aa = 0.d0
    CALL qs_00_M (mesh, 1.d0, mass%ia, mass%ja, mass%aa)

    !===lumped
    ALLOCATE(lumped(mesh%np))
    DO i = 1, mesh%np
       lumped(i) = SUM(mass%aa(mass%ia(i):mass%ia(i+1)-1))
    END DO

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

end program plt_to_vtk
