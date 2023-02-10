MODULE  update_euler_poisson
  USE space_dim
  USE mesh_handling
  USE matrix_type
  USE input_data
  IMPLICIT NONE
  TYPE(matrice_bloc)                   :: stiff, mass, uu_phi
  TYPE(matrice_bloc), DIMENSION(k_dim) :: cij
  REAL(KIND=8), DIMENSION(:), POINTER  :: lumped
  INTEGER,      DIMENSION(:), POINTER  :: diag
CONTAINS

  SUBROUTINE construct_euler_poisson_matrices
    USE st_matrix
    USE fem_s_M
    IMPLICIT NONE
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

    !===diag
    ALLOCATE(diag(mesh%np))
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          IF (i==mass%ja(p)) THEN
             diag(i) = p
             EXIT
          END IF
       END DO
    END DO
    
    !===stiff
    CALL st_csr(mesh%jj, stiff%ia, stiff%ja)
    ALLOCATE(stiff%aa(SIZE(stiff%ja)))
    stiff%aa = 0.d0
    CALL qs_11_M (mesh, 1.d0, stiff%ia, stiff%ja, stiff%aa)
    
    !===csr for pp x phi system
    CALL st_csr_uux2_ppx1(mesh, mesh, uu_phi)
    CALL mat_uu_pp(mesh,uu_phi)
    
  END SUBROUTINE construct_euler_poisson_matrices

  SUBROUTINE mat_uu_pp(mesh,uu_phi,rho,redo)
    USE def_type_mesh
    USE matrix_type
    IMPLICIT NONE
    TYPE(mesh_type),                  INTENT(IN) :: mesh
    REAL(KIND=8), DIMENSION(mesh%np), INTENT(IN) :: rho
    LOGICAL,                          INTENT(IN) :: redo
    TYPE(matrice_bloc),               INTENT(OUT):: uu_phi
    INTEGER :: i, i_b, k, p, ps, pe
    
    lt_over_two = inputs%lambda*inputs%dt/2
    np = mesh%np
    
    !===uu x phi block
    DO i = 1, mesh%np 
       DO k = 1, k_dim
          i_b = i + (k-1)*np
          ps = uu_phi%ia(i_b) + k_dim*(mass%ia(i+1)-mass%ia(i))
          !test
          IF (uu_phi%ia(i_b+1) -ps .NE. mass%ia(i+1)-mass%ia(i)) THEN
             WRITE(*,*) ' bug'
             STOP
          END IF
          !test
          DO p = mass%ia(i), mass%ia(i+1) - 1
             !test
             IF (p+ps>uu_phi%ia(i_b+1)-1)  THEN
                WRITE(*,*) ' bug p+ps>uu_phi%ia(i_b+1)-1'
                STOP
             END IF
             !test
             uu_phi%aa(ps+p) = lt_over_two*rho(i)*cij(k)%aa(p)
          END DO
       END DO
    END DO
    
    !===Return if redo is true. Rest of the matrix is already computed
    IF (redo) RETURN 
    
    !===uu x uu block
    DO i = 1, mesh%np
       DO k = 1, k_dim
          i_b = i + (k-1)*np
          p = uu_phi%ia(i_b) + diag(i) - 1
          !test
          IF (uu_phi%ja(p) .NE. i_b) THEN
             WRITE(*,*) ' Bug uu_phi%ja(p) .ne. i_b'
             STOP
          END IF
          !test
          uu_phi%aa(p) = lumped(i)*(1.d0+lt_over_two)
       END DO
    END DO
    
    !===phi x u block
    DO i = 1, mesh%np 
       i_b = k_dim*np + 1
       DO k = 1, k_dim
          ps = uu_phi%ia(i_b) + (k-1)*(mass%ia(i+1)-mass%ia(i))
          pe = ps + (mass%ia(i+1)-mass%ia(i)) - 1
          uu_phi%aa(ps:pe) = lt_over_two*cij(k)%aa
       END DO
    END DO
    
    !===phi x u block
    DO i = 1, mesh%np 
       i_b = k_dim*np + 1
       ps = uu_phi%ia(i_b) + (k_dim-1)*(mass%ia(i+1)-mass%ia(i))
       pe = ps + (mass%ia(i+1)-mass%ia(i)) - 1
       uu_phi%aa(ps:pe) = stiff%aa
    END DO
  END SUBROUTINE mat_uu_pp
  
  SUBROUTINE st_csr_uu_pp(uu_mesh, pp_mesh, mat)
    USE def_type_mesh
    USE matrix_type
    IMPLICIT NONE
    INTEGER :: nparm=80*(k+dim+1)
    TYPE(mesh_type),         INTENT(IN)  :: uu_mesh, pp_mesh
    TYPE(matrice_bloc)                   :: mat
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d
    INTEGER :: nmax
    INTEGER :: m, ni, nj, np_u, np_p, np_tot, i, i_b, j, j_b, ki, kj, n_a_d
    IF (SIZE(uu_mesh%jj, 2) .NE. SIZE(pp_mesh%jj, 2)) THEN
       WRITE(*,*) 'BUG st_csr_uux2_ppx1: mesh sizes are incompatible'
       STOP
    END IF
    np_u = uu_mesh%np
    np_p = pp_mesh%np
    np_tot = k_dim*np_u+np_p
    ALLOCATE (ja_work(np_tot,nparm), a_d(nparm), nja(np_tot))
    ja_work = 0
    nja = 0
    DO m = 1, uu_mesh%me
       DO ni = 1, uu_mesh%gauss%n_w
          i = uu_mesh%jj(ni,m)
          DO ki = 1, k_dim
             i_b = i + (ki-1)*np_u
             !===Bloc uxu
             DO nj = 1, uu_mesh%gauss%n_w
                j = uu_mesh%jj(nj,m)
                DO kj = 1, k_dim
                   j_b = j + (kj-1)*np_u
                   IF (nja(i_b)==0) THEN
                      nja(i_b) = nja(i_b) + 1
                      ja_work(i_b,nja(i_b)) = j_b
                   ELSE IF (MINVAL(ABS(ja_work(i_b,1:nja(i_b))-j_b)) /= 0) THEN 
                      nja(i_b) = nja(i_b) + 1
                      ja_work(i_b,nja(i_b)) = j_b
                   END IF
                END DO
             END DO
             !===Bloc uxp
             DO nj = 1, pp_mesh%gauss%n_w
                j = pp_mesh%jj(nj,m)
                j_b = j + k_dim*np_u
                IF (nja(i_b)==0) THEN
                   nja(i_b) = nja(i_b) + 1
                   ja_work(i_b,nja(i_b)) = j_b
                ELSE IF (MINVAL(ABS(ja_work(i_b,1:nja(i_b))-j_b)) /= 0) THEN 
                   nja(i_b) = nja(i_b) + 1
                   ja_work(i_b,nja(i_b)) = j_b
                END IF
             END DO
          END DO
       END DO
       DO ni = 1, pp_mesh%gauss%n_w
          i = pp_mesh%jj(ni,m)
          i_b = i + k_dim*np_u
          !===Bloc pxu
          DO nj = 1, uu_mesh%gauss%n_w
             j = uu_mesh%jj(nj,m)
             DO kj = 1, k_dim
                j_b = j + (kj-1)*np_u
                IF (nja(i_b)==0) THEN
                   nja(i_b) = nja(i_b) + 1
                   ja_work(i_b,nja(i_b)) = j_b
                ELSE IF (MINVAL(ABS(ja_work(i_b,1:nja(i_b))-j_b)) /= 0) THEN 
                   nja(i_b) = nja(i_b) + 1
                   ja_work(i_b,nja(i_b)) = j_b
                END IF
             END DO
          END DO
          !===Bloc pxp
          DO nj = 1, pp_mesh%gauss%n_w
             j = pp_mesh%jj(nj,m)
             j_b = j + k_dim*np_u
             IF (nja(i_b)==0) THEN
                nja(i_b) = nja(i_b) + 1
                ja_work(i_b,nja(i_b)) = j_b
             ELSE IF (MINVAL(ABS(ja_work(i_b,1:nja(i_b))-j_b)) /= 0) THEN 
                nja(i_b) = nja(i_b) + 1
                ja_work(i_b,nja(i_b)) = j_b
             END IF
          END DO
       END DO
    END DO
    IF (MAXVAL(nja)>nparm) THEN 
       WRITE(*,*) 'ST_SPARSEKIT: dimension of ja must be >= ',nparm
       STOP
    END IF
    nmax = 0
    DO i = 1, np_tot
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(mat%ja(nmax))
    ALLOCATE(mat%ia(np_tot+1))
    mat%ia(1) = 1
    DO i = 1, np_tot
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_p1_CSR'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i) 
          STOP
       END IF
       mat%ia(i+1) = mat%ia(i) + nja(i)
       mat%ja(mat%ia(i):mat%ia(i+1)-1) = a_d(1:nja(i))
    END DO
    DEALLOCATE (ja_work, nja, a_d)
  END SUBROUTINE st_csr_uu_pp

END MODULE update_euler_poisson

