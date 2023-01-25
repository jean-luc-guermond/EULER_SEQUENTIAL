MODULE lin_alg_tools
CONTAINS
    SUBROUTINE compute_cij(mesh,cij)
    USE def_type_mesh
    USE matrix_type
    USE space_dim
    USE st_matrix
    USE euler_bc_arrays
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    TYPE(matrice_bloc), DIMENSION(k_dim):: cij
    REAL(KIND=8), DIMENSION(k_dim) :: xx
    INTEGER :: k, m, i, j, ni, nj, p
    INTEGER :: js, ms, nsi, nsj
 
    DO k = 1, k_dim
       CALL st_csr(mesh%jj, cij(k)%ia, cij(k)%ja)
       ALLOCATE(cij(k)%aa(SIZE(cij(k)%ja)))
       cij(k)%aa = 0.d0
    END DO

    IF (k_dim==1) THEN
       DO m = 1, mesh%me
          DO ni = 1, mesh%gauss%n_w  
             i = mesh%jj(ni, m)
             DO nj = 1, mesh%gauss%n_w  
                j = mesh%jj(nj, m)
                DO k = 1, k_dim
                   xx(k) = SUM(mesh%gauss%dw(k,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
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
    ELSE
       !===\int (\grad phi_j) phi_i dx + bdy term 
       DO m = 1, mesh%me
          DO ni = 1, mesh%gauss%n_w  
             i = mesh%jj(ni, m)
             DO nj = 1, mesh%gauss%n_w  
                j = mesh%jj(nj, m)
                DO k = 1, k_dim
                   xx(k) = SUM(mesh%gauss%dw(k,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
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
       !===Boundary term
       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(mesh%sides(ms)-inputs%udotn_zero_list)) == 0) THEN
             DO nsi = 1, mesh%gauss%n_ws
                i = mesh%jjs(nsi, ms)
                DO nsj = 1, mesh%gauss%n_ws
                   j = mesh%jjs(nsj, ms)
                   js = mesh%iis(nsj, ms)
                   DO k = 1, k_dim
                      xx(k) = SUM(mesh%gauss%wws(nsi,:)*mesh%gauss%wws(nsj,:)*mesh%gauss%rjs(:,ms) &
                           * (surf_normal_vtx(k,js)-mesh%gauss%rnorms(k,:,ms)))
                   END DO
                   DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
                      IF (cij(1)%ja(p) == j) THEN 
                         DO k = 1, k_dim
                            cij(k)%aa(p) = cij(k)%aa(p) + xx(k)
                         END DO
                         EXIT
                      ENDIF
                   ENDDO
                END DO
             END DO
          END IF
       END DO
    END IF
  END SUBROUTINE compute_cij

  SUBROUTINE Axb(a,b,x)
    USE matrix_type
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(IN) :: a
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: b
    REAL(KIND=8), DIMENSION(:), INTENT(OUT):: x
    INTEGER :: i, ps, pe
    DO i = 1, SIZE(a%ia)-1
       ps = a%ia(i)
       pe = a%ia(i+1)-1
       x(i) = SUM(a%aa(ps:pe)*b(a%ja(ps:pe)))
    END DO
  END SUBROUTINE Axb

  SUBROUTINE Axb_ij(a,b,x)
    USE matrix_type
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(IN) :: a
    TYPE(matrice_bloc), INTENT(OUT) :: x
    REAL(KIND=8), DIMENSION(:), INTENT(IN):: b
    INTEGER :: i, ps, pe
    DO i = 1, SIZE(a%ia)-1
       ps = a%ia(i)
       pe = a%ia(i+1)-1
       x%aa(ps:pe) = a%aa(ps:pe)*(b(a%ja(ps:pe))-b(i))
    END DO
  END SUBROUTINE Axb_ij


  SUBROUTINE compute_mass(mesh,mass)
    USE matrix_type
    USE def_type_mesh
    USE st_matrix
    USE fem_s_M
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    TYPE(matrice_bloc)                  :: mass
    CALL st_csr(mesh%jj, mass%ia, mass%ja)
    ALLOCATE(mass%aa(SIZE(mass%ja)))
    mass%aa = 0.d0
    CALL qs_00_M (mesh, 1.d0, mass%ia, mass%ja, mass%aa)
  END SUBROUTINE compute_mass

  SUBROUTINE compute_stiffness(mesh,stiff)
    USE matrix_type
    USE def_type_mesh
    USE st_matrix
    USE fem_s_M
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    TYPE(matrice_bloc)                  :: stiff
    CALL st_csr(mesh%jj, stiff%ia, stiff%ja)
    ALLOCATE(stiff%aa(SIZE(stiff%ja)))
    stiff%aa = 0.d0
    CALL qs_11_M (mesh, 1.d0, stiff%ia, stiff%ja, stiff%aa)
  END SUBROUTINE compute_stiffness

  SUBROUTINE lumped_mass(mesh,mass,lumped)
    USE matrix_type
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    TYPE(matrice_bloc)                  :: mass
    REAL(KIND=8), DIMENSION(:), POINTER :: lumped
    INTEGER :: i, np
    ALLOCATE(lumped(mesh%np))
    np = SIZE(mass%ia)-1
    ALLOCATE(lumped(np))
    DO i = 1, np
       lumped(i) = SUM(mass%aa(mass%ia(i):mass%ia(i+1)-1))
    END DO
  END SUBROUTINE lumped_mass

  SUBROUTINE diag_mat(ia,ja,diag)
    IMPLICIT NONE
    INTEGER, DIMENSION(:)          :: ia,ja
    INTEGER, DIMENSION(:), POINTER :: diag
    INTEGER :: np, i, p
    np = SIZE(ia)-1
    ALLOCATE(diag(np))
    DO i = 1, np
       DO p = ia(i), ia(i+1) - 1
          IF (i==ja(p)) THEN
             diag(i) = p
             EXIT
          END IF
       END DO
    END DO
  END SUBROUTINE diag_mat


END MODULE lin_alg_tools
