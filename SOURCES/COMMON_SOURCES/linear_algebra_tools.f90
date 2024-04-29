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

  SUBROUTINE st_csr_vect(uu_jj, n_b, ia, ja)
    USE sorting
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(IN)  :: uu_jj
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja
    INTEGER,                 INTENT(IN)  :: n_b 
    INTEGER :: nparm=180
    INTEGER :: me, nw_u, nmax, np_u
    INTEGER :: m, ni, nj, i, j, n_a_d, ib, jb, ki, kj, np_tot
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d
    nw_u = SIZE(uu_jj, 1)
    me  = SIZE(uu_jj, 2)
    np_u = MAXVAL(uu_jj)
    np_tot = n_b*np_u
    ALLOCATE (ja_work(np_tot,nparm), ia(np_tot+1), a_d(nparm), nja(np_tot))
    ja_work = 0
    nja = 0

    !===Bloc uu x uu
    DO m = 1, me
       DO ni = 1, nw_u
          i = uu_jj(ni,m)
          DO nj = 1, nw_u
             j = uu_jj(nj,m)
             DO ki = 1, n_b
                ib = (i-1)*n_b + ki
                DO kj = 1, n_b
                   jb = (j-1)*n_b + kj
                   IF (nja(ib)==0) THEN
                      nja(ib) = nja(ib) + 1
                      ja_work(ib,nja(ib)) = jb
                   ELSE IF (MINVAL(ABS(ja_work(ib,1:nja(ib))-jb)) /= 0) THEN
                      nja(ib) = nja(ib) + 1
                      ja_work(ib,nja(ib)) = jb
                   END IF
                END DO
             END DO
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
    ALLOCATE(ja(nmax))
    ia(1) = 1
    DO i = 1, np_tot
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_csr_vect'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i)
          STOP
       END IF
       ia(i+1) = ia(i) + nja(i)
       ja(ia(i):ia(i+1)-1) = a_d(1:nja(i))
    END DO
    DEALLOCATE ( ja_work, nja, a_d )
  END SUBROUTINE st_csr_vect

  SUBROUTINE elasticity_M (mesh, alpha, beta, stiff)
    !=== 2*alpha(<< (Dw), (D_) >> + << (Dw), (D_)^t >>)/2 + beta << (D.w), (D._) >>
    USE def_type_mesh
    USE space_dim
    USE matrix_type
    IMPLICIT NONE
    TYPE(mesh_type)                             :: mesh
    REAL(KIND=8),                 INTENT(IN)    :: alpha, beta
    TYPE(matrice_bloc),           INTENT(OUT)   :: stiff
    INTEGER, DIMENSION(mesh%gauss%n_w)          :: jj_loc

    INTEGER ::  m, l, p, ni, nj, ki, k1, kj, i, j, iloc, jloc, n_b
    REAL(KIND=8) :: x, y
    n_b = k_dim
    CALL st_csr_vect(mesh%jj, n_b, stiff%ia, stiff%ja)
    ALLOCATE(stiff%aa(SIZE(stiff%ja)))
    stiff%aa = 0.d0
    DO m = 1, mesh%me
       jj_loc=mesh%jj(:,m)
       DO ni = 1, mesh%gauss%n_w
          iloc = jj_loc(ni)
          DO ki = 1, k_dim
             i = (iloc-1)*n_b + ki
             DO nj = 1, mesh%gauss%n_w
                jloc = jj_loc(nj)
                DO kj = 1, k_dim
                   j = (jloc-1)*n_b + kj
                   x = 0.d0
                   DO l = 1, mesh%gauss%l_G
                      y =  alpha*mesh%gauss%dw(kj,ni,l,m)*mesh%gauss%dw(ki,nj,l,m) &
                           + beta*mesh%gauss%dw(ki,ni,l,m)*mesh%gauss%dw(kj,nj,l,m)
                      IF (kj.EQ.ki) THEN
                         DO k1 = 1, k_dim
                            y = y + alpha*mesh%gauss%dw(k1,ni,l,m)*mesh%gauss%dw(k1,nj,l,m)
                         END DO
                      END IF
                      x = x + y * mesh%gauss%rj(l,m)
                   END DO
                   DO p = stiff%ia(i),  stiff%ia(i+1) - 1
                      IF (stiff%ja(p) == j) THEN
                         stiff%aa(p) = stiff%aa(p) + x
                         EXIT
                      ENDIF
                   ENDDO

                ENDDO
             ENDDO
          ENDDO
       ENDDO
    END DO
  END SUBROUTINE elasticity_M

  SUBROUTINE stiff_11_M (mesh, stiff)
    USE def_type_mesh
    USE space_dim
    USE matrix_type
    IMPLICIT NONE
    TYPE(mesh_type)                   :: mesh
    TYPE(matrice_bloc)                :: stiff
    INTEGER,      DIMENSION(mesh%gauss%n_w) :: j_loc
    REAL(KIND = 8), DIMENSION(mesh%gauss%k_d, mesh%gauss%n_w, mesh%gauss%l_G) :: dw_loc
    REAL(KIND=8) :: x
    INTEGER :: m, l, ni, nj, i, j, p
    DO m = 1, mesh%me
       j_loc = mesh%jj(:,m)
       dw_loc = mesh%gauss%dw(:,:,:,m)
       DO ni = 1, mesh%gauss%n_w
          i = j_loc(ni)
          DO nj = 1, mesh%gauss%n_w
             j = j_loc(nj)
             x = 0.d0
             DO l = 1, mesh%gauss%l_G
                x = x + SUM(dw_loc(:,nj,l)*dw_loc(:,ni,l))*mesh%gauss%rj(l,m)
             END DO
             DO p = stiff%ia(i),  stiff%ia(i+1) - 1
                IF (stiff%ja(p) == j) THEN
                   stiff%aa(p) = stiff%aa(p) + x
                   EXIT
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE stiff_11_M

END MODULE lin_alg_tools
