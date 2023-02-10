MODULE update
  USE matrix_type
  USE space_dim
  USE input_data
  USE mesh_handling
  PUBLIC:: euler, construct_matrices, compute_dij
  PRIVATE
  LOGICAL, PUBLIC                              :: plot
  TYPE(matrice_bloc), DIMENSION(k_dim), PUBLIC :: cij, testcij
  TYPE(matrice_bloc), PUBLIC                   :: dij
  TYPE(matrice_bloc)                           :: dijH, lij
  TYPE(matrice_bloc), PUBLIC                   :: mass, stiff
  TYPE(matrice_bloc)                           :: mc_minus_ml
  TYPE(matrice_bloc), DIMENSION(k_dim+2)       :: fctmat
  TYPE(matrice_bloc)                           :: resij
  INTEGER, DIMENSION(:), POINTER, PUBLIC       :: diag
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC  :: lumped
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC  :: urelaxi, drelaxi
  INTEGER                                      :: isolve
  INTEGER                                      :: max_nb_pt_stencil
  REAL(KIND=8), PARAMETER                      :: small=1.d-8
  REAL(KIND=8), PARAMETER                      :: urelax=1.001, drelax=.999d0   
CONTAINS 

  SUBROUTINE construct_matrices
    USE st_matrix
    USE fem_s_M
    IMPLICIT NONE
    INTEGER :: m, p, ni, nj, i, j, k, nsi, nsj, ms
    REAL(KIND=8), DIMENSION(k_dim) :: xx
    REAL(KIND=8) :: dd
    !===mass
    CALL st_csr(mesh%jj, mass%ia, mass%ja)
    ALLOCATE(mass%aa(SIZE(mass%ja)))
    mass%aa = 0.d0
    CALL qs_00_M (mesh, 1.d0, mass%ia, mass%ja, mass%aa)

    !===fctmat
    DO k = 1, k_dim+2
       CALL st_csr(mesh%jj, fctmat(k)%ia, fctmat(k)%ja)
       ALLOCATE(fctmat(k)%aa(SIZE(fctmat(k)%ja)))
       fctmat(k)%aa = 0.d0
    END DO

    !===lumped
    ALLOCATE(lumped(mesh%np))
    DO i = 1, mesh%np
       lumped(i) = SUM(mass%aa(mass%ia(i):mass%ia(i+1)-1))
    END DO

    !===stiff
    CALL st_csr(mesh%jj, stiff%ia, stiff%ja)
    ALLOCATE(stiff%aa(SIZE(stiff%ja)))
    stiff%aa = 0.d0
    CALL qs_11_M (mesh, 1.d0, stiff%ia, stiff%ja, stiff%aa)
    !DO i = 1, mesh%np
    !   stiff%aa(mass%ia(i):mass%ia(i+1)-1) = stiff%aa(mass%ia(i):mass%ia(i+1)-1)*lumped(i)
    !END DO

    !===urelaxi, drelaxi
    ALLOCATE(urelaxi(mesh%np), drelaxi(mesh%np))
    DD = SUM(lumped)
    IF (k_dim==1) THEN
       !urelaxi = urelax !LEBLANC
       urelaxi = MIN(1.d0 + 2*(lumped/DD)*SQRT(lumped/DD),urelax)
       !drelaxi = drelax !LEBLANC
       drelaxi = MAX(1.d0 - 2*(lumped/DD)*SQRT(lumped/DD),drelax)
    ELSE IF (k_dim==2) THEN
       urelaxi = MIN(1.d0 + 2*SQRT(SQRT(lumped/DD))**3,urelax)
       drelaxi = MAX(1.d0 - 2*SQRT(SQRT(lumped/DD))**3,drelax)
    END IF

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

    !===Mass - lumped
    CALL st_csr(mesh%jj, mc_minus_ml%ia, mc_minus_ml%ja)
    ALLOCATE(mc_minus_ml%aa(SIZE(mc_minus_ml%ja)))
    mc_minus_ml%aa = mass%aa
    mc_minus_ml%aa(diag) = mc_minus_ml%aa(diag) - lumped

    !===dij
    CALL st_csr(mesh%jj, dij%ia, dij%ja)
    ALLOCATE(dij%aa(SIZE(dij%ja)))
    dij%aa = 0.d0

    !===dijH
    CALL st_csr(mesh%jj, dijH%ia, dijH%ja)
    ALLOCATE(dijH%aa(SIZE(dijH%ja)))
    dijH%aa = 0.d0

    !===lij
    CALL st_csr(mesh%jj, lij%ia, lij%ja)
    ALLOCATE(lij%aa(SIZE(lij%ja)))
    lij%aa = 0.d0

    !===cij(1), cij(2)
    DO k = 1, k_dim
       CALL st_csr(mesh%jj, cij(k)%ia, cij(k)%ja)
       ALLOCATE(cij(k)%aa(SIZE(cij(k)%ja)))
       cij(k)%aa = 0.d0
    END DO

    !IF (k_dim==1 .or. k_dim==2) THEN
    IF (k_dim==1) THEN
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
    ELSE
       !===Change \int (\grad phi_i) phi_j dx 
       DO m = 1, mesh%me
          DO ni = 1, mesh%gauss%n_w  
             i = mesh%jj(ni, m)
             DO nj = 1, mesh%gauss%n_w  
                j = mesh%jj(nj, m)
                DO k = 1, k_dim
                   xx(k) = - SUM(mesh%gauss%dw(k,ni,:,m) * mesh%gauss%ww(nj,:)*mesh%gauss%rj(:,m))
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

       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(mesh%sides(ms)-inputs%udotn_zero_list)) == 0) THEN
             DO nsi = 1, mesh%gauss%n_ws
                i = mesh%jjs(nsi, ms)
                DO nsj = 1, mesh%gauss%n_ws
                   j = mesh%jjs(nsj, ms)
                   DO k = 1, k_dim
                      xx(k) = SUM(mesh%gauss%wws(nsi,:)*mesh%gauss%wws(nsj,:)*mesh%gauss%rjs(:,ms))*(normal_v(k,j))
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
          ELSE
             DO nsi = 1, mesh%gauss%n_ws
                i = mesh%jjs(nsi, ms)
                DO nsj = 1, mesh%gauss%n_ws
                   j = mesh%jjs(nsj, ms)
                   DO k = 1, k_dim
                      xx(k) = SUM(mesh%gauss%wws(nsi,:)*mesh%gauss%wws(nsj,:)*mesh%gauss%rjs(:,ms)*(mesh%gauss%rnorms(k,:,ms)))
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

    !TEST=================1D
    !cij(2)%aa = 1.d-30
    !TEST=================1D

    !===entropy viscosity matrix
    CALL st_csr(mesh%jj, resij%ia, resij%ja)
    ALLOCATE(resij%aa(SIZE(resij%ja)))
    resij%aa = 0.d0

    !===Maximum number of points in stencil
    max_nb_pt_stencil = 0
    DO i = 1, mesh%np
       max_nb_pt_stencil = MAX(max_nb_pt_stencil, mass%ia(i+1)-mass%ia(i))
    END DO

  END SUBROUTINE construct_matrices

  SUBROUTINE euler(un,unext)
    USE mesh_handling
    USE boundary_conditions
    USE input_data
    USE lambda_module
    USE lambda_module_full
    USE pardiso_solve
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np)              :: ulow, du, rkgal
    REAL(KIND=8), DIMENSION(mesh%np)                      :: rhomax, rhomin, Emax, Emin, ff
    REAL(KIND=8), DIMENSION(mesh%np)                      :: ccmin, ccmax, cc, smin, Esmall, rho_s_small
    REAL(KIND=8), DIMENSION(mesh%np)                      :: rho_e, rho_e_min, rho_e_max, rho_e_small
    REAL(KIND=8), DIMENSION(mesh%np)                      :: alpha_rho, alpha_rho_e, denom_rho, denom_rho_e
    REAL(KIND=8), DIMENSION(k_dim+2,max_nb_pt_stencil)    :: ubar
    REAL(KIND=8), DIMENSION(max_nb_pt_stencil)            :: rho_e_bar, rhou2bar
    REAL(KIND=8), DIMENSION(k_dim+2,k_dim,mesh%np)  :: vv
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np)        :: rk
    INTEGER,      DIMENSION(mesh%gauss%n_w) :: jloc
    REAL(KIND=8), DIMENSION(k_dim+2)        :: ul, ulp
    REAL(KIND=8) rhobar, Ebar, rkloc, sl, slp
    INTEGER :: comp, p, p_start, p_end, i, j, ij, k, l, m, it
    LOGICAL, SAVE :: once=.TRUE.
    IF (once) THEN
       isolve=-1
       once=.FALSE.
    END IF

    !===Viscosity
    SELECT CASE(inputs%viscosity_type)
    CASE('galerkin') 
       dij%aa = 0.d0
    CASE('Roe')
       CALL compute_Roe_dij(un)
    CASE('viscous') 
       CALL compute_dij(un)
    CASE DEFAULT
       WRITE(*,*) ' BUG in euler'
       STOP
    END SELECT

    !===First-order viscosity
    CALL maxmin(un(1,:),dij,rhomax,rhomin)
    CALL maxmin(un(k_dim+2,:),dij,Emax,Emin)
    CALL estimate_rho_e_min_max(un,rho_e,rho_e_min,rho_e_max)
    !TEST Bojan's proposal
    rhomax = un(1,:)
    rhomin = un(1,:)
    !TEST
    
    DO i = 1, mesh%np
       cc(i) = (un(k_dim+2,i)- SUM(un(2:k_dim+1,i)**2)/(2*un(1,i)))/(un(1,i)**gamma)
    END DO
    DO i = 1, mesh%np
       ccmin(i) = MINVAL(cc(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
    END DO
    DO comp = 1, k_dim+2
       vv(comp,:,:)=flux(comp,un)
    END DO
    alpha_rho=0.d0
    alpha_rho_e=0.d0
    denom_rho=0.d0
    denom_rho_e=0.d0
    rk = 0.d0
    rkgal=0.d0
    DO i = 1, mesh%np
       p_start = mass%ia(i)
       p_end   = mass%ia(i+1)-1
       ubar = 0.d0
       ij = 0
       DO p = p_start, p_end
          ij = ij + 1
          j = mass%ja(p)
          DO comp = 1, k_dim+2
             rkloc = 0.d0
             DO k = 1, k_dim
                rkloc  = rkloc - cij(k)%aa(p)*(vv(comp,k,j)-vv(comp,k,i))
             END DO
             rkgal(comp,i) = rkgal(comp,i) + rkloc
             rk(comp, i) = rk(comp, i) + rkloc + dij%aa(p)*(un(comp,j)-un(comp,i))
             IF (inputs%if_fct_limit) THEN
                ubar(comp,ij) = (rkloc/dij%aa(p) + un(comp,j)+un(comp,i))/2
             END IF
          END DO
          alpha_rho(i)   = alpha_rho(i)   + un(1,j) - un(1,i)
          alpha_rho_e(i) = alpha_rho_e(i) + rho_e(j) - rho_e(i)
          denom_rho(i) = un(1,j) + un(1,i)
          denom_rho_e(i) = denom_rho_e(i) + ABS(rho_e(j)) + ABS(rho_e(i))
       END DO
       IF (inputs%if_fct_limit) THEN
          rhomin(i) = MIN(rhomin(i),MINVAL(ubar(1,:)))
          rhomax(i) = MAX(rhomax(i),MAXVAL(ubar(1,:)))
          Emin(i)   = MIN(Emin(i),  MINVAL(ubar(k_dim+2,:)))
          Emax(i)   = MAX(Emax(i),  MAXVAL(ubar(k_dim+2,:)))
          ij = p_end-p_start+1
          DO k = 1, ij
             rho_e_bar(k) = ubar(k_dim+2,k) - SUM(ubar(2:k_dim+1,k)**2)/(2*ubar(1,k))
          END DO
          rho_e_max(i) = MAX(rho_e_max(i), MAXVAL(rho_e_bar(1:ij)))
          rho_e_min(i) = MIN(rho_e_min(i), MINVAL(rho_e_bar(1:ij)))
          !rhomin(i) = MIN(rhomin(i), MINVAL(un(1,mass%ja(p_end:p_start))))
          !rhomax(i) = MAX(rhomax(i), MAXVAL(un(1,mass%ja(p_end:p_start))))
          !TEST Bojan's proposal (negligible effects compared to the presence of the bar states)
          !rhomin(i) = MIN(rhomin(i), MINVAL(un(1,i)+un(1,mass%ja(p_end:p_start)))/2)
          !rhomax(i) = MAX(rhomax(i), MAXVAL(un(1,i)+un(1,mass%ja(p_end:p_start)))/2)
          !TEST
       END IF
    END DO
    !===Solve and update
    DO comp = 1, k_dim+2
       IF (inputs%viscosity_type=='galerkin') THEN
          CALL solve_pardiso(mass%aa,mass%ia,mass%ja,rk(comp,:),ff,isolve,1)
          isolve=ABS(isolve)
          unext(comp,:) = un(comp,:)+inputs%dt*ff
       ELSE
          unext(comp,:) = un(comp,:)+inputs%dt*rk(comp,:)/lumped
       END IF
    END DO

    IF (inputs%if_fct_limit) CALL check_positivity(un)

    IF (.NOT.inputs%if_entropy_viscosity .OR. inputs%viscosity_type=='galerkin') THEN
       RETURN
    END IF

    !===Continue with entropy viscosity
    ulow = unext
    CALL entropy_residual(un,rkgal)
    !CALL smoothness_viscosity(un)

    rk = rkgal
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          rk(:,i) = rk(:,i) + dijH%aa(p)*(un(:,j)-un(:,i))
       END DO
    END DO
    !===Solve and update (no limiting, lumped mass matrix)
    !IF (.NOT.inputs%if_fct_limit) THEN
    !   DO comp = 1, k_dim+2
    !      unext(comp,:) = un(comp,:)+inputs%dt*rk(comp,:)/lumped
    !   END DO
    !   RETURN
    !END IF
    
    !===Solve and update (limiting and consistent mass matrix)
    IF (inputs%high_solve=='high_consist') THEN 
       DO comp = 1, k_dim+2
          CALL solve_pardiso(mass%aa,mass%ia,mass%ja,rk(comp,:),ff,isolve,1)
          isolve=ABS(isolve)
          unext(comp,:) = un(comp,:)+inputs%dt*ff
       END DO
    ELSE IF (inputs%high_solve=='high_lumped') THEN
       DO comp = 1, k_dim+2
          unext(comp,:) = un(comp,:)+inputs%dt*rk(comp,:)/lumped
       END DO
    ELSE
       WRITE(*,*) 'BUG: inputs%high_solve not well deffined'
       STOP
    END IF
       
    !===FCT correction
    IF (inputs%if_fct_limit) THEN
       IF (inputs%high_solve=='high_lumped') THEN
          du = 0.d0
       ELSE
          du = unext-un
       END IF
       IF (inputs%if_relax_rho) THEN
          CALL relax(un(1,:),rhomin,rhomax)
       END IF
       CALL relax(un(k_dim+2,:),Emax,Emin)
       !CALL relax(rho_e,rho_e_min,rho_e_max)
       CALL relax_cmin(un,ccmin)
       !CALL estimate_cmin(un,ccmin)
       smin = LOG(ccmin)
       DO i = 1, mesh%np
          Esmall(i)= small*MINVAL(un(k_dim+2,mass%ja(mass%ia(i):mass%ia(i+1)-1)))
       END DO
       rho_e_small = small*ABS(rho_e_max)
       rho_s_small = small*ABS(rhomax*smin)
       !===Fct matrix, viscosity + mass matrix correction
       CALL compute_fct_matrix_full(un,du)
       DO it = 1, 2 ! 2 and more is good for large CFL
          !===Limit density
          !CALL FCT_generic(ulow,rhomax,rhomin,fctmat(1),lumped,1)
          CALL LOCAL_limit(ulow,rhomax,rhomin,fctmat(1),lumped,1) !===Works best
          !===Limit Total energy (Kills second order in max norm on 1D SMOOTH problem)
          !CALL FCT_generic(ulow,Emax,Emin,fctmat(k_dim+2),lumped,k_dim+2)
          !CALL LOCAL_limit(ulow,Emax,Emin,fctmat(k_dim+2),lumped,k_dim+2)
          !===Limit rho*e - (rho*e)_min
          !CALL convex_limiting(ulow,un,rho_e_min,rho_e_small,psi_rho_e,psi_rho_e_prime)
          !===Limit rho*e - e^(s_min) rho^gamma
          CALL limit_specific_entropy(ulow,un,ccmin) !===Works best
          !CALL convex_limiting(ulow,un,ccmin,Esmall,rho_e_minus_c_rho_gamma,rho_e_minus_c_rho_gamma_prime)
          !===Limit rho (s - s_min) 
          !CALL convex_limiting(ulow,un,smin,rho_s_small,rho_s,rho_s_prime)
          !===Tranpose lij
          CALL transpose_op(lij,'min')
          CALL update_fct_full(ulow,ulow)
          DO k = 1, k_dim+2
             fctmat(k)%aa = (1-lij%aa)*fctmat(k)%aa
          END DO
       END DO
       !WRITE(*,*) MINVAL(rhomax-ulow(1,:)), MINVAL(ulow(1,:)-rhomin) 
       unext = ulow
    END IF

    !===End of computation
    RETURN

  END SUBROUTINE euler

  SUBROUTINE LOCAL_limit_scal(unext,maxn,minn,mat,mass)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: mass, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)   :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT):: unext
    REAL(KIND=8), DIMENSION(SIZE(unext))     :: Qplus, Qminus, Pplus, Pminus, Rplus, Rminus
    INTEGER,      DIMENSION(SIZE(lij%ia))    :: iao
    LOGICAL,      DIMENSION(SIZE(unext))     :: check
    REAL(KIND=8), PARAMETER :: smallplus = 1.d-15, smallminus = -1.d-15
    REAL(KIND=8) :: x, fij, maxni, minni, ui, uij, xij, lambdai, umax, usmall
    INTEGER      :: i, j, p,  next

    !===Compute lij
    umax = MAX(MAXVAL(maxn),-MINVAL(minn))
    usmall = umax*smallplus
    lij%aa = 1.d0
    DO i = 1, SIZE(unext)
       lambdai = 1.d0/(mat%ia(i+1) - 1.d0 - mat%ia(i))
       maxni = maxn(i)
       minni = minn(i)
       ui = unext(i)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          xij = mat%aa(p)/(mass(i)*lambdai)
          uij = ui + xij
          IF (uij>maxni) THEN
             lij%aa(p) = MIN(ABS(maxni - ui)/(ABS(xij)+usmall),1.0)
          ELSE IF (uij<minni) THEN
             lij%aa(p) = MIN(ABS(minni - ui)/(ABS(xij)+usmall),1.d0)
          END IF
       END DO
    END DO
    lij%aa(diag) = 0.d0
    !===lij = min(lij,lij^T)
    CALL transpose_op(lij,'min')
    !===Update u
    DO i = 1, SIZE(unext)
       x = 0.d0
       DO p = mat%ia(i), mat%ia(i+1) - 1
          x = x + lij%aa(p)*mat%aa(p)
       END DO
       unext(i) = unext(i) + x/mass(i)
    END DO
  END SUBROUTINE LOCAL_limit_scal

  SUBROUTINE relax(un,minn,maxn)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:)              :: un
    REAL(KIND=8), DIMENSION(:)              :: minn
    REAL(KIND=8), DIMENSION(:)              :: maxn
    REAL(KIND=8), DIMENSION(SIZE(un))       :: alpha, denom
    INTEGER      :: i, j, p, ps, pe
    REAL(KIND=8) :: x, mx, mn
    
    alpha = 0.d0
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          alpha(i) = alpha(i) + un(i) - un(j)
       END DO
    END DO

    IF (inputs%limiter_type=='median') THEN !===Median
       DO i = 1, SIZE(un)
          ps = mass%ia(i)
          pe = mass%ia(i+1) - 1
          denom(i) = Median(alpha(mass%ja(ps:pe)), pe-ps+1)
       END DO
       maxn = MIN(urelaxi*maxn,maxn + ABS(denom)/2)
       minn = MAX(drelaxi*minn,minn - ABS(denom)/2) 
    ELSE IF(inputs%limiter_type=='avg') THEN !===Average
       denom = 0.d0
       DO i = 1, SIZE(un)
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = mass%ja(p)
             IF (i==j) CYCLE
             denom(i) = denom(i) + alpha(j) + alpha(i)
          END DO
       END DO
       DO i = 1, SIZE(un)
          denom(i) = denom(i)/(mass%ia(i+1)-mass%ia(i))/2
       END DO
       maxn = MIN(urelaxi*maxn,maxn + ABS(denom)/2)
       minn = MAX(drelaxi*minn,minn - ABS(denom)/2) 
    ELSE IF(inputs%limiter_type=='minmod') THEN !===Minmod
       denom = alpha    
       DO i = 1, SIZE(un)
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = mass%ja(p)
             IF (i==j) CYCLE
             IF (denom(i)*alpha(j).LE.0.d0) THEN
                denom(i) = 0.d0
             ELSE IF (ABS(denom(i)) > ABS(alpha(j))) THEN
                denom(i) = alpha(j)
             END IF
          END DO
       END DO
       maxn = MIN(urelaxi*maxn,maxn + ABS(denom)/2)
       minn = MAX(drelaxi*minn,minn - ABS(denom)/2)
    END IF
  END SUBROUTINE RELAX

  SUBROUTINE relax_cmin(un,cmin)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:)            :: un
    REAL(KIND=8), DIMENSION(:)              :: cmin
    REAL(KIND=8), DIMENSION(SIZE(cmin))     :: dc
    REAL(KIND=8), DIMENSION(k_dim+2)        :: ul
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: jloc
    INTEGER      :: i, j, m, k, p
    REAL(KIND=8) :: cl

    dc = 0.d0
!!$    DO m = 1, mesh%me
!!$       jloc = mesh%jj(:,m)
!!$       DO  k = 1, k_dim+2
!!$          ul(k) = SUM(un(k,jloc))/mesh%gauss%n_w
!!$       END DO
!!$       cl = (ul(k_dim+2)-SUM(ul(2:k_dim+1)**2)/(2.d0*ul(1)))/(ul(1)**gamma)
!!$       dc(jloc) = max(dc(jloc),ABS(cmin(jloc)-cl))
!!$    END DO
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          ul(:) = (un(:,i)+un(:,j))/2
          cl = (ul(k_dim+2)-SUM(ul(2:k_dim+1)**2)/(2.d0*ul(1)))/(ul(1)**gamma)
          dc(i) = MAX(dc(i),cl-cmin(i))
       END DO
    END DO
    cmin =  MAX(drelaxi*cmin, cmin - dc)
    
  END SUBROUTINE RELAX_cmin

  SUBROUTINE limit_specific_entropy(ulow,un,cmin)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN) :: ulow, un
    REAL(KIND=8), DIMENSION(mesh%np)                     :: cmin
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: test
    REAL(KIND=8), DIMENSION(mesh%np)  :: centrop, Esmall
    REAL(KIND=8), DIMENSION(k_dim+2)  :: ul, ur, Pij, x
    REAL(KIND=8) :: lambdai, coeff, psir, psil, ll, lr, llold, lrold, Qplus, dQplus, Budget
    INTEGER      :: i, j, p, k, Card

    !===
    DO i = 1, mesh%np
       Esmall(i)= small*MINVAL(un(k_dim+2,mass%ja(mass%ia(i):mass%ia(i+1)-1)))
       lambdai = 1.d0/(mass%ia(i+1) - 1.d0 - mass%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))

       !===Budget
       Qplus = 0.d0
       Card  = 0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          ul = ulow(:,i)
          DO k = 1 , k_dim+2
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(k,i) + lij%aa(p)*Pij(k) !===Density must be positive
          END DO
          dQplus = MIN(psi_func(ul,cmin(i),0.d0),psi_func(ur,cmin(i),0.d0))
          IF (dQplus>0.d0) THEN
             Qplus = Qplus + dQplus
          ELSE
             Card  = Card + 1
          END IF
       END DO
       IF (Card.NE.0) THEN
          Budget = -Qplus/Card
       ELSE
          Budget = -1d15*Qplus
       END IF
       Budget =0.d0
       !===End Budget

       DO p = mass%ia(i), mass%ia(i+1) - 1
          j =  mass%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p)
          DO k = 1, k_dim+2
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(k,i) + lr*Pij(k) !===Density must be positive
          END DO
          psir = psi_func(ur,cmin(i),Budget)
          IF (psir.GE.-Esmall(i)) THEN
             lij%aa(p) = lij%aa(p)
             CYCLE
          END IF
          ll = 0.d0
          ul = ulow(:,i)
          psil = psi_func(ul,cmin(i),Budget)
          DO WHILE (ABS(psil-psir) .GT. Esmall(i))
             llold = ll
             lrold = lr
             ll = ll - psil*(lr-ll)/(psir-psil)
             lr = lr - psir/psi_prime_func(Pij,ur,cmin(i))
             IF (ll.GE.lr) THEN
                ll = lr !lold
                EXIT
             END IF
             IF (ll< llold) THEN
                ll = llold
                EXIT
             END IF
             IF (lr > lrold) THEN
                lr = lrold
                EXIT
             END IF
             ul = ulow(:,i) + ll*Pij
             ur = ulow(:,i) + lr*Pij
             psil = psi_func(ul,cmin(i),Budget)
             psir = psi_func(ur,cmin(i),Budget)
          END DO
          IF (psir.GE.-Esmall(i)) THEN
             lij%aa(p) = lr
          ELSE
             lij%aa(p) = ll
          END IF
       END DO
    END DO

  CONTAINS
    FUNCTION psi_func(u,cmin,Budget) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u
      REAL(KIND=8),                     INTENT(IN) :: cmin
      REAL(KIND=8)                                 :: psi, Budget
      psi = u(k_dim+2) - SUM(u(2:k_dim+1)**2)/(2.d0*u(1)) - cmin*u(1)**gamma - Budget
    END FUNCTION psi_func
    FUNCTION psi_prime_func(Pij,u,cmin) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u, Pij
      REAL(KIND=8),                     INTENT(IN) :: cmin
      REAL(KIND=8)                                 :: psi
      psi = Pij(k_dim+2) - SUM(u(2:k_dim+1)*Pij(2:k_dim+1))/u(1) &
           + Pij(1)*SUM(u(2:k_dim+1)**2)/(2*u(1)**2) &
           - cmin*gamma*Pij(1)*u(1)**(gamma-1.d0)
    END FUNCTION psi_prime_func
  END SUBROUTINE limit_specific_entropy

  SUBROUTINE estimate_rho_e_min_max(un,rho_e,rho_e_min,rho_e_max)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(mesh%np)                     :: rho_e_min, rho_e_max
    REAL(KIND=8), DIMENSION(mesh%np)                     :: rho_e
    INTEGER :: i
    DO i = 1, mesh%np
       rho_e(i) = un(k_dim+2,i)-SUM(un(2:k_dim+1,i)**2)/(2.d0*un(1,i))
    END DO
    DO i = 1, mesh%np
       rho_e_min(i)  = MINVAL(rho_e(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
       rho_e_max(i)  = MAXVAL(rho_e(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
    END DO
  END SUBROUTINE estimate_rho_e_min_max

  SUBROUTINE estimate_cmin(un,cmin)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(mesh%np)                     :: cmin, cmax, rk, s, ds, c, dc, alpha
    REAL(KIND=8), DIMENSION(k_dim+2) :: up, ux, ul, ulp
    INTEGER, DIMENSION(mesh%gauss%n_w) :: jloc
    INTEGER :: i, j, it, n, m, j0, j1, l, k, p
    REAL(KIND=8) :: dh, xx, a, b, cl, sl, slp, sl1, sl2, num, denom
    DO i = 1, mesh%np
       c(i) = (un(k_dim+2,i)-SUM(un(2:k_dim+1,i)**2)/(2.d0*un(1,i)))/(un(1,i)**gamma)
    END DO
    DO i = 1, mesh%np
       cmin(i)  = MINVAL(c(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
       cmax(i)  = MAXVAL(c(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
    END DO

   IF (.FALSE.) THEN
       DO m = 1, mesh%me
          jloc = mesh%jj(:,m)
          DO  k = 1, k_dim+2
             ul(k) = SUM(un(k,jloc))/mesh%gauss%n_w
          END DO
          cl = (ul(k_dim+2)-SUM(ul(2:k_dim+1)**2)/(2.d0*ul(1)))/(ul(1)**gamma)
          cmin(jloc) = MIN(cmin(jloc),cmin(jloc)- ABS(cmin(jloc)-cl))
       END DO
    ELSE IF (.TRUE.) THEN
       s = LOG(c)
       dc = -1.d10
       ds = -1.d10
       DO m = 1, mesh%me
          jloc = mesh%jj(:,m)
          DO  k = 1, k_dim+2
             ul(k) = SUM(un(k,jloc))/mesh%gauss%n_w
          END DO
          cl = (ul(k_dim+2)-SUM(ul(2:k_dim+1)**2)/(2.d0*ul(1)))/(ul(1)**gamma)
          dc(jloc) = MAX(dc(jloc),MAXVAL(ABS(c(jloc)-cl)))
       END DO
       cmin = MAX(drelax*cmin, c - dc)
    ELSE
       !===Alternative definition of c_min
       DO i = 1, mesh%np
          denom =0.d0
          num  = 0.d0
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = dij%ja(p)
             num = num + c(j) - c(i)
          END DO
          IF (denom.LE.ABS(num)) THEN
             alpha(i) = 1.d0
             num = denom
          ELSE
             alpha(i) = ABS(num)/MAX(ABS(cmax(i)),ABS(cmin(i)))
          END IF
       END DO
       cmin = c*(1.d0-alpha)
       !==End alternative definition of c_min 
    END IF
  END SUBROUTINE estimate_cmin

  SUBROUTINE convex_limiting(ulow,un,cmin,psi_small,psi_func,psi_prime_func)
    USE boundary_conditions
    IMPLICIT NONE
    INTERFACE
       FUNCTION psi_func(u,cmin,Budget) RESULT(psi)
         USE boundary_conditions
         USE space_dim
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(k_dim+2) :: u
         REAL(KIND=8)                     :: cmin
         REAL(KIND=8)                     :: Budget
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_func
       FUNCTION psi_prime_func(Pij,u,cmin) RESULT(psi)
         USE boundary_conditions
         USE space_dim
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(k_dim+2) :: u, Pij
         REAL(KIND=8)                     :: cmin
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_prime_func
    END INTERFACE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN) :: ulow, un
    REAL(KIND=8), DIMENSION(mesh%np)                     :: cmin, psi_small
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: test
    REAL(KIND=8), DIMENSION(mesh%np)  :: centrop
    REAL(KIND=8), DIMENSION(k_dim+2)  :: ul, ur, Pij, x
    REAL(KIND=8) :: lambdai, coeff, psir, psil, ll, lr, llold, lrold, Qplus, dQplus, Budget
    INTEGER      :: i, j, p, k, Card, it
    DO i = 1, mesh%np
       lambdai = 1.d0/(mass%ia(i+1) - 1.d0 - mass%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))
       !===Budget
       Qplus = 0.d0
       Card  = 0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          ul = ulow(:,i)
          DO k = 1 , k_dim+2
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(k,i) + lij%aa(p)*Pij(k) !===Density must be positive
          END DO
          dQplus = MIN(psi_func(ul,cmin(i),0.d0),psi_func(ur,cmin(i),0.d0))
          IF (dQplus>0.d0) THEN
             Qplus = Qplus + dQplus
          ELSE
             Card  = Card + 1
          END IF
       END DO
       IF (Card.NE.0) THEN
          Budget = -Qplus/Card
       ELSE
          Budget = -1d15*Qplus
       END IF
!!$       Budget =0.d0
       !===End Budget
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j =  mass%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p) !===Use previous limiter
          DO k = 1, k_dim+2
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(k,i) + lr*Pij(k) 
          END DO
          psir = psi_func(ur,cmin(i),Budget)
          IF (psir.GE.-psi_small(i)) THEN
             CYCLE
          END IF
          ll = 0.d0
          ul = ulow(:,i)
          !psil = max(psi_func(ul,cmin(i),Budget),0.d0) !===To avoid roundoff negative
          psil = psi_func(ul,cmin(i),Budget)
          DO WHILE (ABS(psil-psir) .GT. psi_small(i))
             llold = ll
             lrold = lr
             ll = ll - psil*(lr-ll)/(psir-psil)
             lr = lr - psir/psi_prime_func(Pij,ur,cmin(i))
             IF (ll< llold) THEN
                !write(*,*) ' f1', ll , llold, psil, psir,psi_small(i)
                !stop
                ll = llold
                EXIT
             END IF
             IF (lr > lrold) THEN
                !write(*,*) ' f2', lr , lrold, psil, psir,psi_small(i),psi_prime_func(Pij,ur,cmin(i))
                !write(*,*), 'cmin', cmin(i)
                !stop
                lr = lrold
                EXIT
             END IF
             IF (ll.GE.lr) THEN
                ll = lr
                EXIT
             END IF
             ul = ulow(:,i) + ll*Pij
             ur = ulow(:,i) + lr*Pij
             !psil = max(psi_func(ul,cmin(i),Budget),0.d0) !===To avoid roundoff negative
             psil = psi_func(ul,cmin(i),Budget)
             psir = psi_func(ur,cmin(i),Budget)
          END DO
          IF (psir.GE.-psi_small(i)) THEN
             lij%aa(p) = lr
          ELSE
             lij%aa(p) = ll
          END IF
       END DO
    END DO
  CONTAINS
  END SUBROUTINE convex_limiting

  FUNCTION rho_e_minus_c_rho_gamma(u,cmin,Budget) RESULT(psi)
    USE boundary_conditions
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2) :: u
    REAL(KIND=8)                     :: cmin
    REAL(KIND=8)                     :: Budget
    REAL(KIND=8)                     :: psi
    psi = u(k_dim+2) - SUM(u(2:k_dim+1)**2)/(2.d0*u(1)) - cmin*u(1)**gamma - Budget
  END FUNCTION rho_e_minus_c_rho_gamma

  FUNCTION rho_e_minus_c_rho_gamma_prime(Pij,u,cmin) RESULT(psi)
    USE boundary_conditions
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2) :: u, Pij
    REAL(KIND=8)                     :: cmin
    REAL(KIND=8)                     :: psi
    psi = Pij(k_dim+2) - SUM(u(2:k_dim+1)*Pij(2:k_dim+1))/u(1) &
         + Pij(1)*SUM(u(2:k_dim+1)**2)/(2*u(1)**2) &
         - cmin*gamma*Pij(1)*u(1)**(gamma-1.d0)
  END FUNCTION rho_e_minus_c_rho_gamma_prime

  FUNCTION psi_rho_e(u,emin,Budget) RESULT(psi)
    USE boundary_conditions
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2) :: u
    REAL(KIND=8)                     :: emin
    REAL(KIND=8)                     :: Budget
    REAL(KIND=8)                     :: psi
    psi = u(k_dim+2) - SUM(u(2:k_dim+1)**2)/(2.d0*u(1)) - emin - Budget
  END FUNCTION psi_rho_e

  FUNCTION psi_rho_e_prime(Pij,u,emin) RESULT(psi)
    USE boundary_conditions
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2) :: u, Pij
    REAL(KIND=8)                     :: emin
    REAL(KIND=8)                     :: psi
    REAL(KIND=8) :: rho_e
    psi = Pij(k_dim+2) - SUM(Pij(2:k_dim+1)*u(2:k_dim+1))/u(1) &
         + Pij(1)*SUM(u(2:k_dim+1)**2)/(2.d0*u(1)**2)
  END FUNCTION psi_rho_e_prime

  FUNCTION rho_s(u,smin,Budget) RESULT(psi)
    USE boundary_conditions
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2) :: u
    REAL(KIND=8)                     :: smin
    REAL(KIND=8)                     :: Budget
    REAL(KIND=8)                     :: psi
    psi = u(1)*(LOG(u(k_dim+2) - SUM(u(2:k_dim+1)**2)/(2.d0*u(1))) - gamma*LOG(u(1)) - smin) - Budget
  END FUNCTION rho_s

  FUNCTION rho_s_prime(Pij,u,smin) RESULT(psi)
    USE boundary_conditions
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2) :: u, Pij
    REAL(KIND=8)                     :: smin
    REAL(KIND=8)                     :: psi
    REAL(KIND=8) :: rho_e
    rho_e = u(k_dim+2) - SUM(u(2:k_dim+1)**2)/(2.d0*u(1))
    psi = Pij(1)*(LOG(rho_e) -smin -gamma*(1+LOG(u(1)))) &
         + (u(1)/rho_e)*(Pij(k_dim+2) - SUM(Pij(2:k_dim+1)*u(2:k_dim+1))/u(1) &
         + Pij(1)*SUM(u(2:k_dim+1)**2)/(2.d0*u(1)**2))
  END FUNCTION rho_s_prime

  SUBROUTINE compute_fct_matrix_full(un,du)
    USE mesh_handling
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN) :: un, du
    INTEGER :: i, j, k, p
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          DO k = 1 , k_dim+2
             fctmat(k)%aa(p) = inputs%dt*(dijH%aa(p)-dij%aa(p))*(un(k,j)-un(k,i)) &
                  -mc_minus_ml%aa(p)*(du(k,j)-du(k,i))
          END DO
       END DO
    END DO
  END SUBROUTINE compute_fct_matrix_full

  SUBROUTINE compute_fct_matrix_approx(un,rkgal)
    USE mesh_handling
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np)             :: rkgal
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np)             :: X
    INTEGER :: i, j, k, p
    X =  rkgal
    DO i = 1, mesh%np
       DO k = 1 , k_dim+2
          rkgal(k,i) = rkgal(k,i) &
          + SUM(dijH%aa(mass%ia(i):mass%ia(i+1)-1)*(un(k,mass%ia(i):mass%ia(i+1)-1)-un(k,i)))
       END DO
    END DO
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          DO k = 1 , k_dim+2
             fctmat(k)%aa(p) = inputs%dt*(dijH%aa(p)-dij%aa(p))*(un(k,j)-un(k,i)) &
                  +mass%aa(p)*(-X(k,j)/lumped(j)+X(k,i)/lumped(i))
          END DO
       END DO
    END DO
  END SUBROUTINE compute_fct_matrix_approx

  SUBROUTINE update_fct_full(ulow,unext)
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: ulow
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: unext
    REAL(KIND=8), DIMENSION(k_dim+2) :: x
    INTEGER :: i, j, k, p
    DO i = 1, mesh%np
       x = 0.d0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          DO k = 1, k_dim+2
             x(k) = x(k) + lij%aa(p)*fctmat(k)%aa(p)
          END DO
       END DO
       DO k = 1, k_dim+2
          unext(k,i) = ulow(k,i) + x(k)/lumped(i)
       END DO
    END DO
  END SUBROUTINE update_fct_full

  SUBROUTINE update_fct(ulow,un,du,unext)
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: ulow, un, du
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: unext
    REAL(KIND=8), DIMENSION(k_dim+2) :: x
    INTEGER :: i, j, k, p
    DO i = 1, mesh%np
       x = 0.d0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          DO k = 1, k_dim+2
             x(k) = x(k) + lij%aa(p)*(inputs%dt*(dijH%aa(p)-dij%aa(p))*(un(k,j)-un(k,i))&
                  -mc_minus_ml%aa(p)*(du(k,j)-du(k,i)))
          END DO
       END DO
       DO k = 1, k_dim+2
          unext(k,i) = ulow(k,i) + x(k)/lumped(i)
       END DO
    END DO
  END SUBROUTINE update_fct

  SUBROUTINE compute_fct_matrix(un,du,fct)
    USE mesh_handling
    REAL(KIND=8), DIMENSION(mesh%np), INTENT(IN) :: un, du
    TYPE(matrice_bloc),               INTENT(IN) :: fct
    INTEGER :: i, j, p
    DO i = 1, mesh%np
       DO p = fct%ia(i), fct%ia(i+1) - 1
          j = fct%ja(p)
          fct%aa(p) = inputs%dt*(dijH%aa(p)-dij%aa(p))*(un(j)-un(i)) &
               -mc_minus_ml%aa(p)*(du(j)-du(i))
       END DO
    END DO
  END SUBROUTINE compute_fct_matrix

  SUBROUTINE compute_dij(un)
    USE mesh_handling
    USE boundary_conditions
    USE lambda_module
    USE lambda_module_full
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)  :: un
    INTEGER                                       :: i, p, j, k
    REAL(KIND=8)                                  :: norm_cij, lambda_max, lambdal, lambdar, pl, pr
    REAL(KIND=8)                                  :: ml, mr, ul, ur, el, er, rhol, rhor
    REAL(KIND=8), DIMENSION(k_dim)                :: nij


    !IF (inputs%type_test='RPE4') THEN
    !   
    !   RETURN
    !END IF

    
    DO i = 1, mesh%np
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          IF (i.NE.j) THEN
             DO k = 1, k_dim
                nij(k) = cij(k)%aa(p)
             END DO
             norm_cij = SQRT(SUM(nij**2))
             nij=nij/norm_cij
             mr = SUM(un(2:k_dim+1,j)*nij)
             ml = SUM(un(2:k_dim+1,i)*nij)
             rhor = un(1,j)
             rhol = un(1,i)
             er = un(k_dim+2,j) - 0.5d0*(SUM(un(2:k_dim+1,j)**2) - mr**2)/rhor
             el = un(k_dim+2,i) - 0.5d0*(SUM(un(2:k_dim+1,i)**2) - ml**2)/rhol
             ul = ml/rhol
             ur = mr/rhor
             pr = ABS(er-0.5d0*rhor*ur**2)*(gamma-1)
             pl = ABS(el-0.5d0*rhol*ul**2)*(gamma-1)
             CALL lambda(rhol,ul,pl,rhor,ur,pr,lambdal,lambdar)
             !CALL lambda_full(1.d-5,rhol,ul,pl,rhor,ur,pr,lambdal,lambdar,k)
             lambda_max = MAX(ABS(lambdal), ABS(lambdar))
             dij%aa(p) = norm_cij*lambda_max
          ELSE
             dij%aa(p) = 0.d0
          END IF
       END DO
    END DO
    DO i = 1, mesh%np
       dij%aa(diag(i)) = -SUM(dij%aa(dij%ia(i):dij%ia(i+1)-1))
    END DO
  END SUBROUTINE compute_dij

  SUBROUTINE compute_Roe_dij(un)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)  :: un
    INTEGER                                       :: i, p, j, k
    REAL(KIND=8)                                  :: norm_cij, lambda_max, pl, pr
    REAL(KIND=8)                                  :: ml, mr, ul, ur, el, er, rhol, rhor
    REAL(KIND=8), DIMENSION(k_dim)                :: nij
    REAL(KIND=8) :: hl, hr, thetal, thetar, ubar, hbar, cbar
    DO i = 1, mesh%np
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          IF (i.NE.j) THEN
             DO k = 1, k_dim
                nij(k) = cij(k)%aa(p)
             END DO
             norm_cij = SQRT(SUM(nij**2))
             nij=nij/norm_cij
             mr = SUM(un(2:k_dim+1,j)*nij)
             ml = SUM(un(2:k_dim+1,i)*nij)
             rhor = un(1,j)
             rhol = un(1,i)
             er = un(k_dim+2,j) - 0.5d0*(SUM(un(2:k_dim+1,j)**2) - mr**2)/rhor
             el = un(k_dim+2,i) - 0.5d0*(SUM(un(2:k_dim+1,i)**2) - ml**2)/rhol
             ul = ml/rhol
             ur = mr/rhor
             pr = ABS(er-0.5d0*rhor*ur**2)*(gamma-1)
             pl = ABS(el-0.5d0*rhol*ul**2)*(gamma-1)
             hl = (el+pl)/rhol
             hr = (er+pr)/rhor
             thetal = SQRT(rhol)/(SQRT(rhol)+SQRT(rhor))
             thetar = 1.d0 - thetal
             ubar = thetal*ul+thetar*ur
             hbar = thetal*hl+thetar*hr
             cbar = (gamma-1)*(hbar-0.5d0*ubar**2)
             lambda_max = MAX(ABS(ubar-cbar),ABS(ubar+cbar))
             dij%aa(p) = norm_cij*lambda_max
          ELSE
             dij%aa(p) = 0.d0
          END IF
       END DO
    END DO
    DO i = 1, mesh%np
       dij%aa(diag(i)) = -SUM(dij%aa(dij%ia(i):dij%ia(i+1)-1))
    END DO
  END SUBROUTINE compute_Roe_dij

  SUBROUTINE check_positivity(un)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np)                      :: u2, e
    INTEGER :: k
    REAL(KIND=8) emax, rhomax
    u2 = 0.d0
    DO k = 1, k_dim 
       u2 = u2 + (un(k+1,:)/un(1,:))**2
    END DO
    e = un(k_dim+2,:)/un(1,:) - u2/2
    rhomax = maxval(un(1,:))
    emax = maxval(e)
    IF (MINVAL(un(1,:))<-small*rhomax) THEN
       WRITE(*,*) ' Density is negative', MINVAL(un(1,:))
       STOP
    END IF
    IF (MINVAL(e)<-small*emax) THEN
       WRITE(*,*) ' Internal energy is negative', MINVAL(e)
       DO k = 1, mesh%np
          IF (e(k)<0.d0) THEN
             WRITE(*,*) ' location', mesh%rr(1:k_dim,k)
          END IF
       END DO
       STOP
    END IF
  END SUBROUTINE check_positivity
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
  
!!$  SUBROUTINE entropy_residual_bis(un,rkgal)
!!$    USE mesh_handling
!!$    USE input_data
!!$    USE boundary_conditions
!!$    USE sub_plot
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)  :: un, rkgal
!!$    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np)              :: DS
!!$    REAL(KIND=8), DIMENSION(mesh%np)  :: s, e, u2, res, en, press, pg, maxn, minn, rescale
!!$    INTEGER :: comp, k, i, j, p
!!$    REAL(KIND=8) :: zz
!!$
!!$    u2 = 0.d0
!!$    DO k = 1, k_dim 
!!$       u2 = u2 + (un(k+1,:)/un(1,:))**2
!!$    END DO
!!$    
!!$    IF (inputs%entropy_type=='p') THEN
!!$       press = (gamma-1.d0)*(un(k_dim+2,:) - un(1,:)*u2/2)
!!$       pg = (1/gamma)*press**(1.d0/gamma-1.d0)
!!$       en = gamma*press*pg  !S=Pressure*(1/gamma)
!!$       s = en/un(1,:)
!!$       DS(1,:) = pg*(gamma-1.d0)*u2/2
!!$       DO k = 1, k_dim 
!!$          DS(k+1,:) = -pg*(gamma-1.d0)*(un(k+1,:)/un(1,:))
!!$       END DO
!!$       DS(k_dim+2,:) = pg*(gamma-1.d0)
!!$    ELSE IF (inputs%entropy_type=='s') THEN
!!$       e = un(k_dim+2,:)/un(1,:) - u2/2
!!$       s = (1.d0/(gamma-1d0))*log(e) - log(un(1,:))
!!$       DS(1,:) = s + (1.d0/(gamma-1d0))*(u2/(2*e)-gamma)
!!$       DO k = 1, k_dim 
!!$          DS(k+1,:) = -un(k+1,:)/((gamma-1)*un(1,:)*e)
!!$       END DO
!!$       DS(k_dim+2,:) = (1.d0/(gamma-1d0))/e
!!$       en = un(1,:)*s  !S = rho*s
!!$    ELSE
!!$       WRITE(*,*) ' Bug: entropy_type not defined'
!!$       STOP
!!$    END IF
!!$
!!$    CALL maxmin(en,dij,maxn,minn)
!!$    !===Rescale
!!$    rescale = ABS(maxn-minn)/2
!!$    rescale = rescale + small*max(ABS(maxn),ABS(minn)) !Rescale 1
!!$    !rescale = max(ABS(maxn),ABS(minn)) !Rescale 2
!!$    !rescale = SQRT(rescale*max(ABS(maxn),ABS(minn))) !Rescale 3
!!$    !===Rescale
!!$    res = rkgal(1,:)*DS(1,:)
!!$    DO comp = 2, k_dim+2
!!$       res = res + rkgal(comp,:)*DS(comp,:)
!!$    END DO
!!$
!!$    DO i = 1, mesh%np
!!$       DO p = mass%ia(i), mass%ia(i+1) - 1
!!$          j = mass%ja(p)
!!$          zz = 0.d0
!!$          DO k = 1, k_dim
!!$             zz = zz + cij(k)%aa(p)*un(k+1,j)*s(j) !m*s=u*rho*s=u*S
!!$          END DO
!!$          res(i) = res(i) + zz
!!$       END DO
!!$    END DO
!!$    res = ABS(res)/rescale
!!$    !===Plot
!!$    IF (inputs%time+inputs%dt>inputs%Tfinal) THEN
!!$       !IF (k_dim==1) THEN
!!$       !   write(*,*) SUM(res)/mesh%np !mean (res*h/lumped)
!!$       !ELSE
!!$       !   write(*,*) SUM(res/SQRT(lumped))/mesh%np !mean (res*h/lumped)
!!$       !END IF
!!$       DO i = 1, mesh%np
!!$          s(i) = inputs%ce*res(i)*maxval(1.d0/dij%aa(mass%ia(i):mass%ia(i+1) - 1))
!!$       END DO
!!$       write(*,*) ' mean val |res/dij|', SUM(s)/mesh%np
!!$       SELECT CASE(k_dim)
!!$       CASE(1)
!!$          CALL plot_1d(mesh%rr(1,:),s,'res.plt')
!!$       CASE DEFAULT
!!$          CALL plot_scalar_field(mesh%jj, mesh%rr, s, 'res.plt')
!!$       END SELECT
!!$    END IF
!!$    !===Plot
!!$    
!!$    DO i = 1, mesh%np
!!$       DO p = resij%ia(i), resij%ia(i+1) - 1
!!$          j = resij%ja(p)
!!$          resij%aa(p) = max(res(i),res(j))
!!$       END DO
!!$    END DO
!!$    dijH%aa= MIN(dij%aa, inputs%ce*resij%aa)
!!$  END SUBROUTINE entropy_residual_bis
  
  SUBROUTINE entropy_residual(un,rkgal)
    USE mesh_handling
    USE input_data
    USE boundary_conditions
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)  :: un, rkgal
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np)              :: DS
    REAL(KIND=8), DIMENSION(k_dim+2,k_dim,mesh%np)        :: vv
    INTEGER,      DIMENSION(mesh%gauss%n_w)               :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w):: dw_loc
    REAL(KIND=8), DIMENSION(k_dim+2)  :: divful, dsl, yy, xx
    REAL(KIND=8), DIMENSION(mesh%np)  :: s, e, u2, res, en, press, pg, diffs, diffr, &
         maxn, minn, raxn, rinn, rescale, avgs, avgr, absres1, absres2, absres3, absres4, ccmin
    INTEGER :: comp, k, i, j, p, l, m, n, ps, pe
    REAL(KIND=8) :: ss, zz, spe, smin, smax, sbar, abszz

    u2 = 0.d0
    DO k = 1, k_dim 
       u2 = u2 + (un(k+1,:)/un(1,:))**2
    END DO

    IF (inputs%type_test=='RPE4') THEN
       e = un(1,:)*un(k_dim+2,:) - un(1,:)**2*u2/2 + 3*un(1,:)**3
       s =  (1.d0/(gamma-1d0))*LOG(e/un(1,:)**2) + LOG((3-un(1,:))/(3*un(1,:)))
       DS(1,:) = s + un(1,:)* ( (un(k_dim+2,:)+9*un(1,:)**2)/((gamma-1)*e) &
            -2/((gamma-1)*un(1,:)) -1/(3-un(1,:)) -1/un(1,:) )
       DO k = 1, k_dim
          DS(k+1,:) = un(1,:)*(-un(k+1,:)/e)/(gamma-1)
       END DO
       DS(k_dim+2,:) = un(1,:)*(un(1,:)/e)/(gamma-1)
       en = un(1,:)*s !S = rho*s
    ELSE IF (inputs%entropy_type=='p') THEN
       press = (gamma-1.d0)*(un(k_dim+2,:) - un(1,:)*u2/2)
       pg = (1/gamma)*press**(1.d0/gamma-1.d0)
       en = gamma*press*pg  !S=Pressure*(1/gamma)
       s = en/un(1,:)
       DS(1,:) = pg*(gamma-1.d0)*u2/2
       DO k = 1, k_dim 
          DS(k+1,:) = -pg*(gamma-1.d0)*(un(k+1,:)/un(1,:))
       END DO
       DS(k_dim+2,:) = pg*(gamma-1.d0)
    ELSE IF (inputs%entropy_type=='s') THEN
       e = un(k_dim+2,:)/un(1,:) - u2/2
       s = (1.d0/(gamma-1d0))*LOG(e) - LOG(un(1,:))
       DS(1,:) = s + (1.d0/(gamma-1d0))*(u2/(2*e)-gamma)
       DO k = 1, k_dim 
          DS(k+1,:) = -un(k+1,:)/((gamma-1)*un(1,:)*e)
       END DO
       DS(k_dim+2,:) = (1.d0/(gamma-1d0))/e
       en = un(1,:)*s  !S = rho*s
    ELSE
       WRITE(*,*) ' Bug: entropy_type not defined'
       STOP
    END IF

    !===Rescale
    !IF (inputs%type_test=='RPE4') THEN
    IF (.FALSE.) THEN   
       !sbar = (maxval(s)+minval(s))/2
!!$       CALL maxmin(un(1,:),mass,raxn,rinn)
!!$       CALL maxmin(s,dij,maxn,minn)
!!$       rescale = abs(maxn-minn)*(raxn+rinn)/2 +  abs(raxn-rinn)*abs(maxn-sbar)/2 + small*abs(maxn)      
       !CALL avg_diff(un(1,:),mass,diffr,avgr)
       !CALL avg_diff(s,mass,diffs,avgs)
       !rescale = diffs*avgr/2 + diffr*(ABS(avgs-sbar))/2 + small*avgr*abs(avgs) !BEST

       DO i = 1, mesh%np
          ps = mass%ia(i)
          pe = mass%ia(i+1)-1
          zz =  SUM(DS(1,mass%ja(ps:pe)))/(pe-ps+1)
          DS(1,i) = DS(1,i) - zz        !RIGHT SCALING WITH differences
          pg(i) = SUM(DS(:,i)*un(:,i))   
       END DO
       CALL maxmin(pg,dij,maxn,minn)
       CALL maxmin(s,dij,raxn,rinn)
       rescale = (ABS(maxn)+ABS(minn)) + small*(ABS(raxn)+ABS(rinn))    ! ABS(maxn-minn)/2 + small*maxn

    ELSE     
       IF (inputs%rescale==0) THEN !===(0)
          smax = MAXVAL(s)
          smin = MINVAL(s)
          sbar = (smax+smin)/2
          !IF (inputs%entropy_type=='s') THEN
          !   sbar = (smax+smin)/2
          !ELSE
          !   sbar = 0.d0
          !END IF
!!$          CALL maxmin(s,dij,maxn,minn)
!!$          CALL maxmin(un(1,:),dij,raxn,rinn)
!!$          rescale = (abs(maxn-minn)*raxn +  abs(raxn-rinn)*abs(maxn-sbar) + small*abs(maxn))/2
          CALL avg_diff(un(1,:),mass,diffr,avgr)
          CALL avg_diff(s,mass,diffs,avgs)
          rescale = diffs*avgr + diffr*(ABS(avgs-sbar)) + small*avgr*ABS(avgs) !BEST
       ELSE IF (inputs%rescale==1) THEN !===(1)
          smax = MAXVAL(s)
          smin = MINVAL(s)
          sbar = (smax+smin)/2
          DO i = 1, mesh%np
             ps = mass%ia(i)
             pe = mass%ia(i+1)-1
             !sbar =  SUM(DS(1,mass%ja(ps:pe)))/(pe-ps+1)
             pg(i) = (DS(1,i) - sbar)*un(1,i) + SUM(DS(2:,i)*un(2:,i))   !RIGHT SCALING WITH differences
          END DO
          CALL maxmin(pg,dij,maxn,minn)
          CALL maxmin(s,dij,raxn,rinn)
          rescale = (ABS(maxn)+ABS(minn)) + small*(ABS(raxn)+ABS(rinn))

          !==OLD
          CALL maxmin(en,dij,maxn,minn)
          rescale = ABS(maxn-minn)/2 + small*MAX(ABS(maxn),ABS(minn))

       ELSE IF(inputs%rescale==2) THEN !===(2)
          DO i = 1,mesh%np
             en(i) = un(1,i)
          END DO
          CALL maxmin(en,dij,maxn,minn)
          rescale = ABS(maxn-minn)/2 + small*maxn
       ELSE IF(inputs%rescale==3) THEN !===(3)  
          DO i = 1,mesh%np
             e(i) = un(k_dim+2,i)/un(1,i) - u2(i)/2
             en(i)=0.d0
             DO p = mass%ia(i), mass%ia(i+1) - 1
                j = mass%ja(p)
                yy(1) = un(1,j)
                yy(2:k_dim+1) = (-un(2:k_dim+1,i)*yy(1)/un(1,i) + un(2:k_dim+1,j))/un(1,i)
                yy(k_dim+2) = (-un(k_dim+2,i)*un(1,i)+SUM(un(2:k_dim+1,i)**2))*yy(1)/un(1,i)**3 &
                     - SUM(un(2:k_dim+1,i)*un(2:k_dim+1,j))/un(1,i)**2 + un(k_dim+2,j)/un(1,i)
                en(i) = en(i)+yy(1)**2/un(1,i)+SUM(yy(2:k_dim+1)**2)*un(1,i)/((gamma-1)*e(i)**2)&
                     + yy(k_dim+2)**2*un(1,i)/((gamma-1)*e(i))
             END DO
             en(i)=en(i)/(mass%ia(i+1)-mass%ia(i))
          END DO
          CALL maxmin(en,dij,maxn,minn)
          rescale = ABS(maxn-minn)/2 + small*maxn
       ELSE IF(inputs%rescale==4) THEN !===(4)  
          DO i = 1, mesh%np
             rescale(i)=0.d0
             DO p = mass%ia(i), mass%ia(i+1)-1
                j = mass%ja(p)
                rescale(i) = rescale(i) + ABS(SUM(DS(:,i)*(un(:,i)-un(:,j))))
             END DO
             rescale(i) = rescale(i)/(mass%ia(i+1)-mass%ia(i)) + small*ABS(un(1,i))
          END DO
       ELSE IF(inputs%rescale==5) THEN !===(5) 
          !===(5)
          DO i = 1, mesh%np
             e(i) = un(k_dim+2,i)/un(1,i) - u2(i)/2
             en(i)=0.d0
             DO p = mass%ia(i), mass%ia(i+1) - 1
                j = mass%ja(p)
                xx = un(:,i) - un(:,j)
                yy(1) = xx(1)
                yy(2:k_dim+1) = (-un(2:k_dim+1,i)*xx(1)/un(1,i) + xx(2:k_dim+1))/un(1,i)
                yy(k_dim+2) = (-un(k_dim+2,i)*un(1,i)+SUM(un(2:k_dim+1,i)**2))*xx(1)/un(1,i)**3 &
                     - SUM(un(2:k_dim+1,i)*xx(2:k_dim+1))/un(1,i)**2 + xx(k_dim+2)/un(1,i)
                en(i) = en(i)+yy(1)**2/un(1,i)+SUM(yy(2:k_dim+1)**2)*un(1,i)/((gamma-1)*e(i)**2)&
                     + yy(k_dim+2)**2*un(1,i)/((gamma-1)*e(i))
             END DO
             en(i)=en(i)/(mass%ia(i+1)-mass%ia(i))
             rescale(i) = SQRT(ABS(SUM(DS(:,i)*un(:,i)))*ABS(en(i)))/2 + small*ABS(en(i))
          END DO
       END IF
    END IF
    !===Rescale

    IF(inputs%rescale==6) THEN !===(6)  
       !================TEST IN PAPER DONE WITH THIS SETTING
       !================DO NOT CHANGE
       avgs = s
       !===Here I am testing
       !avgs = sum(s)/mesh%np !TO BE COMMENTED FOR TEST IN EULER FCT PAPER
       !===Here I am testing
       res = rkgal(1,:)*(DS(1,:)- avgs)
       absres1 = ABS(res)
       DO comp = 2, k_dim+2
          res = res + rkgal(comp,:)*DS(comp,:)
          absres1 =  absres1 + ABS(rkgal(comp,:)*DS(comp,:))
       END DO
       absres3=ABS(res)
       !===It is essential to take the absolute value on each components
       !===to get correct convergence in 1D on the expansion wave.
       !===This also gives the correct scaling for VdW in 1D.

       DO i = 1, mesh%np
          ss = 0.d0
          zz = 0.d0
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = mass%ja(p)
             DO k = 1, k_dim
                zz = zz + cij(k)%aa(p)*un(k+1,j)*(s(j)-avgs(i))
                ss = ss + cij(k)%aa(p)*un(k+1,j)*(s(j))
             END DO
          END DO
          res(i) = ABS(res(i) + zz)
          absres2(i) = ABS(zz) !==Essential to have this normalization in 2D
          absres4(i) = abs(ss)
       END DO

       !================TEST IN PAPER DONE WITH THIS SETTING
       !================DO NOT CHANGE
       !s = MIN(1.d0,inputs%ce*res/(absres1+absres2+small*MAXVAL(res)))
       !================TEST IN PAPER DONE WITH THIS SETTING
       !================DO NOT CHANGE
       !===Here I am testing TO BE COMMENTED FOR TESTS IN EULER FCT PAPER
       !u2 = 0.d0
       !DO k = 1, k_dim 
       !   u2 = u2 + (un(k+1,:)/un(1,:))**2
       !END DO
       !e = un(k_dim+2,:)/un(1,:) - u2/2
       !s = (1.d0/(gamma-1d0))*LOG(e) - LOG(un(1,:))
       !DO i = 1, mesh%np
       !   write(*,*) absres1(i),absres2(i), absres3(i), absres4(i), s(i)!-avgs(i)
       !end DO
       s = max(res-1.d-8*absres4,0.d0) 
       s = MIN(1.d0,inputs%ce*res/max(absres1+absres2,small*MAXVAL(res)))
       !===Here I am testing
       DO i = 1, mesh%np
          DO p = resij%ia(i), resij%ia(i+1) - 1
             j = resij%ja(p)
             dijH%aa(p) = dij%aa(p)*((s(i)+s(j))/2)
          END DO
       END DO
  
       !================TEST IN PAPER DONE WITH THIS SETTING.
       !================DO NOT CHANGE
       IF (inputs%time+inputs%dt.GE.inputs%Tfinal) THEN
          SELECT CASE(k_dim)
          CASE(1)
             CALL plot_1d(mesh%rr(1,:),s,'res.plt')
          CASE DEFAULT
             CALL plot_scalar_field(mesh%jj, mesh%rr, s, 'res.plt')
          END SELECT
       END IF
       RETURN    
    END IF
    !===(6)




    
    !===Compute residual
    res = rkgal(1,:)*DS(1,:)
    DO comp = 2, k_dim+2
       res = res + rkgal(comp,:)*DS(comp,:)
    END DO

!!$    res = 0.d0
!!$    DO comp = 1, k_dim+2
!!$       vv(comp,:,:)=flux(comp,un)
!!$    END DO
!!$    DO m = 1, mesh%me
!!$       j_loc = mesh%jj(:,m)
!!$       DO l = 1, mesh%gauss%l_G
!!$          dw_loc = mesh%gauss%dw(:,:,l,m)
!!$          divful = 0.d0
!!$          dsl = 0.d0
!!$          DO comp = 1, k_dim+2
!!$             dsl(comp) = dsl(comp) + SUM(DS(comp,j_loc)*mesh%gauss%ww(:,l))
!!$             DO k = 1, k_dim
!!$                divful(comp) =  divful(comp) + SUM(vv(comp,k,j_loc)*dw_loc(k,:))
!!$             END DO
!!$          END DO
!!$          DO n = 1, mesh%gauss%n_w
!!$             i = j_loc(n)
!!$             res(i) = res(i) - SUM(dsl*divful)*mesh%gauss%ww(n,l)*mesh%gauss%rj(l,m)
!!$          END DO
!!$       END DO
!!$    END DO
    !===End compute residual

    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          zz = 0.d0
          DO k = 1, k_dim
             zz = zz + cij(k)%aa(p)*un(k+1,j)*s(j) !m*s=u*rho*s=u*S
          END DO
          res(i) = res(i) + zz
       END DO
    END DO
    res = ABS(res)/rescale

    !===Plot
    IF (inputs%time+inputs%dt.GE.inputs%Tfinal) THEN
       DO i = 1, mesh%np
          s(i) = inputs%ce*res(i)*MAXVAL(1.d0/dij%aa(mass%ia(i):mass%ia(i+1) - 1))
       END DO
       !write(*,*) ' mean val |res/dij|', SUM(s)/mesh%np
       SELECT CASE(k_dim)
       CASE(1)
          CALL plot_1d(mesh%rr(1,:),s,'res.plt')
          CALL plot_1d(mesh%rr(1,:),rescale,'rescale.plt')
          DO i = 1, mesh%np
             s(i) = MAXVAL(dij%aa(mass%ia(i):mass%ia(i+1) - 1))
          END DO
          CALL plot_1d(mesh%rr(1,:),s,'dij.plt')
          CALL plot_1d(mesh%rr(1,:),res,'residual.plt')
       CASE DEFAULT
          CALL plot_scalar_field(mesh%jj, mesh%rr, s, 'res.plt')
          CALL plot_scalar_field(mesh%jj, mesh%rr, rescale, 'rescale.plt')
       END SELECT
    END IF
    !===Plot

    DO i = 1, mesh%np
       s(i) = SUM(res(resij%ja(resij%ia(i):resij%ia(i+1) - 1)))/(resij%ia(i+1)-resij%ia(i))
    END DO
    DO i = 1, mesh%np
       DO p = resij%ia(i), resij%ia(i+1) - 1
          j = resij%ja(p)
          resij%aa(p) = (s(i)+s(j))/2
       END DO
    END DO

!!$    DO i = 1, mesh%np
!!$       DO p = resij%ia(i), resij%ia(i+1) - 1
!!$          j = resij%ja(p)
!!$          resij%aa(p) = max(res(i),res(j))
!!$       END DO
!!$    END DO

    dijH%aa= MIN(dij%aa, inputs%ce*resij%aa)
  END SUBROUTINE entropy_residual

  SUBROUTINE smoothness_viscosity(un)
    USE mesh_handling
    USE boundary_conditions
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np) :: alpha, bbeta, u2, e, s, ent, press, pg, num, ddenom, denom
    REAL(KIND=8) :: nmax, nmin, ddmax, ddmin
    INTEGER      :: i, j, k, p

    u2 = 0.d0
    DO k = 1, k_dim 
       u2 = u2 + (un(k+1,:)/un(1,:))**2
    END DO

    IF (inputs%type_test=='RPE4') THEN
       e = un(1,:)*un(k_dim+2,:) - un(1,:)**2*u2/2 + 3*un(1,:)**3
       s =  (1.d0/(gamma-1d0))*LOG(e/un(1,:)**2) + LOG((3-un(1,:))/(3*un(1,:)))
       ent = un(1,:)*s !S = rho*s
    ELSE IF (inputs%entropy_type=='p') THEN
       press = (gamma-1.d0)*(un(k_dim+2,:) - un(1,:)*u2/2)
       pg = (1/gamma)*press**(1.d0/gamma-1.d0)
       ent = gamma*press*pg  !S=Pressure*(1/gamma)
    ELSE IF (inputs%entropy_type=='s') THEN
       e = un(k_dim+2,:)/un(1,:) - u2/2
       s = (1.d0/(gamma-1d0))*LOG(e) - LOG(un(1,:))
       ent = un(1,:)*s  !S = rho*s
    ELSE
       WRITE(*,*) ' Bug: entropy_type not defined'
       STOP
    END IF
    
    ent = un(k_dim+2,:)
    DO i = 1, mesh%np
       denom(i) =0.d0
       ddenom(i) =0.d0
       num(i)  = 0.d0
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          IF (i==j) CYCLE
          !num(i) = num(i) + ent(j) - ent(i)
          !denom(i) = denom(i) + ABS(ent(j) - ent(i))
          !ddenom(i) = ddenom(i) + ABS(ent(j)) + ABS(ent(i))
          num(i) = num(i) + stiff%aa(p)*(ent(j) - ent(i))
          denom(i) = denom(i) + ABS(stiff%aa(p)*(ent(j) - ent(i)))
          ddenom(i) = ddenom(i) + ABS(stiff%aa(p))*(ABS(ent(j)) + ABS(ent(i)))
       END DO
       num(i) = num(i)/lumped(i)
       denom(i) = denom(i)/lumped(i)
       ddenom(i) = ddenom(i)/lumped(i)
       IF (denom(i).GE.1.d-7*ddenom(i)) THEN
          alpha(i) = ABS(num(i))/(denom(i))
          bbeta(i) = ABS(num(i))/(ddenom(i))
       ELSE
          alpha(i) = 0.d0
          bbeta(i) = 0.d0
          num(i) = small*ddenom(i)
       END IF
       alpha(i) = alpha(i)**2
       !alpha(i) = (MAX(alpha(i)-0.25,0.d0)*4.d0/3.d0)**2 !New formulation, June 16, 2017
    END DO
 
!!$    DO i = 1, mesh%np
!!$       nmax = MAXVAL(num(dij%ja(dij%ia(i):dij%ia(i+1)-1)))
!!$       nmin = MINVAL(num(dij%ja(dij%ia(i):dij%ia(i+1)-1)))
!!$       ddmax = MAXVAL(ddenom(dij%ja(dij%ia(i):dij%ia(i+1)-1)))
!!$       ddmin = MINVAL(ddenom(dij%ja(dij%ia(i):dij%ia(i+1)-1)))
!!$       IF (denom(i).GE.1d-7*ddenom(i)) THEN
!!$          bbeta(i) = ABS(nmax-nmin)/(ABS(nmax)+ABS(nmin))
!!$          !bbeta(i) = ABS(nmax-nmin)/(2*ABS(ent(i)))
!!$          !bbeta(i) = ABS(nmax-nmin)/(ddenom(i))
!!$          !bbeta(i) = ABS(nmax-nmin)/(ddmax+ddmin)
!!$       ELSE
!!$          bbeta(i) = ABS(nmax-nmin)/(ddmax+ddmin)
!!$       END IF
!!$    END DO
    DO i = 1, mesh%np
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          dijH%aa(p)= dij%aa(p)*MAX(MIN(alpha(i),bbeta(i)),MIN(alpha(j),bbeta(j)))
       END DO
    END DO

    !===Plot
    IF (inputs%time+inputs%dt.GE.inputs%Tfinal) THEN
       SELECT CASE(k_dim)
       CASE(1)
          CALL plot_1d(mesh%rr(1,:),bbeta,'beta.plt')
       CASE DEFAULT
          CALL plot_scalar_field(mesh%jj, mesh%rr, bbeta, 'beta.plt')
       END SELECT
    END IF
    !===Plot
    
  END SUBROUTINE smoothness_viscosity


  SUBROUTINE maxmin(un,mat,maxn,minn)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: un
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: maxn, minn
    INTEGER      :: i
    DO i = 1, SIZE(un)
       maxn(i) = MAXVAL(un(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
       minn(i) = MINVAL(un(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
    END DO
  END SUBROUTINE maxmin

  SUBROUTINE avg_diff(un,mat,diff,avg)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: un
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: diff, avg
    INTEGER      :: i
    DO i = 1, SIZE(un)
       diff(i) = SUM(ABS(un(mat%ja(mat%ia(i):mat%ia(i+1)-1))-un(i)))/(mat%ia(i+1)-mat%ia(i)-1)
       avg(i) = SUM(ABS(un(mat%ja(mat%ia(i):mat%ia(i+1)-1))))/(mat%ia(i+1)-mat%ia(i)-1)
    END DO
  END SUBROUTINE avg_diff


  SUBROUTINE LOCAL_limit(unext,maxn,minn,mat,mass,comp)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: mass, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)   :: mat
    INTEGER,                    INTENT(IN)   :: comp
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np) :: unext
    LOGICAL,      DIMENSION(mesh%np)         :: check
    REAL(KIND=8), PARAMETER :: smallplus = 1.d-15
    REAL(KIND=8) :: maxni, minni, ui, uij, xij, lambdai, umax, usmall
    INTEGER      :: i, j, p

    !===Compute lij
    umax = MAX(MAXVAL(maxn),-MINVAL(minn))
    usmall = umax*smallplus
    !check=.TRUE.
    !IF (comp==1) check(rho_js_D) = .FALSE.
    lij%aa = 1.d0
    DO i = 1, mesh%np
       lambdai = 1.d0/(mat%ia(i+1) - 1.d0 - mat%ia(i))
       maxni = maxn(i)
       minni = minn(i)
       ui = unext(comp,i)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          xij = mat%aa(p)/(mass(i)*lambdai)
          uij = ui + xij
          IF (uij>maxni) THEN
             lij%aa(p) = MIN(ABS(maxni - ui)/(ABS(xij)+usmall),1.0)
          ELSE IF (uij<minni) THEN
             lij%aa(p) = MIN(ABS(minni - ui)/(ABS(xij)+usmall),1.d0)
          END IF
       END DO
    END DO
    lij%aa(diag) = 0.d0
  END SUBROUTINE LOCAL_limit

  SUBROUTINE transpose_op(mat,TYPE)
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(INOUT):: mat
    CHARACTER(LEN=3),  INTENT(IN)   :: TYPE
    INTEGER, DIMENSION(SIZE(mat%ia)) :: iao
    INTEGER:: i, j, p, next
    IF  (TYPE/='min' .AND. TYPE/='max') THEN
       WRITE(*,*) ' BUG in tanspose_op'
       STOP
    END IF
    iao = mat%ia
    DO i = 1, SIZE(mat%ia)-1
       DO p = mat%ia(i), mat%ia(i+1)-1 
          j = mat%ja(p)
          next = iao(j)
          iao(j) = next+1
          IF (j.LE.i) CYCLE
          IF (TYPE=='min') THEN
             mat%aa(next) = MIN(mat%aa(p),mat%aa(next))
             mat%aa(p) = mat%aa(next)
          ELSE IF (TYPE=='max') THEN
             mat%aa(next) = MAX(mat%aa(p),mat%aa(next))
             mat%aa(p) = mat%aa(next)
          END IF
       END DO
    END DO
  END SUBROUTINE transpose_op


  !===JUNK
  SUBROUTINE FCT_generic(ulow,maxn,minn,mat,dg,comp)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: dg, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    INTEGER,                    INTENT(IN)  :: comp
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np):: ulow
    REAL(KIND=8), DIMENSION(mesh%np)        :: Qplus, Qminus, Pplus, Pminus, Rplus, Rminus
    REAL(KIND=8), PARAMETER :: smallplus = small, smallminus = -small
    REAL(KIND=8) :: fij
    INTEGER      :: i, j, p, jp, jm
    Qplus  = dg*(maxn-ulow(comp,:))
    Qminus = dg*(minn-ulow(comp,:))
    Pplus  = smallplus
    Pminus = smallminus
    DO i = 1, mesh%np

       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          jp = 0
          jm =0
          IF (fij.GE.0.d0) THEN
             jp = jp + 1
             Pplus(i)  = Pplus(i) + fij
          ELSE
             jm = jm + 1
             Pminus(i) = Pminus(i) + fij
          END IF
       END DO
       IF (jp>0) THEN
          Rplus(i)  =  MIN(Qplus(i)/Pplus(i),1.d0)
       ELSE
          RPLUS(i) = 1.d0
       END IF
       IF (jm>0) THEN
          Rminus(i) =  MIN(Qminus(i)/Pminus(i),1.d0)
       ELSE
          Rminus(i) =1.d0
       END IF
    END DO

    DO i = 1, mesh%np
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             lij%aa(p) = MIN(Rplus(i),Rminus(j))
          ELSE
             lij%aa(p) = MIN(Rminus(i),Rplus(j))
          END IF
       END DO
    END DO
  END SUBROUTINE FCT_generic


  REAL FUNCTION  Median(X, N)
    IMPLICIT  NONE
    REAL(KIND=8), DIMENSION(1:), INTENT(IN) :: X
    INTEGER, INTENT(IN)                :: N
    REAL(KIND=8), DIMENSION(1:N)            :: Temp
    INTEGER                            :: i

    DO i = 1, N                       ! make a copy
       Temp(i) = X(i)
    END DO
    CALL  Sort(Temp, N)               ! sort the copy
    IF (MOD(N,2) == 0) THEN           ! compute the median
       Median = (Temp(N/2) + Temp(N/2+1)) / 2.0
    ELSE
       Median = Temp(N/2+1)
    END IF
  END FUNCTION  Median
  ! --------------------------------------------------------------------
  ! INTEGER FUNCTION  FindMinimum():
  !    This function returns the location of the minimum in the section
  ! between Start and End.
  ! --------------------------------------------------------------------

  INTEGER FUNCTION  FindMinimum(x, Start, END)
    IMPLICIT  NONE
    REAL(KIND=8), DIMENSION(1:), INTENT(IN) :: x
    INTEGER, INTENT(IN)                :: Start, END
    REAL(KIND=8)                       :: Minimum
    INTEGER                            :: Location
    INTEGER                            :: i

    Minimum  = x(Start)		! assume the first is the min
    Location = Start			! record its position
    DO i = Start+1, END		! start with next elements
       IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
          Minimum  = x(i)		!      Yes, a new minimum found
          Location = i                !      record its position
       END IF
    END DO
    FindMinimum = Location        	! return the position
  END FUNCTION  FindMinimum

  ! --------------------------------------------------------------------
  ! SUBROUTINE  Swap():
  !    This subroutine swaps the values of its two formal arguments.
  ! --------------------------------------------------------------------

  SUBROUTINE  Swap(a, b)
    IMPLICIT  NONE
    REAL(KIND=8), INTENT(INOUT) :: a, b
    REAL(KIND=8)                :: Temp

    Temp = a
    a    = b
    b    = Temp
  END SUBROUTINE  Swap

  ! --------------------------------------------------------------------
  ! SUBROUTINE  Sort():
  !    This subroutine receives an array x() and sorts it into ascending
  ! order.
  ! --------------------------------------------------------------------

  SUBROUTINE  Sort(x, Size)
    IMPLICIT  NONE
    REAL(KIND=8), DIMENSION(1:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN)                   :: Size
    INTEGER                               :: i
    INTEGER                               :: Location

    DO i = 1, Size-1			! except for the last
       Location = FindMinimum(x, i, Size)	! find min from this to last
       CALL  Swap(x(i), x(Location))	! swap this and the minimum
    END DO
  END SUBROUTINE  Sort



 
END MODULE update
