MODULE limiting_euler
  USE input_data
  USE mesh_handling
  USE matrix_type
  PUBLIC :: compute_fct_matrix_full
  REAL(KIND=8), PARAMETER, PRIVATE :: small=1.d-8
CONTAINS

  SUBROUTINE compute_fct_matrix_full(scalar,un,du,fctmat,dijH,dijL,mc_minus_ml)
    IMPLICIT NONE
    TYPE(matrice_bloc)                           :: dijL, dijH, mc_minus_ml
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size)       :: fctmat
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: un, du
    REAL(KIND=8) :: scalar
    INTEGER :: i, j, k, p
    DO i = 1, mesh%np
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          DO k = 1, inputs%syst_size
             fctmat(k)%aa(p) = scalar*(dijH%aa(p)-dijL%aa(p))*(un(k,j)-un(k,i)) &
                  -mc_minus_ml%aa(p)*(du(k,j)-du(k,i))
          END DO
       END DO
    END DO
  END SUBROUTINE compute_fct_matrix_full

  SUBROUTINE update_fct_full(ulow,unext,fctmat,lij,lumped)
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np) :: ulow
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np) :: unext
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size)   :: fctmat
    TYPE(matrice_bloc)                                :: lij
    REAL(KIND=8), DIMENSION(mesh%np)                  :: lumped
    REAL(KIND=8), DIMENSION(inputs%syst_size) :: x
    INTEGER :: i, k, p
    DO i = 1, mesh%np
       x = 0.d0
       DO p = lij%ia(i), lij%ia(i+1) - 1
          DO k = 1, inputs%syst_size
             x(k) = x(k) + lij%aa(p)*fctmat(k)%aa(p)
          END DO
       END DO
       DO k = 1, inputs%syst_size
          unext(k,i) = ulow(k,i) + x(k)/lumped(i)
       END DO
    END DO
  END SUBROUTINE update_fct_full

   SUBROUTINE LOCAL_limit(unext,maxn,minn,mat,mass,comp,diag,lij)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: mass, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)   :: mat
    INTEGER,                    INTENT(IN)   :: comp
    INTEGER, DIMENSION(:),     INTENT(IN)    :: diag
    TYPE(matrice_bloc),         INTENT(OUT)  :: lij
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np) :: unext
    REAL(KIND=8), PARAMETER :: smallplus = 1.d-15
    REAL(KIND=8) :: maxni, minni, ui, uij, xij, lambdai, umax, usmall
    INTEGER      :: i, j, p

    !===Compute lij
    umax = MAX(MAXVAL(maxn),-MINVAL(minn))
    usmall = umax*smallplus
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
  
 SUBROUTINE FCT_generic(ulow,maxn,minn,mat,dg,comp,lij)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: dg, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    INTEGER,                    INTENT(IN)  :: comp
    TYPE(matrice_bloc),         INTENT(OUT)  :: lij
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np):: ulow
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
  
  
  SUBROUTINE limit_specific_entropy(ulow,un,cmin,lumped,fctmat,lij)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: ulow, un
    REAL(KIND=8), DIMENSION(mesh%np)                     :: cmin
    REAL(KIND=8), DIMENSION(mesh%np)                     :: lumped
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size)      :: fctmat
    TYPE(matrice_bloc),                       INTENT(OUT):: lij
    REAL(KIND=8), DIMENSION(mesh%np)  :: Esmall
    REAL(KIND=8), DIMENSION(inputs%syst_size)  :: ul, ur, Pij
    REAL(KIND=8) :: lambdai, coeff, psir, psil, ll, lr, llold, lrold, Budget
    INTEGER      :: i, j, p, k

    !===
    DO i = 1, mesh%np
       Esmall(i)= small*MINVAL(un(inputs%syst_size,lij%ja(lij%ia(i):lij%ia(i+1)-1)))
       lambdai = 1.d0/(lij%ia(i+1) - 1 - lij%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))

       Budget =0.d0
       !===End Budget

       DO p = lij%ia(i), lij%ia(i+1) - 1
          j =  lij%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p)
          DO k = 1, inputs%syst_size
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(k,i) + lr*Pij(k) !===Density must be positive
          END DO
          if (ur(1)<0.d0) then
             lij%aa(p) = 0.d0  !===CFL is too large
             CYCLE
          end if
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
      REAL(KIND=8), DIMENSION(inputs%syst_size), INTENT(IN) :: u
      REAL(KIND=8),                     INTENT(IN) :: cmin
      REAL(KIND=8)                                 :: psi, Budget
      psi = u(inputs%syst_size) - SUM(u(2:inputs%syst_size-1)**2)/(2.d0*u(1)) - cmin*u(1)**gamma - Budget
    END FUNCTION psi_func
    FUNCTION psi_prime_func(Pij,u,cmin) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(inputs%syst_size), INTENT(IN) :: u, Pij
      REAL(KIND=8),                     INTENT(IN) :: cmin
      REAL(KIND=8)                                 :: psi
      psi = Pij(inputs%syst_size) - SUM(u(2:inputs%syst_size-1)*Pij(2:inputs%syst_size-1))/u(1) &
           + Pij(1)*SUM(u(2:inputs%syst_size-1)**2)/(2*u(1)**2) &
           - cmin*gamma*Pij(1)*u(1)**(gamma-1.d0)
    END FUNCTION psi_prime_func
  END SUBROUTINE limit_specific_entropy

  SUBROUTINE relax(un,minn,maxn,stiff,urelaxi,drelaxi)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:)              :: un
    REAL(KIND=8), DIMENSION(:)              :: minn
    REAL(KIND=8), DIMENSION(:)              :: maxn
    TYPE(matrice_bloc)                      :: stiff
    REAL(KIND=8), DIMENSION(:)              :: urelaxi, drelaxi
    REAL(KIND=8), DIMENSION(SIZE(un))       :: alpha, denom
    INTEGER      :: i, j, p
    REAL(KIND=8) :: norm
    alpha = 0.d0
    DO i = 1, mesh%np
       norm = 0.d0
       DO p = stiff%ia(i), stiff%ia(i+1) - 1
          j = stiff%ja(p)
          IF (i==j) CYCLE
          alpha(i) = alpha(i) + stiff%aa(p)*(un(i) - un(j))
          norm = norm + stiff%aa(p)
       END DO
       alpha(i) = alpha(i)/norm
    END DO
    SELECT CASE(inputs%limiter_type)
    CASE('avg') !==Average
       denom = 0.d0
       DO i = 1, SIZE(un)
          DO p = stiff%ia(i), stiff%ia(i+1) - 1
             j = stiff%ja(p)
             IF (i==j) CYCLE
             denom(i) = denom(i) + alpha(j) + alpha(i)
          END DO
       END DO
       DO i = 1, SIZE(un)
          denom(i) = denom(i)/(2*(stiff%ia(i+1)-stiff%ia(i)-1))
       END DO
       maxn = MIN(urelaxi*maxn,maxn + ABS(denom)/2)
       minn = MAX(drelaxi*minn,minn - ABS(denom)/2)
    CASE ('minmod') !===Minmod
       denom = alpha    
       DO i = 1, SIZE(un)
          DO p = stiff%ia(i), stiff%ia(i+1) - 1
             j = stiff%ja(p)
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
    CASE DEFAULT
       WRITE(*,*) ' BUG in relax'
       STOP
    END SELECT
  END SUBROUTINE RELAX

  SUBROUTINE relax_cmin(un,cmin,mass,drelaxi)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:)            :: un
    TYPE(matrice_bloc)                      :: mass
    REAL(KIND=8), DIMENSION(:)              :: drelaxi
    REAL(KIND=8), DIMENSION(:)              :: cmin
    REAL(KIND=8), DIMENSION(SIZE(cmin))     :: dc
    REAL(KIND=8), DIMENSION(inputs%syst_size):: ul
    INTEGER      :: i, j, p
    REAL(KIND=8) :: cl
    dc = 0.d0
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          ul(:) = (un(:,i)+un(:,j))/2
          cl = (ul(inputs%syst_size)-SUM(ul(2:inputs%syst_size-1)**2)/(2.d0*ul(1)))/(ul(1)**gamma)
          dc(i) = MAX(dc(i),cl-cmin(i))
       END DO
    END DO
    cmin = MAX(drelaxi*cmin, cmin - dc)
  END SUBROUTINE RELAX_cmin
END MODULE limiting_euler
