MODULE limiters
  USE matrix_type
  USE input_data
  IMPLICIT NONE
  PUBLIC ::  LOCAL_limit, FCT_generic, fct_visc, maxmin, relax, compute_relax_coeff
  TYPE(matrice_bloc), PUBLIC :: lij, fctmat
  PRIVATE
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC  :: up_relax, down_relax
  REAL(KIND=8), PARAMETER                      :: up_relax_hard=1.001, down_relax_hard=.999d0   
CONTAINS
  SUBROUTINE LOCAL_limit(diag,lumped,ulow,unext,maxn,minn)
    USE mesh_handling
    USE CSR_transpose
    IMPLICIT NONE
    INTEGER,      DIMENSION(:), INTENT(IN)   :: diag
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: lumped
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: maxn, minn
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: ulow
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)  :: unext
    REAL(KIND=8) :: x, maxni, minni, ui, uij, pij, lambdai, usmall, small=0.d-7
    INTEGER      :: i, j, p,  ps, pe

    !===Compute lij
    usmall = small*MAXVAL(ABS(ulow))
    DO i = 1, mesh%np -1
       !nblambda = (mass%ia(i+1) - 1.d0 - mass%ia(i))
       !lambdai = lumped(i)/nblambda
       ps = fctmat%ia(i)
       pe = fctmat%ia(i+1) - 1
       x = 0.d0
       DO p = ps, pe
          x  = x + ABS(fctmat%aa(p))
       END DO
       maxni = maxn(i)
       minni = minn(i)
       ui = ulow(i)
       DO p = ps, pe
          j = fctmat%ja(p)
          IF(i==j) CYCLE
          lambdai = lumped(i)*ABS(fctmat%aa(p))/x !Best choice
          pij = fctmat%aa(p)/lambdai
          uij = ui + pij
          IF (uij>maxni) THEN
             lij%aa(p) = MIN(ABS(maxni - ui)/(ABS(pij)+usmall),1.d0)
          ELSE IF (uij<minni) THEN
             lij%aa(p) = MIN(ABS(minni - ui)/(ABS(pij)+usmall),1.d0)
          ELSE
             lij%aa(p) =1.d0
          END IF
       END DO
    END DO
    lij%aa(diag) = 0.d0
    !===lij = min(lij,lij^T)
    CALL transpose_op(lij,'min')

    !===Update u
    DO i = 1, mesh%np - 1 
       x = 0.d0
       DO p = fctmat%ia(i), fctmat%ia(i+1) - 1
          x = x + lij%aa(p)*fctmat%aa(p)
       END DO
       unext(i) = ulow(i) + x/lumped(i)
    END DO

  END SUBROUTINE LOCAL_limit

  SUBROUTINE FCT_generic(diag,ulow,unext,maxn,minn,mat,dg)
    IMPLICIT NONE
    INTEGER,      DIMENSION(:), INTENT(IN)   :: diag
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: dg, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ulow
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(SIZE(ulow))       :: Qplus, Qminus, Pplus, Pminus, Rplus, Rminus
    REAL(KIND=8), PARAMETER :: smallplus = 1.d-15, smallminus = -1.d-15
    REAL(KIND=8) :: x, fij, lij_loc
    INTEGER      :: i, j, p
    Qplus  = dg*(maxn-ulow)
    Qminus = dg*(minn-ulow)
    Pplus  = smallplus
    Pminus = smallminus
    DO i = 1, SIZE(ulow)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             Pplus(i)  = Pplus(i) + fij
          ELSE
             Pminus(i) = Pminus(i) + fij
          END IF
       END DO
    END DO
    Rplus  =  MIN(Qplus,Pplus)/Pplus
    Rminus =  MIN(ABS(Qminus),ABS(Pminus))/ABS(Pminus)
    DO i = 1, SIZE(ulow)
       x = 0.d0
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             lij_loc = MIN(Rplus(i),Rminus(j))
          ELSE
             lij_loc = MIN(Rminus(i),Rplus(j))
          END IF
          x = x + lij_loc*fij
          lij%aa(p) = lij_loc
       END DO
       unext(i) = ulow(i) + x/dg(i)
    END DO
    lij%aa(diag) = 0.d0
  END SUBROUTINE FCT_generic

  SUBROUTINE FCT_visc(diag,unext,un,maxn,minn,mat,dg,nonconvex)
    IMPLICIT NONE
    INTEGER,      DIMENSION(:), INTENT(IN)   :: diag
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: un, dg, maxn, minn
    TYPE(matrice_bloc),         INTENT(INOUT):: mat, nonconvex
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: unext
    REAL(KIND=8), DIMENSION(SIZE(un))       :: Qplus, Qminus, Pplus, Pminus, Rplus, Rminus
    REAL(KIND=8), PARAMETER :: smallplus = 1.d-15, smallminus = -1.d-15
    REAL(KIND=8) :: x, fij, lij
    INTEGER      :: i, j, p 
    Qplus  = dg*(maxn-unext)
    Qminus = dg*(minn-unext)
    Pplus  = smallplus
    Pminus = smallminus
    DO i = 1, SIZE(un)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = -mat%aa(p)*(un(j)-un(i))
          IF (fij.GE.0.d0) THEN
             Pplus(i)  = Pplus(i) + fij
          ELSE
             Pminus(i) = Pminus(i) + fij
          END IF
       END DO
       Rplus(i)  =  MIN(Qplus(i)/Pplus(i),1.d0)
       Rminus(i) =  MIN(Qminus(i)/Pminus(i),1.d0)
    END DO
    DO i = 1, SIZE(un)
       x = 0.d0
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = -mat%aa(p)*(un(j)-un(i))
          IF (fij.GE.0.d0) THEN
             lij = MIN(Rplus(i),Rminus(j))
          ELSE
             lij = MIN(Rminus(i),Rplus(j))
          END IF
          mat%aa(p) = MAX((1.d0-lij),nonconvex%aa(p))*mat%aa(p)
       END DO
    END DO
    DO i = 1, SIZE(un)
       mat%aa(diag(i)) = -SUM(mat%aa(mat%ia(i):mat%ia(i+1)-1))
    END DO
  END SUBROUTINE FCT_visc

  SUBROUTINE maxmin(l_max,un,mat,maxn,minn)
    IMPLICIT NONE
    INTEGER                                 :: l_max
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: un
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: maxn, minn
    REAL(KIND=8), DIMENSION(SIZE(un))       :: max_in, min_in
    INTEGER      :: i, l
    max_in=un
    min_in=un
    DO l = 1, l_max
       DO i = 1, SIZE(un)
          maxn(i) = MAXVAL(max_in(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
          minn(i) = MINVAL(min_in(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
       END DO
       IF (l==l_max) EXIT
       max_in = maxn
       min_in = minn
    END DO
    
  END SUBROUTINE maxmin

  SUBROUTINE relax(stiff,un,maxn,minn, opt_EV_RES)
    USE mesh_handling
    USE matrix_type
    IMPLICIT NONE
    TYPE(matrice_bloc),         INTENT(IN)  :: stiff
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(:), OPTIONAL    :: opt_EV_RES
    REAL(KIND=8), DIMENSION(:)              :: minn
    REAL(KIND=8), DIMENSION(:)              :: maxn
    REAL(KIND=8), DIMENSION(SIZE(un))       :: alpha, denom
    INTEGER      :: i, j, p, ps, pe
    CHARACTER(*), PARAMETER :: limiter_type='avg'
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
    IF (PRESENT(opt_EV_RES)) THEN
       alpha = alpha*(1-opt_EV_RES)
    END IF
    SELECT CASE(limiter_type)
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
    CASE DEFAULT
       WRITE(*,*) ' BUG in relax'
       STOP
    END SELECT
    IF (PRESENT(opt_EV_RES)) THEN
       denom = denom*(1-opt_EV_RES)
    END IF
    !maxn = MIN(up_relax  *(maxn+ABS(maxn))/2 + down_relax*(maxn-ABS(maxn))/2,maxn + 2*ABS(denom)/2)
    !minn = MAX(down_relax*(minn+ABS(minn))/2 + up_relax  *(minn-ABS(minn))/2,minn - 2*ABS(denom)/2)
    maxn = maxn + 2*ABS(denom)
    minn = minn - 2*ABS(denom)
    !maxn = MIN(inputs%global_max,maxn)
    !minn = MAX(inputs%global_min,minn)
 
  END SUBROUTINE RELAX

  SUBROUTINE compute_relax_coeff(lumped)
    USE mesh_handling
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: lumped
    REAL(KIND=8), DIMENSION(SIZE(lumped))  :: ratio
    REAL(KIND=8) :: DD
    ALLOCATE(up_relax(mesh%np), down_relax(mesh%np))
    DD = SUM(lumped)
    IF (k_dim==1) THEN
       ratio = 2*(lumped/DD)*SQRT(lumped/DD)
       up_relax   = MIN(1.d0 + 2*ratio,up_relax_hard)
       down_relax = MAX(1.d0 - 2*ratio,down_relax_hard)
    ELSE IF (k_dim==2) THEN
       ratio = 2*SQRT(SQRT(lumped/DD))**3
       up_relax   = MIN(1.d0 + 4*ratio,up_relax_hard)
       down_relax = MAX(1.d0 - 4*ratio,down_relax_hard)
    END IF
  END SUBROUTINE compute_relax_coeff
  
END MODULE limiters

