  SUBROUTINE compute_greedy_dij(un)
    USE boundary_conditions !===> to access to the flux
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8):: cmin, rho_e_min 
    REAL(KIND=8), DIMENSION(k_dim+2,mesh%np), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(k_dim+2,k_dim,mesh%np)       :: vv
    REAL(KIND=8), DIMENSION(k_dim)                       :: nij
    REAL(KIND=8), DIMENSION(k_dim+2)                     :: Plr, usafe
    REAL(KIND=8)                                         :: p_min, rho_min, p_max, rho_max
    REAL(KIND=8)                                         :: romin, romax
    REAL(KIND=8) :: mr, ml, rhor, rhol, er, el, pr, pl, ul, ur, rhosmall, norm_cij, E_small, ttt
    REAL(KIND=8) :: lambda_max, lambda_max_safe, lambda_max_small, lbd1, lbd2
    REAL(KIND=8) :: Sl, Sr, S_small, aaa_entrop, bbb_entrop
    REAL(KIND=8) :: small = 1.d-7
    INTEGER :: i, j, p, k, comp

    DO comp = 1, k_dim+2
       vv(comp,:,:)=flux(comp,un)
    END DO

    DO i = 1, mesh%np
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          IF (i.EQ.j) THEN
             dij%aa(p) = 0.d0
             CYCLE
          END IF
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
          pr = ABS(er-0.5d0*rhor*ur**2)*(gamma-1.d0)
          pl = ABS(el-0.5d0*rhol*ul**2)*(gamma-1.d0)
          !===Compute p_min,rho_min,p_max,rho_max,lambda_max_safe
          CALL ptilde(gamma,rhol,ul,pl,rhor,ur,pr,p_min,rho_min,p_max,rho_max,lambda_max_safe)
          rhosmall = rho_max*small
          lambda_max_small = lambda_max_safe*small
          !===Compute bar states
          DO comp = 1, k_dim+2
             Plr(comp)= - SUM((vv(comp,:,j)-vv(comp,:,i))*nij)/2.d0 
          END DO
          usafe = (un(:,i)+un(:,j))/2.d0 + Plr/lambda_max_safe
          romin = min(rho_min,usafe(1))
          romax = max(rho_max,usafe(1))
          !===First lambda based on density
          IF (ABS(rhol+rhor-2.d0*romin).LE.rhosmall) THEN
             lbd1 = lambda_max_safe
          ELSE
             lbd1 = lambda_max_safe*(rhol+rhor-2.d0*usafe(1))/(rhol+rhor-2.d0*romin)
          END IF
          IF (ABS(rhol+rhor-2.d0*romax).LE.rhosmall) THEN
             lbd2 = lambda_max_safe
          ELSE
             lbd2 = lambda_max_safe*(2.d0*usafe(1)-rhol-rhor)/(2.d0*romax-rhol-rhor)
          END IF
          lambda_max = max(lbd1, lbd2, lambda_max_small)
          IF (lambda_max>lambda_max_safe) THEN
             WRITE(*,*) ' THERE IS A BUG 1', lambda_max, lambda_max_safe
          END IF
          !===Second lambda based on internal energy rho e
!!$          rho_e_min = (1-small)*p_min/(gamma-1)
!!$          ttt = 1.d0/lambda_max
!!$          CALL rho_e_limit(un(:,i),un(:,j),Plr,ttt,rho_e_min)
!!$          lambda_max = MAX(lambda_max,1/ttt)
!!$          IF (lambda_max>lambda_max_safe) THEN
!!$             WRITE(*,*) ' THERE IS A BUG 2'
!!$          END IF
          !===Third lambda based on minimum principle on entropy
          !===Compute cmin
          cmin = (1.d0/(gamma-1.d0))*min(pl/rhol**(gamma),pr/rhor**(gamma)) !===gamma is defined in lambda_module
          cmin = cmin*(1-small)
          E_small = small*(un(k_dim+2,i)+un(k_dim+2,j))
          IF (lambda_max.GE.lambda_max_small) THEN
             ttt = 1.d0/lambda_max
          ELSE
             write(*,*) ' lambda_max seems to be equal to zero', lambda_max
             ttt = 1.d0/lambda_max_small
             stop
          END IF
          CALL Newton_secant(un(:,i),un(:,j),Plr,ttt,psi_entrop,psi_entrop_prime,E_small)
          IF (ttt.le.0.d0) THEN
             WRITE(*,*) ' BUG ttt<0', ttt
             STOP
          END IF
          lambda_max = MAX(lambda_max,1.d0/ttt)
          !===Fourth lambda based on entropy inequality
          Sl = entrop(un(:,i),0d0,0.d0,0.d0)
          Sr = entrop(un(:,j),0d0,0.d0,0.d0)
          aaa_entrop = -(Sl+Sr)/2.d0
          bbb_entrop = (ur*Sr-ul*Sl)/2.d0
          IF (aaa_entrop>0.d0) THEN
             aaa_entrop = aaa_entrop*(1-small)
          ELSE
             aaa_entrop = aaa_entrop*(1+small)
          END IF
          S_small = small*(ABS(Sl)+ABS(Sr))
          ttt = 1.d0/lambda_max
          CALL entrop_ineq(un(:,i),un(:,j),Plr,ttt,aaa_entrop,bbb_entrop,S_small)
          lambda_max = MAX(lambda_max,1.d0/ttt)
          !===Definition of dij
          dij%aa(p) = norm_cij*lambda_max
       END DO
    END DO
    CALL transpose_op(dij,'max')
    DO i = 1, mesh%np
       dij%aa(diag(i)) = -SUM(dij%aa(dij%ia(i):dij%ia(i+1)-1))
    END DO

  CONTAINS
    FUNCTION psi_entrop(u) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u
      REAL(KIND=8)                                 :: psi
      psi = u(k_dim+2) - SUM(u(2:k_dim+1)**2)/(2.d0*u(1)) - cmin*u(1)**gamma
    END FUNCTION psi_entrop

    FUNCTION psi_entrop_prime(Pij,u) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u, Pij
      REAL(KIND=8)                                 :: psi
      psi = Pij(k_dim+2) - SUM(u(2:k_dim+1)*Pij(2:k_dim+1))/u(1) &
           + Pij(1)*SUM(u(2:k_dim+1)**2)/(2*u(1)**2) &
           - cmin*gamma*Pij(1)*u(1)**(gamma-1.d0)
    END FUNCTION psi_entrop_prime
  END SUBROUTINE compute_greedy_dij

  SUBROUTINE Newton_secant(ui, uj, Pij,limiter,psi_func,psi_prime_func,psi_small)
    IMPLICIT NONE
    INTERFACE
       FUNCTION psi_func(u) RESULT(psi)
         USE space_dim
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN)  :: u
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_func
       FUNCTION psi_prime_func(Pij,u) RESULT(psi)
         USE space_dim
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN)  :: u, Pij
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_prime_func
    END INTERFACE
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: Pij
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: ui, uj
    REAL(KIND=8), INTENT(IN)    :: psi_small
    REAL(KIND=8), INTENT(INOUT) :: limiter
    REAL(KIND=8), DIMENSION(k_dim+2) ::  ul, ur, ulr
    REAL(KIND=8) :: psil, psir, ll, lr, llold, lrold
    
    ul = ui
    ur = uj
    ulr = (ul+ur)/2
    
    lr = limiter
    ur = ulr + lr*Pij
    psir = psi_func(ur)
    !IF (ur(1).le.0.d0) THEN
    !   write(*,*) ' density negative in Newton secant', ur(1)
    !   stop
    !END IF
    
    IF (psir.GE.-psi_small) THEN
       !===input limiter is okay
       !IF (ur(1).le.0.d0 .or. ur(k_dim+2).le.0.d0) THEN
       !   write(*,*) ' density negative in Newton secant', ur(1), ur(k_dim+2)
       !   stop
       !end if
       RETURN
    END IF
    ll = 0.d0
    ul = ulr
    psil = max(psi_func(ulr),0.d0)
 
    DO WHILE (ABS(psil-psir) .GT. psi_small)
       llold = ll
       lrold = lr
       ll = ll - psil*(lr-ll)/(psir-psil)
       lr = lr - psir/psi_prime_func(Pij,ur)
       IF (ll.GE.lr) THEN
          ll = lr
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
       ul = ulr + ll*Pij
       ur = ulr + lr*Pij
       psil = psi_func(ul)
       psir = psi_func(ur)
    END DO
    IF (psir.GE.-psi_small) THEN
       limiter = lr
    ELSE
       limiter = ll
    END IF
  END SUBROUTINE Newton_secant

  SUBROUTINE entrop_ineq(ui,uj,Pij,limiter,aaa,bbb,psi_small)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: Pij
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: ui, uj
    REAL(KIND=8), INTENT(IN)    :: aaa, bbb, psi_small
    REAL(KIND=8), INTENT(INOUT) :: limiter
    REAL(KIND=8), DIMENSION(k_dim+2) ::  ul, ur, ulr
    REAL(KIND=8) :: psil, psir, dpsir, ll, lr, llold, lrold
    integer      :: k
    ul = ui
    ur = uj
    ulr = (ul+ur)/2

    lr = limiter
    ur = ulr + lr*Pij
    psir = entrop(ur,aaa,bbb,lr)    
    IF (psir.LE. psi_small) THEN
       !===input limiter is okay
       RETURN
    END IF
    ll = 0.d0
    ul = ulr
    psil = max(entrop(ulr,aaa,bbb,ll),0.d0)
 
    DO WHILE (ABS(psil-psir) .GT. psi_small)
       llold = ll
       lrold = lr
       ll = ll - psil*(lr-ll)/(psir-psil)
       lr = lr - psir/entrop_prime(Pij,ur,bbb)
       IF (ll.GE.lr) THEN
          ll = lr
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
       ul = ulr + ll*Pij
       ur = ulr + lr*Pij
       psil = entrop(ul,aaa,bbb,ll)
       psir = entrop(ur,aaa,bbb,lr)
    END DO
    IF (psir.LE.psi_small) THEN
       limiter = lr
       ul = ur
    ELSE
       limiter = ll
    END IF
    
    if (entrop(ulr+limiter*pij,aaa,bbb,limiter)>2*psi_small) THEN
       WRITE(*,*) ' BUG  entrop ineq. violated',  entrop(ulr+limiter*pij,aaa,bbb,limiter), ABS(psil-psir)
       stop
    END if
  END SUBROUTINE entrop_ineq
  

  SUBROUTINE ptilde(gamma,rhol,ul,pl,rhor,ur,pr,p_min,rho_min,p_max,rho_max,lambda_max)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: rhol, ul, pl, rhor, ur, pr
    REAL(KIND=8), INTENT(OUT):: p_min,rho_min,p_max,rho_max,lambda_max
    REAL(KIND=8)             :: pmin, pmax, rhomin, rhomax, amin, amax, phimin, phimax
    REAL(KIND=8)             :: al, capAl, capBl, ar, capAr, capBr, capAmin, capBmin, p1, p2
    REAL(KIND=8)             :: ratio, ptil, lambdal,  lambdar 
    REAL(KIND=8)             :: gamma, exponent, one_m_exponent, gm1_over_gp1
    REAL                     :: Mgas, Mgasptwo
 
    Mgas = 2/(gamma-1)
    Mgasptwo = Mgas+2
    exponent=(gamma-1)/(2*gamma)
    one_m_exponent = 1-exponent
    gm1_over_gp1 =(gamma-1)/(gamma+1)
    
    al = SQRT(gamma*pl/rhol)
    capAl = 2/((gamma+1)*rhol)
    capBl = pl*(gamma-1)/(gamma+1)

    ar = SQRT(gamma*pr/rhor)
    capAr = 2/((gamma+1)*rhor)
    capBr = pr*(gamma-1)/(gamma+1)

    IF (pl.LE.pr) THEN
       pmin   = pl
       rhomin = rhol
       pmax   = pr
       rhomax = rhor
    ELSE
       pmin   = pr
       rhomin = rhor
       pmax   = pl
       rhomax = rhol
    END IF
    capAmin = 2/((gamma+1)*rhomin)
    capBmin = pmin*(gamma-1)/(gamma+1)
    amin = SQRT(gamma*pmin/rhomin)
    amax = SQRT(gamma*pmax/rhomax)
    ratio = (pmin/pmax)**exponent
    phimin = (2/(gamma-1.d0))*amax*(ratio-1.d0) + ur-ul
    !===Check phi(pmin)>0. Phi(p) defined by (4.1) JLG/BP JCP 321 (2016) 908-926
    IF (phimin.GE.0) THEN !===Two expansions
       p2 = pmin*((amin + amax - (ur-ul)*(gamma-1)/2)/(amin + amax*ratio))**(Mgas+2) !=== =p* (exact)
       rho_min = min(rhol*(p2/pl)*(1/gamma),rhor*(p2/pr)*(1/gamma))
       !TEST
       !rho_min = min((rhol*(p2/pl)*(1/gamma) + rhol)/2,(rhor*(p2/pr)*(1/gamma)+rhor)/2)
       !TEST
       rho_max = max(rhol,rhor)
       p_min = p2
       p_max = pmax
       lambdal = ul - al
       lambdar = ur + ar
       lambda_max = MAX(ABS(lambdal),ABS(lambdar))
       RETURN
    END IF
    !===We continue with phimin>0
    phimax = (pmax-pmin)*SQRT(capAmin/(pmax+capBmin)) + ur-ul !===(3.3) with (3.4) for p=pmax
    ptil = pmin*((amin + amax - (ur-ul)*(gamma-1)/2)/(amin + amax*ratio))**(Mgas+2)
    IF (phimax < 0.d0) THEN
       p1 = pmax
       p2 = ptil
    ELSE
       p1 = pmin
       p2 = MIN(pmax,ptil)
    END IF
    p1 = MAX(p1,p2-phi(p2,ul,pl,ur,pr)/phi_prime(p2,pl,pr))
    !===Assign min and max
    rho_min = min(rhol,rhor)
    rho_max = max(rhol,rhor)
    !===First wave
    !===Phi(pmin)<0: shock
    rho_max = max(rho_max,rhomin*(gm1_over_gp1 + p1/pmin)/(gm1_over_gp1*(p1/pmin) +1.d0))  !===(4.19) Toro, p.122
    !===Second wave
    IF (phimax < 0.d0) THEN
       !===Phi(pmax)<0: shock
       rho_max = max(rho_max,rhomax*(gm1_over_gp1 + p1/pmax)/(gm1_over_gp1*(p1/pmax) +1.d0))
       p_min = pmin
       p_max = p1
    ELSE
       !===Phi(pmax)>0: expansion
       rho_min = min(rho_min,rhomax*(p2/pmax)*(1/gamma))
       p_min = pmin
       p_max = pmax
    END IF

    lambdal = ul - al*SQRT(1.d0+MAX((p2-pl)/pl,0.d0)*one_m_exponent)
    lambdar = ur + ar*SQRT(1.d0+MAX((p2-pr)/pr,0.d0)*one_m_exponent)
    lambda_max = MAX(ABS(lambdal),ABS(lambdar))
  CONTAINS

    FUNCTION phi(p,ul,pl,ur,pr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: p, ul, pl, ur, pr
      REAL(KIND=8)             :: vv, fl, fr
      IF (p>pl) THEN
         fl = (p-pl)*SQRT(capAl/(p+capBl))
      ELSE
         fl = (2*al/(gamma-1))*((p/pl)**exponent-1)
      END IF
      IF (p>pr) THEN
         fr = (p-pr)*SQRT(capAr/(p+capBr))
      ELSE
         fr = (2*ar/(gamma-1))*((p/pr)**exponent-1)
      END IF
      vv = fl + fr + ur - ul
    END FUNCTION phi

    FUNCTION phi_prime(p,pl,pr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: p, pl, pr
      REAL(KIND=8)             :: vv, fl, fr
      IF (p>pl) THEN
         fl = SQRT(capAl/(p+capBl))*(1-(p-pl)/(2*(capBl+p)))
      ELSE
         fl = (al/(gamma*pl))*(p/pl)**(-(gamma+1)/(2*gamma))
      END IF
      IF (p>pr) THEN
         fr = SQRT(capAr/(p+capBr))*(1-(p-pr)/(2*(capBr+p)))
      ELSE
         fr = (ar/(gamma*pr))*(p/pr)**(-(gamma+1)/(2*gamma))
      END IF
      vv = fl + fr
    END FUNCTION phi_prime
  END SUBROUTINE ptilde
