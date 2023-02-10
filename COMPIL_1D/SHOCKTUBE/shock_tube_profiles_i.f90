MODULE  shock_tube_profiles_i

CONTAINS

SUBROUTINE  shock_tube_profiles (r_l, v_l, P_l,   &
        r_c_l, r_c_r, v_c, P_c,  r_r, v_r, P_r,   &
                                 x_D, t, xx,      &
                                 rr,  vv,  PP,  ee)

!  ...............
!                !
!                !
!                !.............
!
!  xa --------- x_D --------- xb

!  The left-going wave is always a rarefaction wave while the
!  right-going wave is always a shock wave

   USE  ideal_polytropic_gas

   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: r_l, v_l, P_l,  &
      r_c_l, r_c_r, v_c, P_c,  r_r, v_r, P_r,  &
                               x_D, t

   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: xx
   REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rr,  vv,  PP,  ee

   REAL(KIND=8), DIMENSION(SIZE(xx)) :: rr0, vv0, PP0, ee0

   REAL(KIND=8) :: c_l, c_r, c_c_l, v_f_l, v_f_r, v_sw, c, &
                   x, dx, y_f_l, y_f_r, y_cd, y_sw, y

   INTEGER :: i, n, m_points, N_time_levels = 1


   OPEN (UNIT = 30, FILE = 'shock_profile.m' )

   m_points = SIZE(xx)

   c_l = SQRT(gamma_ * P_l/r_l)
   c_r = SQRT(gamma_ * P_r/r_r)

   c_c_l = SQRT(gamma_ * P_c/r_c_l)

!  determine  v_f_l,  v_f_r,  v_sw

   v_f_l = v_l  -  c_l

   v_f_r = v_c  -  c_c_l

   v_sw = (r_r*v_r - r_c_r*v_c)/(r_r - r_c_r)

   !WRITE (*,*) 'v_sw = ' , v_sw

!  initial discontinuous profiles of constant states

   DO i = 1, m_points;   x = xx(i)

      IF (x <= x_D) THEN ! left constant state
         rr0(i) = r_l
         vv0(i) = v_l
         PP0(i) = P_l
         ee0(i) = P_l/gamma_m1_  +  r_l * v_l**2 /2
      ELSE ! right constant state
         rr0(i) = r_r
         vv0(i) = v_r
         PP0(i) = P_r
         ee0(i) = P_r/gamma_m1_  +  r_r * v_r**2 /2
      ENDIF

   ENDDO


!  Subsequent time level(s)

   DO n = 1, N_time_levels;   !!! t = n * dt

      y_f_l = v_f_l * t
      y_f_r = v_f_r * t
      y_cd  = v_c   * t
      y_sw  = v_sw  * t

      DO i = 1, m_points

         y = xx(i) - x_D

         IF (y <= y_f_l) THEN ! to the left of the rarefaction fan
            rr(i) = r_l
            vv(i) = v_l
            PP(i) = P_l
            ee(i) = P_l/gamma_m1_  +  r_l * v_l**2 /2
            CYCLE
         ENDIF

         IF (y_f_l < y  .AND.  y <= y_f_r) THEN ! rarefaction fan
            vv(i) = v_l  +  (v_c - v_l) * (y - y_f_l) / (y_f_r - y_f_l)
            c = c_l  -  alpha_ * (vv(i) - v_l)
            rr(i) = r_l * (c/c_l)**(1/alpha_)
            PP(i) = P_l * (c/c_l)**(1/tau_)
            ee(i) = PP(i)/gamma_m1_  +  rr(i) * vv(i)**2 /2
            CYCLE
         ENDIF

         IF (y_f_r < y  .AND.  y <= y_cd) THEN ! between rarefaction and contact
            rr(i) = r_c_l
            vv(i) = v_c
            PP(i) = P_c
            ee(i) = P_c/gamma_m1_  +  r_c_l * v_c**2 /2
            CYCLE
         ENDIF

         IF (y_cd < y  .AND.  y <= y_sw) THEN ! between contact and shock
            rr(i) = r_c_r
            vv(i) = v_c
            PP(i) = P_c
            ee(i) = P_c/gamma_m1_  +  r_c_r * v_c**2 /2
            CYCLE
         ENDIF

         IF (y > y_sw) THEN ! to the right of the shock
            rr(i) = r_r
            vv(i) = v_r
            PP(i) = P_r
            ee(i) = P_r/gamma_m1_  +  r_r * v_r**2 /2
         ENDIF

      ENDDO


!      WRITE (30,*) 'clear'

!      WRITE (30,*) 'xx=['
!      DO i = 1, m_points;  WRITE (30,*) xx(i);   ENDDO
!      WRITE (30,*) '];'

!      WRITE (30,*) 'rr0=['
!      DO i = 1, m_points;  WRITE (30,*) rr0(i);  ENDDO
!      WRITE (30,*) '];'

!      WRITE (30,*) 'vv0=['
!      DO i = 1, m_points;  WRITE (30,*) vv0(i);  ENDDO
!      WRITE (30,*) '];'

!      WRITE (30,*) 'PP0=['
!      DO i = 1, m_points;  WRITE (30,*) PP0(i);  ENDDO
!      WRITE (30,*) '];'

!      WRITE (30,*) 'ee0=['
!      DO i = 1, m_points;  WRITE (30,*) ee0(i);  ENDDO
!      WRITE (30,*) '];'

!      WRITE (30,*) 'rr=['
!      DO i = 1, m_points;  WRITE (30,*) rr(i);   ENDDO
!      WRITE (30,*) '];'

!      WRITE (30,*) 'vv=['
!      DO i = 1, m_points;  WRITE (30,*) vv(i);   ENDDO
!      WRITE (30,*) '];'

!      WRITE (30,*) 'PP=['
!      DO i = 1, m_points;  WRITE (30,*) PP(i);   ENDDO
!      WRITE (30,*) '];'

!      WRITE (30,*) 'ee=['
!      DO i = 1, m_points;  WRITE (30,*) ee(i);   ENDDO
!      WRITE (30,*) '];'


!      WRITE (30,*) 'figure(1) '
!      WRITE (30,*) ' assex=[ '
!      WRITE (30,*) xa
!      WRITE (30,*) xb
!      WRITE (30,*) '0.'
!      WRITE (30,*) 1.25*MAXVAL(rr)
!      WRITE (30,*)  ' ]; '
!      WRITE (30,*) ' plot(xx, rr0,''g:'', xx, rr,''r'') '
!      WRITE (30,*) ' axis(assex) '
!      WRITE (30,*) ' title(''Density'') '
!      WRITE (30,*) ' print Density.ps -f1'

!      WRITE (30,*) ' figure(2) '
!      WRITE (30,*) ' assex=[ '
!      WRITE (30,*) xa
!      WRITE (30,*) xb
!      WRITE (30,*) '0.'
!      WRITE (30,*) 1.25*MAXVAL(vv)
!      WRITE (30,*)  ' ]; '
!      WRITE (30,*) ' plot(xx, vv0,''g:'', xx, vv,''r'') '
!      WRITE (30,*) ' axis(assex) '
!      WRITE (30,*) ' title(''Velocity'') '
!      WRITE (30,*) ' print Velocity.ps -f2'

!      WRITE (30,*) ' figure(3) '
!      WRITE (30,*) ' assex=[ '
!      WRITE (30,*) xa
!      WRITE (30,*) xb
!      WRITE (30,*) '0.'
!      WRITE (30,*) 1.25*MAXVAL(PP)
!      WRITE (30,*)  ' ]; '
!      WRITE (30,*) ' plot(xx, PP0,''g:'', xx, PP,''r'') '
!      WRITE (30,*) ' axis(assex) '
!      WRITE (30,*) ' title(''Pressure'') '
!      WRITE (30,*) ' print Pressure.ps -f3'


!      WRITE (30,*) ' figure(4) '
!      WRITE (30,*) ' assex=[ '
!      WRITE (30,*) xa
!      WRITE (30,*) xb
!      WRITE (30,*) '0.'
!      WRITE (30,*) 1.25*MAXVAL(ee)
!      WRITE (30,*)  ' ]; '
!      WRITE (30,*) ' plot(xx, ee0,''g:'', xx, ee,''r'') '
!      WRITE (30,*) ' axis(assex) '
!      WRITE (30,*) ' title(''Total Energy Density'') '
!      WRITE (30,*) ' print Total_Energy_Density.ps -f4'


      CLOSE (30)


   ENDDO


END SUBROUTINE  shock_tube_profiles

!//////////////////////////////////////////////////////////////////////////////

SUBROUTINE  Riemann_tube_profiles (r_l, v_l, P_l,   &
          r_c_l, r_c_r, v_c, P_c,  r_r, v_r, P_r,   &
                           l_r_w,  x_D, t, xx,      &
                                   rr,  vv,  PP,  ee)

!  ...............
!                !
!                !
!                !.............
!
!  xa --------- x_D --------- xb

!  The left-going wave can be a rarefaction wave or a shock wave and
!  the right-going wave can be a rarefaction wave or a shock wave
!  The different cases are distingushed by the values in the
!  two-character parameter variable  l_r_w(1:2) = 'r'  or  's'

!  The program deals also with the degenerate cases when one or
!  two or all three waves are absent from the solution of the
!  Riemann problem.
!  In all these cases, the exact Riemann solver considers the lacking
!  first or third wave or both waves as a shock wave and returns
!
!  l_r_w(1) = 's'  or/and  l_r_w(2) = 's'
!
!  Therefore, the special treatment for dealing with an absent first
!  or third wave (and also with uniform flow with no wave at all)
!  is done under the case of a shock wave.

!  Design and coding by L. Quartapelle and L. Vigevano

!  Modifications for the degenerate cases by Marica Pelanti

   USE  ideal_polytropic_gas

   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: r_l, v_l, P_l,  &
      r_c_l, r_c_r, v_c, P_c,  r_r, v_r, P_r,  &
                               x_D, t

   CHARACTER,    DIMENSION(:), INTENT(IN)  :: l_r_w
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: xx
   REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rr,  vv,  PP,  ee

   REAL(KIND=8) :: c_l,  c_c_l,  c_c_r,  c_r,  &
                   v_fl_l,  v_fr_l,  v_sw_l,   v_fl_r,  v_fr_r,  v_sw_r,   &
                   y_fl_l,  y_fr_l,  y_sw_l,   y_fl_r,  y_fr_r,  y_sw_r,   &
                   y_l,  y_l_cd,  y_cd,  y_r_cd,  y_r,   c,  x,  dx,  y

   INTEGER :: i, m_points

   m_points = SIZE(xx)

   c_l   = SQRT(gamma_ * P_l/r_l)
   c_c_l = SQRT(gamma_ * P_c/r_c_l)
   c_c_r = SQRT(gamma_ * P_c/r_c_r)
   c_r   = SQRT(gamma_ * P_r/r_r)

   IF (l_r_w(1) == 'r')  THEN ! the left-going wave is a rarefaction

      v_fl_l  =  v_l  -  c_l;      y_fl_l  =  v_fl_l * t
      v_fr_l  =  v_c  -  c_c_l;    y_fr_l  =  v_fr_l * t
      y_l  =  y_fl_l;   y_l_cd  =  y_fr_l

   ELSE ! the left-going wave is a shock

      IF (ABS(r_l - r_c_l) <= 1.0E-9) THEN
         ! the left wave is absent
         v_sw_l = v_c
      ELSE
         v_sw_l  =  (r_l*v_l - r_c_l*v_c)/(r_l - r_c_l)
      ENDIF

      y_sw_l  =  v_sw_l * t
      y_l  =  y_sw_l;   y_l_cd  =  y_sw_l
      !WRITE(*,*) 'v_sw_l', v_sw_l

   ENDIF

   y_cd  =  v_c * t

   IF (l_r_w(2) == 'r')  THEN ! the right-going wave is a rarefaction

      v_fl_r  =  v_c  +  c_c_r;    y_fl_r  =  v_fl_r * t
      v_fr_r  =  v_r  +  c_r;      y_fr_r  =  v_fr_r * t
      y_r_cd  =  y_fl_r;   y_r  =  y_fr_r

   ELSE ! the right-going wave is a shock

      IF (ABS(r_r - r_c_r) <= 1.0E-9) THEN
         ! the right wave is absent
         v_sw_r = v_c
      ELSE
         v_sw_r  =  (r_r*v_r - r_c_r*v_c)/(r_r - r_c_r)
      ENDIF

      y_sw_r  =  v_sw_r * t
      y_r_cd  =  y_sw_r;   y_r  =  y_sw_r
      !WRITE(*,*) 'v_sw_r' , v_sw_r

   ENDIF


   DO i = 1, m_points

      y = xx(i) - x_D

      IF (y <= y_l) THEN ! to the left of the left-going wave
         rr(i) = r_l
         vv(i) = v_l
         PP(i) = P_l
         ee(i) = P_l/gamma_m1_  +  r_l * v_l**2 /2
         CYCLE
      ENDIF




      IF (l_r_w(1) == 'r'  .AND.  y_fl_l < y  .AND.  y <= y_fr_l) THEN ! rarefaction fan
         vv(i) = v_l  +  (v_c - v_l) * (y - y_fl_l) / (y_fr_l - y_fl_l)
         c = c_l  -  alpha_ * (vv(i) - v_l)
         rr(i) = r_l * (c/c_l)**(1/alpha_)
         PP(i) = P_l * (c/c_l)**(1/tau_)
         ee(i) = PP(i)/gamma_m1_  +  rr(i) * vv(i)**2 /2
         CYCLE
      ENDIF




      IF (y_l_cd < y  .AND.  y <= y_cd) THEN ! between left wave and contact
         rr(i) = r_c_l
         vv(i) = v_c
         PP(i) = P_c
         ee(i) = P_c/gamma_m1_  +  r_c_l * v_c**2 /2
         CYCLE
      ENDIF

      IF (y_cd < y  .AND.  y <= y_r_cd) THEN ! between contact and right wave
         rr(i) = r_c_r
         vv(i) = v_c
         PP(i) = P_c
         ee(i) = P_c/gamma_m1_  +  r_c_r * v_c**2 /2
         CYCLE
      ENDIF

      IF (l_r_w(2) == 'r'  .AND.  y_fl_r < y  .AND.  y <= y_fr_r) THEN ! rarefaction fan
         vv(i) = v_r  +  (v_c - v_r) * (y - y_fr_r) / (y_fl_r - y_fr_r)
         c = c_r  +  alpha_ * (vv(i) - v_r)
         rr(i) = r_r * (c/c_r)**(1/alpha_)
         PP(i) = P_r * (c/c_r)**(1/tau_)
         ee(i) = PP(i)/gamma_m1_  +  rr(i) * vv(i)**2 /2
         CYCLE
      ENDIF

      IF (y > y_r) THEN ! to the right of the right-going wave
         rr(i) = r_r
         vv(i) = v_r
         PP(i) = P_r
         ee(i) = P_r/gamma_m1_  +  r_r * v_r**2 /2
      ENDIF

   ENDDO


END SUBROUTINE  Riemann_tube_profiles

END MODULE  shock_tube_profiles_i
