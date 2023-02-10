PROGRAM vp_Diagram

 USE vdw_gas

 !-------------------------------------------------------------------------
 IMPLICIT NONE 

 REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v, v_gas, v_liq, &
                                            p, p_sat, p_curve
 
 REAL(KIND=8) :: delta, t, sigma, dumm

 REAL(KIND=8), PARAMETER :: d_delta = 0.003, d_t = 0.01, d_sigma = 0.0007

 INTEGER :: i, j, k, nj
 INTEGER, PARAMETER :: idf = 11, n = 1000
 !-------------------------------------------------------------------------

  ! Loading Saturation Curve data from file
  ALLOCATE (v_gas(n), v_liq(n), p_sat(n), p_curve(n))

  OPEN (UNIT=idf, FILE='curva_di_saturazione.dat', STATUS='unknown')

    READ(idf,*)
    
    DO i = 1, n
      READ(idf,*) dumm, p_sat(i), v_gas(i), v_liq(i)       
    ENDDO

  CLOSE(idf)



  ALLOCATE (p(n), v(n))

  OPEN (UNIT=idf, FILE='vp_diagram.plt', STATUS='unknown')
    
    WRITE(idf,*) 'TITLE = "G=0 Locus"'		         
    WRITE(idf,*) 'VARIABLES = "v", "p"' 	         



    ! Writing Saturation Curve

    WRITE(idf,*) 'ZONE T ="saturazione", I =', 2*n-1 

    DO i = 1, n-1
      WRITE(idf,*) v_liq(n-i), p_sat(n-i)
    ENDDO 

    DO i = 1, n
      WRITE(idf,*) v_gas(i), p_sat(i)
    ENDDO 



    DO j = 1, 1  !21

       delta = 0.03d0   !0.d0 + (j-1)*d_delta

       ! Writing Locus of G = 0. Polytropic Van der Waals
       ! gas 'c_v' is assumed constant.    

       DO i = 1, n

         v(i) = 1.d0 + (i-1)*0.002

         p(i) = Gzero_Locus(v(i), delta)

         p_curve(i) = Linear_Interpolation(p_sat, v_gas, v(i))

       ENDDO

       nj = COUNT(p > p_curve)

       WRITE(idf,*) ''                                        
       WRITE(idf,*) 'ZONE T="d=', delta, '", I=', nj      

       DO i = 1, n
         IF (p(i) > p_curve(i))  WRITE(idf,*) v(i), p(i)
       ENDDO



       ! Writing Van der Waals Isentropes for 20 different values
       ! of the parameter sigma.

       DO k = 1, 51
    
         sigma = 0.300 !+ (k-1)*d_sigma

         DO i = 1, n

           v(i) = 0.35 + (i-1)*0.00465

           p(i) = p_vdW__vs(v(i), sigma, delta)

         ENDDO

         WRITE(idf,*) ''                                        
         WRITE(idf,*) 'ZONE T="delta=', delta, ' sigma=', sigma, '", I=', n
    
         DO i = 1, n
           WRITE(idf,*) v(i), p(i)
         ENDDO

       ENDDO
    
    ENDDO


    
    ! Writing Van der Waals Isotherms for 20 different values
    ! of reduced temprature t,  starting from critical  value
    ! t = 1.

    DO j = 1, 21

      t = 1   !0.9962458333333334d0 + (j-1)*d_t

      DO i = 1, n

        v(i) = 0.35 + (i-1)*0.00465

    	p(i) = p_vdW__vt(v(i), t)

      ENDDO

      WRITE(idf,*) ''				             
      WRITE(idf,*) 'ZONE T="t=', t, '", I=', n

      DO i = 1, n
        WRITE(idf,*) v(i), p(i)
      ENDDO

    ENDDO


  CLOSE(idf)



  CONTAINS


  FUNCTION Linear_Interpolation(gg, xx, x) RESULT(gx)
 
  ! Compute the value of function gg at point x. gg  is  a vector  of  the
  ! values of the function evaluated at the points specified  in  xx.  The
  ! value g(x) is computed by mean of linear interpolation.
  !
  ! This algorith assumes that the grid values xx(i) are ardered in incre-
  ! asing manner
 
  !------------------------------------------------------------------------
  IMPLICIT NONE
 
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: gg, xx
  REAL(KIND=8),               INTENT(IN) :: x
 
  LOGICAL :: check
  INTEGER :: i
 
  REAL(KIND=8) :: x_l, x_r,  &
                  g_l, g_r, gx
  !------------------------------------------------------------------------
 
  check = .FALSE.   ! Manage the possibility that xp is not included
                    ! in the range xx(1) - xx( SIZE(xx) )

  DO i = 1, SIZE(xx) - 1

    IF ( (x >= xx(i)  .AND.  x <  xx(i+1)) .OR. &
         (x <  xx(i)  .AND.  x >= xx(i+1)) ) THEN

        x_l = xx(i);   x_r = xx(i+1)
        g_l = gg(i);   g_r = gg(i+1)

        check = .TRUE.

        EXIT

    ENDIF

  ENDDO

  IF (.NOT. check)  THEN

    PRINT*, ' The point x is external to the interval of the grid'
    PRINT*, ' Linear_Interpolation function is STOPPED.'
    PRINT*, ''
    STOP

  ENDIF

  ! Linear interpolation for computing gg(x)
  gx =  g_l + (x - x_l) * (g_r - g_l) / (x_r - x_l)


  END FUNCTION Linear_Interpolation


END PROGRAM vp_Diagram
