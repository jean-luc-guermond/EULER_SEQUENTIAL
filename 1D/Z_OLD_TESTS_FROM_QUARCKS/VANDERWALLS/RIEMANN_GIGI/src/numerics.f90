
! =========================================================================
!   Description: Numerical Procedures.
!
!   	 Author: Marco Fossati
!   		 Department of Aerospace Engineering
!   		 Politecnico di Milano
!   		 Via La Masa 34, 20156 Milano, ITALY
!   		 e-mail: fossati@aero.polimi.it
!
!          Year: 2005 November
! =========================================================================

MODULE numerics

   PRIVATE :: RK_qs,  RK_ck

   PRIVATE :: Newton_Sys2_Diag_Jac ! OLD IMPROPER IMPLEMENTATION

CONTAINS


FUNCTION  Cubic_Equation(A1, B1, C1,  single_real_root)  RESULT(real_roots)
 
   ! Soluzione dell'equazione algebrica di terzo grado
   !								  			       
   !	 X^3  +  A1 X^2  +  B1 X  +  C1  =  0
   !								  			       
   ! con A, B e C coefficienti REALI
   !
   ! Il programma calcola le radici  REALI (1 o 3)  dell'equazione
   ! di terzo grado e segnala  nella  variabile  logica  opzionale
   ! single_real_root quando esiste una sola radice reale.
   !
   ! Quando c'e' una sola radice reale, essa e' messa nella  prima
   ! componente real_roots(1), mentre le altre due componenti sono
   ! poste artificialmente uguali a zero
   !
   ! Trasformazione  X  -->  x = X + A1/3
   !
   !
   !	 x^3  +  b x  +  c  =  0
   !
   !
   !	 b  =  - A1^2/3  +  B1
   !
   !	 c  =  (2/27) A1^3  -  A1*B1 / 3  +  C1
   !
   !
   !	 x^3  =  3p x  +  2q
   !
   !
   !	 p  =  - b/3  =   A1^2/9   -  B1/3
   !
   !	 q  =  - c/2  = - A1^3/27  +  A1*B1 / 6  -  C1/2
 
   !-----------------------------------------------------------------------
   IMPLICIT NONE
 
   REAL(KIND=8), 	    INTENT(IN)  :: A1, B1, C1
   LOGICAL,	OPTIONAL, INTENT(OUT) :: single_real_root
 
   REAL(KIND=8), DIMENSION(3) :: real_roots
 
   REAL(KIND=8) :: p, q, S_q_p, r, s, t, phi,  &
                   one = 1,  pi = 3.14159265358979323846d0
   !-----------------------------------------------------------------------

   p = A1**2/9  -  B1/3
 
   q = - A1**3/27  +  A1*B1/6  -  C1/2
 
   s = q**2 - p**3
 
   IF (s >= 0) THEN  !  only one real solution
 
      IF (PRESENT (single_real_root)) single_real_root = .TRUE.
  
      S_q_p = SQRT(s)
 
      r = SIGN(one, q + S_q_p) * ABS(q + S_q_p)**(one/3)
      t = SIGN(one, q - S_q_p) * ABS(q - S_q_p)**(one/3)
 
      real_roots(1) = r + t
 
      !  z2 = - (r + t)/2  +  (r - t) * CMPLX(0, SQRT(3.d0))/2
      !  z3 = CONJG(z2)
 
      real_roots(2:3) = 0
 
   ELSE  !  three real solutions:  Francois Viete
 
      IF (PRESENT (single_real_root)) single_real_root = .FALSE.
 
      phi = ACOS(q/(p*SQRT(p)))
 
      real_roots(1) = COS(phi/3)
      real_roots(2) = COS((phi + 2*pi)/3)
      real_roots(3) = COS((phi + 4*pi)/3)
 
      real_roots = 2 * SQRT(p) * real_roots
  
   ENDIF
 
   real_roots = real_roots - A1/3
  
END FUNCTION  cubic_equation		







FUNCTION  Newton(ff, Df, x0, rel_tol, NM_Iter,  f_shift)  RESULT (root)
     
  ! Se il parametro OPZIONALE f_shift e' assente, il programma
  ! determina uno zero reale di una funzione reale f(x), ovverosia 
  ! una radice (root) dell'equazione 
  !
  !      f(root) = 0
  !
  ! mediante il metodo di Newton, scritto di forma incrementale.
  !
  ! Se invece il parametro OPZIONALE f_shift e' presente, si
  ! calcola la soluzione root dell'equazione ``shiftata''     
  !
  !      f(root) = f_shift 
  !
  ! secondo un suggerimento di Alberto Guardone. 

  ! ff(x) e' la funzione, Df(x) la sua derivata e 
  ! x0 un valore iniziale che deve essere sufficientemente
  ! vicino a root affinche' il metodo possa convergere.
  ! Quanto vicino, in generale nessuno lo sa a priori e
  ! quindi occorre analizzare ogni problema specifico,
  ! ossia la forma determinata della funzione f(x) in esame.  

  ! La radice e' calcolata con una precisione relativa specificata
  ! tramite il parametro rel_tol. 
  !
  ! Il numero di iterazioni da eseguire e' fornito nel parametro
  ! NM_Iter : in INPUT   NM_Iter e' il numero massimo di iterazioni
  !           in OUTPUT  NM_Iter e' il numero di iterazioni effettuate
  !

  !-----------------------------------------------------------------------
   IMPLICIT NONE
   
   INTERFACE 

      FUNCTION ff(x) RESULT(f)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: f
      END FUNCTION ff

      FUNCTION Df(x) RESULT(g)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: g
      END FUNCTION Df

   END INTERFACE


   REAL(KIND=8),           INTENT(IN)    :: x0, rel_tol
   INTEGER,                INTENT(INOUT) :: NM_Iter
   REAL(KIND=8), OPTIONAL, INTENT(IN)    :: f_shift 

   REAL(KIND=8) :: root
   
   REAL(KIND = 8) :: d_x, x, xm, int_rel_tol, f_sh = 0, derivata 
   INTEGER :: i, iter
   !-----------------------------------------------------------------------
  	
   ! Verifica della precisione relativa richiesta
   IF (rel_tol > EPSILON(x)) THEN
      int_rel_tol = rel_tol
   ELSE
      int_rel_tol = EPSILON(xm)*10
   ENDIF

   iter = NM_Iter
    
   IF (PRESENT(f_shift)) f_sh = f_shift  


   ! Metodo di Newton       
   x = x0

   DO i = 1, iter

      NM_Iter = i

      derivata = Df(x)

      IF (ABS(derivata) < EPSILON(x)) THEN

         IF (ABS(ff(x)) < EPSILON(x)) RETURN

         PRINT*, ''
         PRINT*, 'Error. NEWTON:'
         PRINT*, 'Derivative equal zero'
         PRINT*, ''
         OPEN(UNIT=11, FILE= 'stopped')
         WRITE(11,*) 'Newton method has computed a zero derivative'
         CLOSE(11)          
         STOP
    
      ENDIF


      d_x = - (ff(x) - f_sh)/derivata
    
  
      ! Nel caso di sistema la divisione dovra' essere sostituita
      ! dalla soluzione di un sistema lineare che effettueremo
      ! nel modo seguente : 
      ! CALL fattorizzazione ( Df ) ! Matrice Jacobiana
      ! CALL sostituzioni ( Df, -(ff - f_sh) ----> d_x )
      ! Ovviamente ABS( ) sara' sostituito da una norma
      ! opportuna del vettore dell'incognite

      x = x + d_x

      IF (ABS(d_x) < ABS(x) * int_rel_tol) THEN
	       
         root = x
         RETURN   
         
      ENDIF

   ENDDO
    
END FUNCTION  Newton 





FUNCTION  Custom_Newton(ff, Df, x0, rel_tol, NM_Iter,  f_shift)  RESULT (root)
     
  ! Se il parametro OPZIONALE f_shift e' assente, il programma
  ! determina uno zero reale di una funzione reale f(x), ovverosia 
  ! una radice (root) dell'equazione 
  !
  !      f(root) = 0
  !
  ! mediante il metodo di Newton, scritto di forma incrementale.
  !
  ! Se invece il parametro OPZIONALE f_shift e' presente, si
  ! calcola la soluzione root dell'equazione ``shiftata''     
  !
  !      f(root) = f_shift 
  !
  ! secondo un suggerimento di Alberto Guardone. 

  ! ff(x) e' la funzione, Df(x) la sua derivata e 
  ! x0 un valore iniziale che deve essere sufficientemente
  ! vicino a root affinche' il metodo possa convergere.
  ! Quanto vicino, in generale nessuno lo sa a priori e
  ! quindi occorre analizzare ogni problema specifico,
  ! ossia la forma determinata della funzione f(x) in esame.  

  ! La radice e' calcolata con una precisione relativa specificata
  ! tramite il parametro rel_tol. 
  !
  ! Il numero di iterazioni da eseguire e' fornito nel parametro
  ! NM_Iter : in INPUT    NM_Iter e' il numero massimo di iterazioni
  !           in OUTPUT   NM_Iter e' il numero di iterazioni effettuate
    
   USE problem_params
  
   !-----------------------------------------------------------------------
   IMPLICIT NONE
   
   INTERFACE 

      FUNCTION ff(x) RESULT(f)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: f
      END FUNCTION ff

      FUNCTION Df(x) RESULT(g)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: g
      END FUNCTION Df

   END INTERFACE


   REAL(KIND=8),           INTENT(IN)    :: x0, rel_tol
   INTEGER,      OPTIONAL, INTENT(INOUT) :: NM_Iter
   REAL(KIND=8), OPTIONAL, INTENT(IN)    :: f_shift

   REAL(KIND=8) :: root
   
   REAL(KIND = 8) :: d_x, x, xm, int_rel_tol, &
                     f_sh = 0, derivata

   INTEGER :: i, iter

   LOGICAL :: forced = .FALSE.
   !-----------------------------------------------------------------------

   WRITE(12,*) '   CUSTOM_NEWTON'

   ! Verifica della precisione relativa richiesta
   IF (rel_tol > EPSILON(x)) THEN
      int_rel_tol = rel_tol
   ELSE
      int_rel_tol = EPSILON(xm)*10
   ENDIF
    
   iter = NM_Iter
    
   IF (PRESENT(f_shift)) f_sh = f_shift  


   ! Newton Method
   x = x0
   WRITE(12,*) '     guess = ', x0

   DO i = 1, iter

      NM_Iter = i     

      !---Corrections to Newton Method due to LQ and MF---------------------------
      IF (x > 1.d+5) THEN

         WRITE(12,*) '     Root has exceeded 100000.'
         WRITE(12,*) '     continues...'
         RETURN

      ENDIF 

      ! Redefining  x since it has become not physical (less than 1).
      ! Correction suggested by Luigi Quartapelle.
      IF (x <= 1) THEN
       
         forced = .TRUE.
         x = 1.03

      ENDIF
      !---------------------------------------------------------------------------


      derivata = Df(x)

      IF (ABS(derivata) < EPSILON(x)) THEN

         IF (ABS(ff(x)) < EPSILON(x)) RETURN

         PRINT*, ''
         PRINT*, 'Error. CUSTOM_NEWTON:'
         PRINT*, 'Derivative equals zero.'
         PRINT*, ''
         OPEN(UNIT=11, FILE='stopped')
         WRITE(11,*) 'Newton Forced method has computed a zero derivative'
         CLOSE(11)          
         STOP
    
      ENDIF


      d_x = - (ff(x) - f_sh)/derivata


      !-Fixed Relaxation factor for Newton method. Marco F.----
      d_x = relax * d_x
      !--------------------------------------------------------


      x = x + d_x


      ! Further relaxation due to the fact that roots may  assume
      ! nonphysical value (<= 1). The overelaxation is made until
      ! roots becomes physical by halving
      ! 
      !IF (x <= 1)  THEN
      !
      !   x = x - d_x
      !   forced = .TRUE.
      !
      !   DO WHILE (x <= 1)
      !
      !      x = x + 0.5*d_x
      !
      !   ENDDO
      !
      !ENDIF 



      IF (ABS(d_x) < ABS(x) * int_rel_tol) THEN
       
         root = x

         IF (forced) WRITE(12,*) '     Newton Method required FORCING.'

         RETURN   
         
      ENDIF

   ENDDO

    
END FUNCTION  Custom_Newton







FUNCTION  Newton_Sys2_Diag_Jac(ff, Df, DDf, x0, rel_tol, NM_Iter)  RESULT (x)
     
  ! Se il parametro OPZIONALE f_shift e' assente, il programma
  ! determina uno zero reale di una funzione reale f(x), ovverosia 
  ! una radice (root) dell'equazione 
  !
  !      f(root) = 0
  !
  ! mediante il metodo di Newton, scritto di forma incrementale.
  !
  ! Se invece il parametro OPZIONALE f_shift e' presente, si
  ! calcola la soluzione root dell'equazione ``shiftata''     
  !
  !      f(root) = f_shift 
  !
  ! secondo un suggerimento di Alberto Guardone. 

  ! ff(x) e' la funzione, Df(x) la sua derivata e 
  ! x0 un valore iniziale che deve essere sufficientemente
  ! vicino a root affinche' il metodo possa convergere.
  ! Quanto vicino, in generale nessuno lo sa a priori e
  ! quindi occorre analizzare ogni problema specifico,
  ! ossia la forma determinata della funzione f(x) in esame.  

  ! La radice e' calcolata con una precisione relativa specificata
  ! tramite il parametro rel_tol. 
  !
  ! Il numero di iterazioni da eseguire e' fornito nel parametro
  ! NM_Iter : in INPUT   NM_Iter e' il numero massimo di iterazioni
  !           in OUTPUT  NM_Iter e' il numero di iterazioni effettuate
  !

   !-----------------------------------------------------------------------
   IMPLICIT NONE
    
   INTERFACE 

      FUNCTION ff(x) RESULT(f)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: f
      END FUNCTION ff

      FUNCTION Df(x) RESULT(g)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: g
      END FUNCTION Df

      FUNCTION DDf(x) RESULT(g)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: g
      END FUNCTION DDf

   END INTERFACE


   REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: x0
   REAL(KIND=8),               INTENT(IN)    :: rel_tol
   INTEGER,                    INTENT(INOUT) :: NM_Iter


   REAL(KIND=8), DIMENSION(SIZE(x0)) :: x
   
   REAL(KIND = 8), DIMENSION(SIZE(x0)) :: d_x

   REAL(KIND = 8) :: int_rel_tol, f_1, f_2, j_1, j_2

   INTEGER :: i, iter

   LOGICAL :: forced = .FALSE.
   !-----------------------------------------------------------------------

   WRITE(12,*) '   NEWTON_SYS2_DIAG_JAC'

   ! Verifica della precisione relativa richiesta
   IF (rel_tol > EPSILON(rel_tol)) THEN
      int_rel_tol = rel_tol   
   ELSE
      int_rel_tol = EPSILON(rel_tol)*10
   ENDIF
    
   iter = NM_Iter

   ! Newton  Method       
   x = x0
   WRITE(12,*) '     guess = ', x0

   DO i = 1, iter

       NM_Iter = i

       ! Test over physical meaning of intermediate values for specific 
       ! volume inside newton method
       IF (x(1) <= 1) THEN

          forced = .TRUE.
          x(1) = 1.03

       ENDIF


       IF (x(2) <= 1) THEN

          forced = .TRUE.
          x(2) = 25

       ENDIF


       j_1 = DDf(x(1)) * (x(2) - x(1))
       j_2 = DDf(x(2)) * (x(2) - x(1))


       IF (ABS(j_1) < EPSILON(rel_tol)) THEN

          f_1 = ff(x(1)) + Df(x(1)) * (x(2) - x(1)) - ff(x(2))

          IF (ABS(f_1) < EPSILON(rel_tol)) THEN
              
             d_x(1) = 0

             j_1 = 1.d100

             WRITE(12,*) '     First Diagonal element equal to zero...'
          
          ELSE  

             PRINT*, ''
             PRINT*, 'Error. NEWTON_SYS2_DIAG_JAC:'
             PRINT*, 'Derivative of the  first diagonal  element'
             PRINT*, 'equals to zero, and the right hand side is'
             PRINT*, 'incompatible.'
             PRINT*, ''
             OPEN(UNIT=11, FILE='stopped')
             WRITE(11,*) 'Newton method Diag_Jac: the derivative of the'
             WRITE(11,*) 'first diagonal element equal to zero, and the right'
             WRITE(11,*) 'hand side is incompatible.   STOP'
             CLOSE(11)
             STOP

          ENDIF     

      ENDIF


      IF (ABS(j_2) < EPSILON(rel_tol)) THEN

          f_2 = ff(x(1)) + Df(x(2)) * (x(2) - x(1)) - ff(x(2))

          IF (ABS(f_2) < EPSILON(rel_tol)) THEN
              
             d_x(2) = 0

             j_2 = 1.d100

             WRITE(12,*)  '     Second Diagonal element equal to zero...'
          
          ELSE  

             PRINT*, ''
             PRINT*, 'Error. NEWTON_SYS2_DIAG_JAC:'
             PRINT*, 'Derivative of the  second diagonal element'
             PRINT*, 'equals to zero, and the right hand side is'
             PRINT*, 'incompatible.'
             PRINT*, ''
             OPEN(UNIT=11, FILE='stopped')
             WRITE(11,*) 'Newton method Diag_Jac: the derivative of the'
             WRITE(11,*) 'second diagonal element equal to zero, and the right'
             WRITE(11,*) 'hand side is incompatible.   STOP'
             CLOSE(11)
             STOP

          ENDIF     

       ENDIF


       f_1 = ff(x(1)) + Df(x(1)) * (x(2) - x(1)) - ff(x(2))
       f_2 = ff(x(1)) + Df(x(2)) * (x(2) - x(1)) - ff(x(2))

       d_x(1) = - f_1 / j_1
       d_x(2) = - f_2 / j_2

       x = x + d_x

      IF (SQRT(SUM(d_x**2)) < SQRT(SUM(x**2)) * int_rel_tol)  RETURN
  	
   ENDDO


END FUNCTION  Newton_Sys2_Diag_Jac






FUNCTION  Newton_Sys2_Jac(ff, Df, DDf, x0, rel_tol, NM_Iter)  RESULT (x)
     
   ! Se il parametro OPZIONALE f_shift e' assente, il programma
   ! determina uno zero reale di una funzione reale f(x), ovverosia 
   ! una radice (root) dell'equazione 
   !
   !      f(root) = 0
   !
   ! mediante il metodo di Newton, scritto di forma incrementale.
   !
   ! Se invece il parametro OPZIONALE f_shift e' presente, si
   ! calcola la soluzione root dell'equazione ``shiftata''     
   !
   !      f(root) = f_shift 
   !
   ! secondo un suggerimento di Alberto Guardone. 

   ! ff(x) e' la funzione, Df(x) la sua derivata e 
   ! x0 un valore iniziale che deve essere sufficientemente
   ! vicino a root affinche' il metodo possa convergere.
   ! Quanto vicino, in generale nessuno lo sa a priori e
   ! quindi occorre analizzare ogni problema specifico,
   ! ossia la forma determinata della funzione f(x) in esame.  

   ! La radice e' calcolata con una precisione relativa specificata
   ! tramite il parametro rel_tol. 
   !
   ! Il numero di iterazioni da eseguire e' fornito nel parametro
   ! NM_Iter : in INPUT   NM_Iter e' il numero massimo di iterazioni
   !           in OUTPUT  NM_Iter e' il numero di iterazioni effettuate
   !

   !-----------------------------------------------------------------------
   IMPLICIT NONE
   
   INTERFACE 

      FUNCTION ff(x) RESULT(f)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: f
      END FUNCTION ff

      FUNCTION Df(x) RESULT(g)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: g
      END FUNCTION Df

      FUNCTION DDf(x) RESULT(g)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: g
      END FUNCTION DDf

   END INTERFACE


   REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: x0
   REAL(KIND=8),               INTENT(IN)    :: rel_tol
   INTEGER,                    INTENT(INOUT) :: NM_Iter


   REAL(KIND=8), DIMENSION(SIZE(x0)) :: x
   
   REAL(KIND = 8), DIMENSION(SIZE(x0)) :: d_x

   REAL(KIND = 8) :: int_rel_tol, f1, f2, Df1, Df2, dx,  &
                     j1, j2, j12, det, f_1, f_2 

   INTEGER :: i, iter

   LOGICAL :: forced = .FALSE.
   !-----------------------------------------------------------------------

   WRITE(12,*) '   NEWTON_SYS2_JAC'

   ! Verifica della precisione relativa richiesta
   IF (rel_tol > EPSILON(rel_tol)) THEN
      int_rel_tol = rel_tol   
   ELSE
      int_rel_tol = EPSILON(rel_tol)*10
   ENDIF
    
   iter = NM_Iter

   ! Newton  Method       
   x = x0
   WRITE(12,*) '     guess = ', x0

   DO i = 1, iter

      NM_Iter = i

      ! Test over physical meaning of intermediate values for specific 
      ! volume inside newton method
      IF (x(1) <= 1) THEN

         forced = .TRUE.
         x(1) = 1.03

      ENDIF


      IF (x(2) <= 1) THEN

         forced = .TRUE.
         x(2) = 25

      ENDIF


      f1 = ff(x(1));   Df1 = Df(x(1))
      f2 = ff(x(2));   Df2 = Df(x(2))

      dx = x(2) - x(1)


      j1 = DDf(x(1)) * dx
      j2 = DDf(x(2)) * dx

      j12 = Df1 - Df2 

      det = j1 * j2  -  j12**2  !  Symmetric Jacobian matrix

      IF (ABS(det) < EPSILON(rel_tol)) THEN
         
         PRINT*, ''
         PRINT*, 'Error. NEWTON_SYS2_JAC:'
         PRINT*, 'Determinant of the Jacobian matrix'
         PRINT*, 'less than the rel_tol.   STOP'
         PRINT*, ''
         OPEN(UNIT=11, FILE='stopped')
         WRITE(11,*) 'Newton method Jac: Determinant of the Jacobian'
         WRITE(11,*) 'matrix less than the rel_tol.   STOP'
         CLOSE(11)
         STOP

      ENDIF
      
      f_1 = f1 - f2  +  Df1 * dx  
      f_2 = f1 - f2  +  Df2 * dx  
           
      d_x(1) = - (j2*f_1  -  j12*f_2)/det 
      d_x(2) = - (j1*f_2  -  j12*f_1)/det 

      x = x + d_x

      IF (SQRT(SUM(d_x**2)) < SQRT(SUM(x**2)) * int_rel_tol)  RETURN
  
   ENDDO


END FUNCTION  Newton_Sys2_Jac




FUNCTION  bisection(ff, a, b, rel_tol,  f_shift)  RESULT (root)
   
   ! Se il parametro OPZIONALE f_shift e' assente, il programma
   ! determina uno zero reale di una funzione reale f(x), 
   ! ovverosia una radice (root) dell'equazione 
   !
   !      f(root) = 0
   !
   ! mediante il metodo di bisezione (ricerca dicotomica).
   !
   ! Se invece il parametro OPZIONALE f_shift e' presente, si
   ! calcola la soluzione root dell'equazione ``shiftata''    
   !
   !      f(root) = f_shift 
   !
   ! secondo un suggerimento di Alberto Guardone. 

   ! Lo zero e' ricercato nell'intervallo (reale) avente per estremi
   ! a e b: e' necessario che a e b siano a cavallo di un solo 
   ! zero, ma a e b possono essere tali che a < b oppure a > b.
   !
   ! La radice e' calcolata con una precisione relativa assegnata 
   ! tramite il parametro rel_tol. 
   !
   ! Il numero di iterazioni e' uguale a log2(1/rel_tol).

   !-----------------------------------------------------------------------
   IMPLICIT NONE
   
   INTERFACE 
      FUNCTION ff(x) RESULT(f)
         REAL(KIND=8), INTENT(IN) :: x
         REAL(KIND=8) :: f
      END FUNCTION ff
   END INTERFACE
   
   REAL(KIND=8), INTENT(IN)           :: a, b, rel_tol
   REAL(KIND=8), INTENT(IN), OPTIONAL :: f_shift

   REAL(KIND=8) :: root
   
   REAL(KIND=8) :: xa, xb, xm,  &
                   fa, fb, fm,  &
                   int_rel_tol,  f_sh = 0  
   !-----------------------------------------------------------------------   

   WRITE(12,*) '   BISECTION'

   IF (PRESENT(f_shift)) f_sh = f_shift
   
   fa = ff(a) - f_sh
   fb = ff(b) - f_sh


   ! Controlla se un estremo e' la radice
   IF (ABS(fa) < 1.d-15) THEN
      root = a
      RETURN
   ENDIF

   IF (ABS(fb) < 1.d-15) THEN 
      root = b     
      RETURN
   ENDIF
   
   ! Controlla se la funzione cambia segno nell'intervallo
   IF (fa * fb > 0) THEN

      PRINT*, 'fa * fb = ', fa * fb
      PRINT*, ' Error. BISECTION:'        
      PRINT*, ' No root found in interval [a,b].'
      PRINT*, ' '
      OPEN(UNIT=11, FILE='stopped')
         WRITE(11,*) 'L`intervallo [a,b] non contiene una radice.'
      CLOSE(11)          
      STOP   
   
   ENDIF
      
      
   ! Verifica della precisione relativa richiesta
   IF (rel_tol > EPSILON(xm)) THEN
      int_rel_tol = rel_tol   
   ELSE   
      int_rel_tol = EPSILON(xm)*10
   ENDIF
   
   
   ! Metodo di bisezione   
   xa = a;  xb = b
   
   xm = (xa + xb)/2;  fm = ff(xm) - f_sh
   
   DO WHILE (ABS(xb - xa) > ABS(xm) * int_rel_tol  .AND.  fm /= 0)
      
      IF (fm * fa > 0) THEN
         xa = xm 
      ELSE
         xb = xm
      ENDIF
  
      xm = (xa + xb)/2;  fm = ff(xm) - f_sh
  
   ENDDO
   
   root = xm
   
END FUNCTION  Bisection 





              
SUBROUTINE  ODE_int (x1, x2,  y_start,  PHI)

   !-----------------------------------------------------------------------         
   IMPLICIT NONE                                                                
                                                                             
   REAL (KIND=8), INTENT(IN)    :: x1, x2                                       
   REAL (KIND=8), INTENT(INOUT) :: y_start                                      
                                                                               
   INTERFACE                                                                    
      FUNCTION PHI(x, y)                                                         
         IMPLICIT NONE                                                            
         REAL (KIND=8), INTENT(IN) :: x, y                                        
         REAL (KIND=8) :: PHI                                                     
      END FUNCTION PHI                                                           
   END INTERFACE                                                                
                                                                              
   INTEGER, PARAMETER :: max_step = 10000                                       
                                                                             
   REAL (KIND=8), PARAMETER :: eps = 1.d-9,  h_min = 1.d-8,  tiny = 1.d-30      
                                                                              
   REAL (KIND=8) :: h,  h_did,  h_next,  x,  dydx,  y,  y_scal                  
                                                                               
   INTEGER :: n, n_bad, n_ok                                                    
   !-----------------------------------------------------------------------       

   x = x1                                                                                     
   h = (x2 - x1)/10                                                                           
                                                                                             
   n_ok = 0;  n_bad = 0                                                                       
                                                                                            
   y = y_start                                                                                
                                                                                             
   DO n = 1, max_step                                                                         
                                                                                             
      dydx = PHI(x, y)                                                                        

      y_scal = ABS(y) + ABS(h*dydx) + tiny                                                     

      IF ((x+h-x2)*(x+h-x1) > 0) h = x2 - x                                                    

      CALL RK_qs (y, dydx, x,  h, eps, y_scal,  h_did, h_next,  PHI)                           

      IF (h_did == h) THEN                                                                     
         n_ok = n_ok + 1                                                                        
      ELSE                                                                                     
         n_bad = n_bad + 1                                                                      
      ENDIF                                                                                    
                                                                                             
      IF ((x-x2)*(x2-x1) >= 0) THEN                                                            
                                                                                             
         y_start = y                                                                            
                                                                                             
         RETURN                                                                                 
                                                                                             
      ENDIF                                                                                    
                                                                                           
      IF (ABS(h_next) < h_min) WRITE (12,*) '     Stepsize smaller than minimum in ODE_int'    
                                                                                             
      h = h_next                                                                               
                                                                                           
   ENDDO                                                                                      
                                                                                           
   WRITE (12,*) '     Too many steps in ODE_int'                                              

END SUBROUTINE  ODE_int

     

      
      
SUBROUTINE  RK_qs (y, dydx, x,  h_try, eps, y_scal,  h_did, h_next,  PHI)

   !-----------------------------------------------------------------------         
   IMPLICIT NONE                                                    
                                                                    
   REAL (KIND=8), INTENT(INOUT) :: y,  dydx,  x                     
   REAL (KIND=8), INTENT(OUT)   :: h_did,  h_next                   
   REAL (KIND=8), INTENT(IN)    :: h_try,  eps,  y_scal             
                                                                    
   INTERFACE                                                        
      FUNCTION PHI(x, y)                                             
         IMPLICIT NONE                                                
         REAL (KIND=8), INTENT(IN) :: x, y                            
         REAL (KIND=8) :: PHI                                         
      END FUNCTION PHI                                               
   END INTERFACE                                                    

   REAL (KIND=8), PARAMETER :: safety = 0.9,    grow = -0.2,  &     
                               shrink = -0.25,  err_con = 1.89e-4   
                                                                  
   REAL (KIND=8) :: h,  err_max,  x_new,  y_err,  y_temp            
   !-----------------------------------------------------------------------         

   h = h_try                                                             

   DO WHILE (.TRUE.)                                                     

      CALL RK_ck (y, dydx, x, h,  y_temp, y_err, PHI)                     
                                                                        
      err_max = MAX(0.d0, ABS(y_err/y_scal))/eps                          
                                                                        
      IF (err_max > 1) THEN                                               
                                                                         
         h = safety * h * err_max**shrink                                  
                                                                         
         IF (h < 0.1*h) h = h/10                                           
                                                                         
         x_new = x + h                                                     

         IF (x_new == x) WRITE (12,*) '     Stepsize underflow in RK_qs'   
                                                                          
      ELSE                                                                
                                                                         
         IF (err_max > err_con) THEN                                       
            h_next = safety * h * err_max**grow                             
         ELSE                                                              
            h_next = 5*h                                                    
         ENDIF                                                             
                                                                         
         h_did = h                                                         
         x = x + h                                                         
         y = y_temp                                                        
                                                                         
         RETURN                                                            
                                                                         
      ENDIF                                                               
                                                                         
   ENDDO

END SUBROUTINE  RK_qs

      


      
SUBROUTINE  RK_ck (y, dydx, x, h,  y_out, y_err, PHI)

   !-----------------------------------------------------------------------              
   IMPLICIT NONE                                      

   REAL (KIND=8), INTENT(IN)  :: y, dydx, x, h        
   REAL (KIND=8), INTENT(OUT) :: y_out,  y_err        

   INTERFACE                                          
     FUNCTION PHI(x, y)                               
       IMPLICIT NONE                                  
       REAL (KIND=8), INTENT(IN) :: x, y              
       REAL (KIND=8) :: PHI                           
     END FUNCTION PHI                                 
   END INTERFACE                                      

   REAL (KIND=8) :: ak2, ak3, ak4, ak5, ak6, y_temp   
                                                      
   REAL (KIND=8), PARAMETER :: &                      

   A2  = 0.2d0,  A3 = 0.3d0,  A4 = 0.6d0,  A5 = 1.d0,  A6 = 0.875d0,   &          
                                                                            
   B21 = 0.2d0,  B31 = 3.d0/40.,  B32 = 9.d0/40.,                      &    
   B41 = 0.3d0,  B42 = -0.9d0,  B43 = 1.2d0,  B51 = -11.d0/54.,        &    
   B52 = 2.5d0,  B53 = -70.d0/27.,  B54 = 35.d0/27.,                   &    
   B61 = 1631.d0/55296.,  B62 = 175.d0/512.,                           &    
   B63 = 575.d0/13824.,  B64 = 44275.d0/110592.,  B65 = 253.d0/4096.,  &    
                                                                            
   C1 = 37.d0/378.,  C3 = 250.d0/621.,  C4 = 125.d0/594.,              &    
   C6 = 512.d0/1771.,                                                  &    
                                                                            
   DC1 = C1 - 2825.d0/27648.,   DC3 = C3 - 18575.d0/48384.,            &    
   DC4 = C4 - 13525.d0/55296.,  DC5 = -277.d0/14336.,                  &    
   DC6 = C6 - 0.25d0
   !-----------------------------------------------------------------------         

   y_temp = y + B21*h*dydx                                                 

   ak2 = PHI (x + A2*h,  y_temp)                                           

   y_temp = y + h*(B31*dydx + B32*ak2)                                     

   ak3 = PHI (x + A3*h,  y_temp)                                            

   y_temp = y + h*(B41*dydx + B42*ak2 + B43*ak3)                           

   ak4 = PHI (x + A4*h,  y_temp)                                            

   y_temp = y + h*(B51*dydx + B52*ak2 + B53*ak3 + B54*ak4)                 

   ak5 = PHI (x + A5*h,  y_temp)                                            

   y_temp = y + h*(B61*dydx + B62*ak2 + B63*ak3 + B64*ak4 + B65*ak5)       

   ak6 = PHI (x + A6*h,  y_temp)                                            
                                                                            

   y_out = y + h*(C1*dydx + C3*ak3 + C4*ak4 + C6*ak6)                      

   y_err = h*(DC1*dydx + DC3*ak3 + DC4*ak4 + DC5*ak5 + DC6*ak6)

END SUBROUTINE  RK_ck


END MODULE numerics  
