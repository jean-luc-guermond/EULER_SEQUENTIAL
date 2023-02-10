program extract
  IMPLICIT NONE
  REAL(KIND=8) :: err1, err0
  CHARACTER(LEN=40) :: fich, stuff1, stuff2, stuff3, stuff4
  INTEGER :: nb_fich, i, k, i0, i1
  open(unit=10,file='data_file',form='formatted',status='unknown')
  OPEN(unit=12,file='error',form='formatted',status='unknown')
  READ(10,*) nb_fich
  DO i = 1, nb_fich
     READ(10,*) fich
     OPEN(unit=11,file=fich,form='formatted',status='unknown')
     DO k= 1, 1 !===read L1 norm
        READ(11,*) stuff1
     END DO
     READ(11,*) stuff1, stuff2, stuff3, stuff4, err1
     CLOSE(11)
     IF (i==1) THEN
        WRITE(12,*) 100*2**(i-1), err1
     ELSE
        backspace(12)
        read(12,*) i0, err0
        write(*,*) i0, err0
        i1 = 2*i0
        WRITE(12,*) i1, err1, log(err0/err1)/log(i1/float(i0)) 
     END IF
  END DO
  
end program extract
