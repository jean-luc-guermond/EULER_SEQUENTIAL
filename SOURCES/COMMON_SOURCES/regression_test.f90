MODULE regression_test
  USE input_data
CONTAINS
  SUBROUTINE regression(un)
    USE mesh_handling
    USE fem_tn
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(:), POINTER :: curr_reg_test, ref_reg_test
    INTEGER :: it, it_max
    INTEGER,      PARAMETER :: in_unit=10, out_unit=11
    REAL(KIND=8), PARAMETER :: tol = 1.d-12
    IF (SIZE(un,1)==mesh%np) THEN
       it_max = SIZE(un,2)
    ELSE
       it_max = SIZE(un,1)
    END IF
    ALLOCATE(curr_reg_test(it_max), ref_reg_test(it_max))
    OPEN(UNIT = out_unit, FILE = 'current_regression_test', FORM = 'formatted', STATUS = 'unknown')
    DO it = 1, it_max
       IF (SIZE(un,1)==mesh%np) THEN
          CALL ns_0 (mesh, un(:,it), curr_reg_test(it))
       ELSE
          CALL ns_0 (mesh, un(it,:), curr_reg_test(it))
       END IF
       WRITE(out_unit,*) curr_reg_test(it)
    END DO
    CLOSE(out_unit)
    CALL check_regression(it_max,curr_reg_test)
  END SUBROUTINE regression

  SUBROUTINE check_regression(it_max,curr_reg_test)
    IMPLICIT NONE
    CHARACTER(LEN=5) :: err_out
    REAL(KIND=8), DIMENSION(it_max), INTENT(IN) :: curr_reg_test
    REAL(KIND=8), DIMENSION(it_max)             :: ref_reg_test
    INTEGER :: it, it_max
    INTEGER,      PARAMETER :: in_unit=10, out_unit=11
    REAL(KIND=8), PARAMETER :: tol = 1.d-12
    REAL(KIND=8) :: err
    IF (inputs%if_regression_test) THEN
       OPEN(UNIT = in_unit, FILE = 'reference_regression_test', FORM = 'formatted', STATUS = 'unknown')
       DO it = 1, it_max
          READ(in_unit,*) ref_reg_test(it)
       END DO
       CLOSE(in_unit)
       err = MAXVAL(ABS(curr_reg_test-Ref_reg_test))/MAXVAL(ABS(ref_reg_test))
       IF (err.GT.tol) THEN
          err_out = 'fail' !===Test failed
          DO it = 1, it_max
             WRITE(*,*) ref_reg_test(it), curr_reg_test(it)
          END DO
       ELSE
          err_out = 'ok'   !===Test passed
       END IF
       STOP err_out
    ELSE
       RETURN
    END IF
  END SUBROUTINE check_regression

  
END MODULE regression_test
