MODULE post_proc
  USE lin_alg_tools
  USE matrix_type
  USE space_dim
  USE def_type_mesh
  PUBLIC :: schlieren
  PRIVATE
  TYPE(matrice_bloc), DIMENSION(k_dim), SAVE :: cij
  TYPE(matrice_bloc) :: mass
  REAL(KIND=8), DIMENSION(:), POINTER ::lumped
CONTAINS
  SUBROUTINE schlieren(mesh,un,grad)
    implicit none
    TYPE(mesh_type)                              :: mesh
    REAL(KIND=8), DIMENSION(mesh%np)     :: un
    REAL(KIND=8), DIMENSION(mesh%np)     :: grad
    REAL(KIND=8), DIMENSION(k_dim) :: xx
    REAL(KIND=8) :: mx, mn
    integer :: i, j, k, p
    LOGICAL :: once=.true.
    IF (once) THEN
       CALL compute_mass(mesh,mass)
       CALL compute_cij(mesh,cij)
       CALL lumped_mass(mesh,mass,lumped)
       once=.false.
    END IF
    !===Gradient
    DO i = 1, mesh%np
       xx = 0.d0
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          DO k = 1, k_dim
             xx(k) = xx(k) + un(j)*cij(k)%aa(p)
          END DO
       END DO
       grad(i) = SQRT(SUM(xx**2))/lumped(i)
    END DO
    mx=MAXVAL(grad)
    mn=MINVAL(grad)
    grad = EXP(-10*(grad-mn)/(mx-mn))
  END SUBROUTINE schlieren
END MODULE post_proc
