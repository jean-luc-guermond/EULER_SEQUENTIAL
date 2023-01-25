MODULE  grid_1d

!  |                           i                                           |
!  o----o--------0-----O================O--------0----o--------------------o
!  |             n3    n1               n2       n4                        |


!  alternatively and equivalently


!  |                           i                                           |
!  o----o--------0-----O================O--------0----o--------------------o
!  |             n4    n2               n1       n3                        |




!  |  i                                                                    |
!  O====O--------0-----o----------------o--------o----o--------------------o
!n1=n3  n2       n3                                                        |



!  |                                                             i         |
!  o---o--------o-----o----------------o---------0----O====================O
!  |                                             n3   n1                 n2=n4

   IMPLICIT NONE

   PRIVATE
   
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE, PUBLIC :: xx, cell_size
   INTEGER,      DIMENSION(:,:), ALLOCATABLE, PUBLIC :: node_pair

   PUBLIC :: gen_grid_1d 

CONTAINS
!=======

  SUBROUTINE  gen_grid_1d (a, b, np,  uniform)

    IMPLICIT NONE

    REAL(KIND=8),      INTENT(IN) :: a, b
    INTEGER,           INTENT(IN) :: np
    LOGICAL, OPTIONAL, INTENT(IN) :: uniform 

    REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: harvest
    INTEGER :: n, i, n1, n2

    IF (ALLOCATED(xx)) THEN
       DEALLOCATE(xx, cell_size, node_pair)
    END IF
    ALLOCATE (xx(np), cell_size(np), node_pair(4,np-1))

    DO n = 1, np
       xx(n) = a + (n-1)*((b-a)/(np-1))
    ENDDO

    DO i = 1, np-1
       node_pair(1,i) = i;  node_pair(2,i) = i + 1
    ENDDO

    CALL  set_horns (xx, node_pair)

    IF (PRESENT(uniform)) THEN

       IF (.NOT. uniform) THEN
          ! NONuniform grid
          ALLOCATE (harvest(np))     
          CALL  random_number(harvest)
          xx = xx + (harvest - 0.5)/(1.75*np)
          DEALLOCATE (harvest)
       ENDIF

    ENDIF

    cell_size = 0

    DO i = 1, np-1   
       n1 = node_pair(1,i);   n2 = node_pair(2,i)
       cell_size(n1) = cell_size(n1) + ABS(xx(n2) - xx(n1))/2
       cell_size(n2) = cell_size(n2) + ABS(xx(n2) - xx(n1))/2
    ENDDO


  END SUBROUTINE  gen_grid_1d

!=========================================================================

SUBROUTINE  set_horns (xx, node_pair)

!  Extension of the node-pair structure node_pair(3:4,:),
!  form the abscissas of the nodes, fully arbitrarily ordered,
!  and from the node-pair connectivity node_pair(1:2,:)

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: xx
   INTEGER,      DIMENSION(:,:), INTENT(INOUT) :: node_pair

   INTEGER :: i, np, n1, n2

   np = SIZE(xx)

   DO i = 1, SIZE(node_pair, 2)

      n1 = node_pair(1,i)
      n2 = node_pair(2,i)

      IF (xx(n2) > xx(n1)) THEN
         node_pair(3:3,i) = MAXLOC( xx,  xx < SPREAD(xx(n1), 1, np) )
         node_pair(4:4,i) = MINLOC( xx,  xx > SPREAD(xx(n2), 1, np) )
      ELSE
         node_pair(3:3,i) = MINLOC( xx,  xx > SPREAD(xx(n1), 1, np) )
         node_pair(4:4,i) = MAXLOC( xx,  xx < SPREAD(xx(n2), 1, np) )
      ENDIF

      IF (node_pair(3,i) == 0)  node_pair(3,i) = node_pair(1,i)
      IF (node_pair(4,i) == 0)  node_pair(4,i) = node_pair(2,i)

      ! WRITE (*,*) node_pair(:,i)

   ENDDO

   WRITE (*,*)

END SUBROUTINE  set_horns


END MODULE  grid_1d
