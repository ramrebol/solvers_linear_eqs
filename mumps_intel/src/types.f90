MODULE types
  !
  !  Definition of new types for some variables
  !
  USE decimal
  !
  IMPLICIT NONE
  !
  TYPE t_sparse
     !
     ! Morse storage (by rows) of a matrix
     !
     ! nzero       : number of the non-zero elements of the matrix
     ! row, column : pointers to the rows and columns
     ! aa          : array with the non-zeros elements of the matrix
     !
     INTEGER                :: nzero
     INTEGER,       POINTER :: row(:),column(:) 
     REAL(KIND=dp), POINTER :: aa(:)         
  END TYPE t_sparse
  !
END MODULE types
