MODULE decimal
  !
  ! Definition of the precision of the real variables
  ! dp = double  precision
  ! sp = simple  precision
  !
  IMPLICIT NONE
  INTEGER, PARAMETER::dp=KIND(1.D0)
  INTEGER, PARAMETER::sp=KIND(1.0)
  !
END MODULE decimal
