PROGRAM main
  !
  ! Programa para resolver sistema lineal con SuperLU
  !
  ! Creacion           : 11/08/2017
  ! Ultima modificacion: ninguna
  !
  USE decimal       ! declaracion de la precision
  USE types         ! definicion de los diferentes tipos a usar
  USE mumps_system  ! MUMPS
  USE system        ! MKL
  !
  IMPLICIT NONE
  !
  REAL(KIND=dp), ALLOCATABLE :: Afull(:,:),rhs(:),rhs_block(:,:),IPIV(:)
  INTEGER                    :: info
  CHARACTER(LEN=32)          :: fout
  TYPE(t_sparse)             :: A
  !
  ! Test 1:   [1 0 3;0 1 1;0 0 1]\[1;1;2] = [-5;-1;2]
  PRINT*,'Test 1:'
  PRINT*,'        [1 0 3]'
  PRINT*,'  A  =  [0 1 1] ;   rhs = [1 1 2];   exact solution: [-5,-1,2]'
  PRINT*,'        [0 0 1]'
  A%nzero = 5
  NULLIFY( A%aa, A%column, A%row )
  ALLOCATE( A%aa(A%nzero), A%column(A%nzero), A%row(4), rhs(3) )
  !
  A%aa = 0.0_dp;  A%column = 0;  A%row = 0;  rhs = 0.0_dp
  !
  A%aa     = (/ 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 1.0_dp /)
  A%row    = (/ 1, 3, 5, 6 /)
  A%column = (/ 1, 3, 2, 3, 3 /)
  !
  rhs      = (/ 1.0_dp, 1.0_dp, 2.0_dp /)
  !
  CALL system_mumps(A,rhs)   ! solving the sparse matrix with MUMPS
  !CALL system_solution(A,rhs,4) ! 4 is Pardiso (Intel MKL version)
  !
  PRINT*,' '
  PRINT*,' --------------------------------------'
  PRINT*,' Solving the linear system using MUMPS or Pardiso (Intel MKL version)'
  PRINT*,' '
  PRINT*,rhs
  PRINT*,' '
  !
  ALLOCATE( Afull(3,3), rhs_block(3,1) )
  Afull = 0.0_dp
  Afull(1,1) = 1.0_dp;  Afull(1,3) = 3.0_dp
  Afull(2,2) = 1.0_dp;  Afull(2,3) = 1.0_dp
  Afull(3,3) = 1.0_dp
  !
  rhs_block = 1.0_dp;  rhs_block(3,1) = 2.0_dp
  !
  ! Solving the same linear system but with a full matrix 
  ALLOCATE( IPIV(3) )
  CALL dgesv(3,3,Afull,3,IPIV,rhs_block,3,info) ! LU with pivot
  !
  PRINT*,'Solving the linear system with a full matrix (using LU from MKL)'
  PRINT*, rhs_block
  !
  DEALLOCATE( rhs, Afull, rhs_block )
  !
END PROGRAM main
