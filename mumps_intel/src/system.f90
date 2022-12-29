  INCLUDE 'mkl_dss.f90'
MODULE system
  !
  ! Modulo que contiene la resolucion del sistema lineal asociado
  ! usando la biblioteca SuperLU/MKL
  ! 
  !
  USE decimal
  USE types
  USE mkl_dss
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC system_solution, system_solution_by_block
  !
CONTAINS
  !
  SUBROUTINE system_solution(A,bb,method)
    !
    ! Linear system solve driver
    !
    TYPE(t_sparse)        :: A
    REAL(kind=dp)         :: bb(:)
    INTEGER, intent(in)   :: method
    INTEGER               :: precon,ierr
    !
    !
    SELECT CASE(method)       
!!$       !
!!$    CASE(1)
!!$       !
!!$       PRINT*,' '
!!$       PRINT*,' ------------------------------------------------------------------- '
!!$       PRINT*,'              Solving linear system using SuperLU 4.0                '
!!$       PRINT*,' ------------------------------------------------------------------- '
!!$       PRINT*,' '
!!$       !
!!$       ! solving using SuperLU
!!$       !
!!$       CALL superlu(A,bb)
!!$       !
!!$       DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
!!$       IF(ierr/=0) THEN
!!$          PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
!!$          STOP
!!$       END IF
!!$       NULLIFY(A%row,A%column,A%aa)
!!$       !
!!$       !
!!$       PRINT*,' '
!!$       PRINT*,' ------------------------------------------------------------------- '
!!$       !
    CASE(2)
       !
       PRINT*,' '
       PRINT*,' ---------------------------------------------------------------------- '
       PRINT*,'  Solving linear system using GMRES without preconditioning (MKL INTEL) '
       PRINT*,' ---------------------------------------------------------------------- '
       PRINT*,' '
       !
       ! Solving using GMRES without preconditioning
       !
       precon = 0
       !
       CALL gmres(A,bb,precon)
       !
       PRINT*,' '
       PRINT*,' ---------------------------------------------------------------------- '
       !
       ! DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
       ! IF(ierr/=0) THEN
       !    PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
       !    STOP
       ! END IF
       !
    CASE(3)
       !
       ! Solving using GMRES with preconditioning
       !
       PRINT*,' '
       PRINT*,' ------------------------------------------------------------------- '
       PRINT*,'  Solving linear system using GMRES with preconditioning (MKL INTEL) '
       PRINT*,' ------------------------------------------------------------------- '
       PRINT*,' '
       !
       precon = 1
       !
       PRINT*,'dimension sistema',SIZE(bb)
       !
       CALL gmres(A,bb,precon)
       !
       PRINT*,' ------------------------------------------------------------------- '
       PRINT*,' '
       !
       ! DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
       ! IF(ierr/=0) THEN
       !    PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
       !    STOP
       ! END IF
       !
    CASE(4)
       !
       PRINT*,' '
       PRINT*,' ------------------------------------------------------------------- '
       PRINT*,'              Solving linear system using Pardiso              '
       PRINT*,' ------------------------------------------------------------------- '
       PRINT*,' '
       !
       ! solving using SuperLU
       !
       CALL mkl_pardiso(A,bb)
       !
       ! DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
       ! IF(ierr/=0) THEN
       !    PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
       !    STOP
       ! END IF
       ! NULLIFY(A%row,A%column,A%aa)
       ! !
       !
       PRINT*,' '
       PRINT*,' ------------------------------------------------------------------- '
    CASE default
       PRINT*,'ERROR: No method for solving linear system!!'
       STOP
    END SELECT
       !

    !
  END SUBROUTINE system_solution
  !
  !
  SUBROUTINE system_solution_by_block(A,bb,method)
    !
    ! Linear system solve driver, when rhs is a given matrix
    !
    TYPE(t_sparse)      :: A
    REAL(KIND=dp)       :: bb(:,:)
    INTEGER, intent(in) :: method
    INTEGER             :: precon,ierr
    !
    SELECT CASE(method)
!!$       !
!!$    CASE(1)
!!$       !
!!$       PRINT*,' '
!!$       PRINT*,' ------------------------------------------------------------------- '
!!$       PRINT*,'              Solving linear system using SuperLU 4.0                '
!!$       PRINT*,' ------------------------------------------------------------------- '
!!$       PRINT*,' '
!!$       !
!!$       ! solving using SuperLU
!!$       !
!!$       CALL superlu(A,bb)
!!$       !
!!$       DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
!!$       IF(ierr/=0) THEN
!!$          PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
!!$          STOP
!!$       END IF
!!$       NULLIFY(A%row,A%column,A%aa)
!!$       !
!!$       !
!!$       PRINT*,' '
!!$       PRINT*,' ------------------------------------------------------------------- '
!!$       !
    ! CASE(2)
    !    !
    !    PRINT*,' '
    !    PRINT*,' ---------------------------------------------------------------------- '
    !    PRINT*,'  Solving linear system using GMRES without preconditioning (MKL INTEL) '
    !    PRINT*,' ---------------------------------------------------------------------- '
    !    PRINT*,' '
    !    !
    !    ! Solving using GMRES without preconditioning
    !    !
    !    precon = 0
    !    !
    !    CALL gmres(A,bb,precon)
    !    !
    !    PRINT*,' '
    !    PRINT*,' ---------------------------------------------------------------------- '
    !    !
    !    ! DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
    !    ! IF(ierr/=0) THEN
    !    !    PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
    !    !    STOP
    !    ! END IF
    !    !
    ! CASE(3)
    !    !
    !    ! Solving using GMRES with preconditioning
    !    !
    !    PRINT*,' '
    !    PRINT*,' ------------------------------------------------------------------- '
    !    PRINT*,'  Solving linear system using GMRES with preconditioning (MKL INTEL) '
    !    PRINT*,' ------------------------------------------------------------------- '
    !    PRINT*,' '
    !    !
    !    precon = 1
    !    !
    !    PRINT*,'dimension sistema',SIZE(bb)
    !    !
    !    CALL gmres(A,bb,precon)
    !    !
    !    PRINT*,' ------------------------------------------------------------------- '
    !    PRINT*,' '
    !    !
    !    ! DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
    !    ! IF(ierr/=0) THEN
    !    !    PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
    !    !    STOP
    !    ! END IF
       !
    CASE(4)
       !
       PRINT*,' '
       PRINT*,' ------------------------------------------------------------------- '
       PRINT*,'              Solving linear system using Pardiso (rhs a matrix)     '
       PRINT*,' ------------------------------------------------------------------- '
       PRINT*,' '
       !
       !
       CALL mkl_pardiso_by_block(A,bb)
       !
       ! DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
       ! IF(ierr/=0) THEN
       !    PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
       !    STOP
       ! END IF
       ! NULLIFY(A%row,A%column,A%aa)
       ! !
       !
       PRINT*,' '
       PRINT*,' ------------------------------------------------------------------- '
       !
       !
    CASE default
       PRINT*,'ERROR: No method for solving linear system!!'
       STOP
    END SELECT
    !
  END SUBROUTINE system_solution_by_block
  !
  !
!!$  SUBROUTINE superlu(A,bb)
!!$    !
!!$    ! Solving the linear system using SuperLU library (i.e. LU for non-symetric matrices)
!!$    !
!!$    TYPE(t_sparse)   ,INTENT(in)    :: A
!!$    REAL(kind=dp) ,INTENT(inout) :: bb(:)
!!$    INTEGER                      :: nrhs,neq,ldb,info,iopt,factors(8)
!!$    !
!!$    neq     = SIZE(bb)
!!$    ldb     = neq
!!$    nrhs    = 1
!!$    !
!!$    ! First, factorize the matrix. The factors are stored in factor handle.
!!$    !
!!$    iopt = 1
!!$    !
!!$    CALL c_fortran_dgssv(iopt,neq,A%nzero,nrhs,A%aa,A%row,A%column,bb,ldb,factors,info)
!!$    !
!!$    IF (info /= 0) THEN
!!$       WRITE (*,*) 'Problems with the LU factorization (S:superlu M:system) INFO = ',info
!!$       STOP
!!$    ENDIF
!!$    !
!!$    ! Second, solve the system using the existing factors.
!!$    !
!!$    iopt = 2
!!$    !
!!$    CALL c_fortran_dgssv(iopt,neq,A%nzero,nrhs,A%aa,A%row,A%column,bb,ldb,factors,info)
!!$    !
!!$    IF (info /= 0) THEN
!!$       WRITE (*,*) 'Problems with the solution of linear system (S:superlu M:system). INFO = ',info
!!$       STOP
!!$    ENDIF
!!$    !
!!$    ! Last, free the storage allocated inside SuperLU
!!$    !
!!$    iopt = 3
!!$    !
!!$    CALL c_fortran_dgssv(iopt,neq,A%nzero,nrhs,A%aa,A%row,A%column,bb,ldb,factors,info)
!!$    !
!!$    IF (info /= 0) THEN
!!$       WRITE (*,*) 'Problems deallocating arrays of SuperLU (S:superlu M:system). INFO = ',info
!!$    ENDIF
!!$    !
!!$  END SUBROUTINE superlu
!!$  !
  SUBROUTINE gmres(A,bb,precon)
    !
    ! Solve the linear system using GMRES preconditioned or not by ILU
    ! (this is the implementation of GMRES in Intel MKL library)
    !
    !
    TYPE(t_sparse),INTENT(in)   :: A
    REAL(kind=dp),INTENT(inout) :: bb(:)
    INTEGER,INTENT(in)          :: precon 
    !
    INTEGER                     :: ipar(128)
    REAL(kind=dp)               :: dpar(128), xn(SIZE(A%row)-1)
    REAL(kind=dp),ALLOCATABLE   :: tmp(:)
    !
    INTEGER                     :: itercount,ierr,nn,ntmp
    INTEGER                     :: rci_request
    !
    REAL(kind=dp),POINTER       :: bilu0(:),trvec(:)
    !
    nn    = SIZE(bb)
    ntmp  = 25 ! ATTENTION ntmp = ipar(15)
    !
    ntmp = (2*ntmp+1)*nn + (ntmp*(ntmp+9))/2 + 1 !Size of the array tmp
    !
    ALLOCATE(tmp(ntmp))
    !
    xn = 0.0_dp ! intial guest
    !
    PRINT*,'using gmres'
    PRINT*,'neq = ', nn
    !
    !--------------------------------------------------------------------------
    ! initialize the solver
    !--------------------------------------------------------------------------
    !
    CALL dfgmres_init(nn, xn, bb, rci_request, ipar, dpar, tmp)
    IF (rci_request /= 0) THEN
       PRINT*,'Problem with  dfgmres_init (S: gmres,  M: system)!!',rci_request
       STOP
    END IF
    !
    ! If precon=1 then compute the proconditioner
    !
    IF(precon==1) THEN
       !
       ALLOCATE(bilu0(A%nzero),trvec(nn))
       bilu0 = 0.0_dp; trvec = 0.0_dp
       !
       ipar(31) = 1
       dpar(31) = 1.0e-20_dp
       dpar(32) = 1.0e-16_dp
       !
       CALL dcsrilu0(nn, A%aa, A%row, A%column, bilu0, ipar, dpar, ierr)
       IF(ierr/=0) THEN
          PRINT *,' error after calculation of the preconditioner (S: gmres,  M: system) ',ierr
          STOP
       ENDIF
       !
       ipar(11) = 1
       !
    END IF
    !
    ipar(5)  = nn
    ipar(9)  = 1
    ipar(10) = 0
    ipar(12) = 1
    ipar(15) = 25 ! must be the same value than ntmp (first time)
    !
    !--------------------------------------------------------------------------
    ! check the correctness and consistency of the newly set parameters
    !--------------------------------------------------------------------------
    !
    CALL dfgmres_check(nn, xn, bb, rci_request,ipar, dpar, tmp)
    !
    IF (rci_request /= 0) THEN
       PRINT*,'Problem with  dfgmres_check (S: gmres **  M: system) !!',rci_request
       STOP
    END IF
    !
    IF (ipar(8)/=0) THEN
       WRITE(*,'(a,i1,a)') 'as ipar(8)=',ipar(8),', the automatic test for the maximal number of iterations will be'
       PRINT *,'performed'
    ELSE
       WRITE(*,'(a,i1,a)') 'as ipar(8)=',ipar(8),', the automatic test for the maximal number of iterations will be'
       PRINT *,'skipped'
    ENDIF
    !
    PRINT *,'+++'
    IF (ipar(9)/=0) THEN
       WRITE(*,'(a,i1,a)') 'as ipar(9)=',ipar(9),', the automatic residual test will be performed'
    ELSE
       WRITE(*,'(a,i1,a)') 'as ipar(9)=',ipar(9),', the automatic residual test will be skipped'
    ENDIF
    PRINT *,'+++'
    !
    PRINT *,'+++'
    IF (ipar(11)/=0) THEN
       WRITE(*,'(a,i1,a)') 'as ipar(11)=',ipar(11),', the preconditioned fgmres iterations will be performed, thus,'
       PRINT *,'the preconditioner action will be requested via rci_request=3'
    ELSE
       WRITE(*,'(a,i1,a)') 'as ipar(11)=',ipar(11),', the preconditioned fgmres iterations will not be performed,'
       PRINT *,'thus, rci_request will not take the value 3'
    ENDIF
    PRINT *,'+++'
    !
    IF (ipar(12)/=0) THEN
       WRITE(*,'(a,i1,a)')'as ipar(12)=',ipar(12),', the automatic test for the norm of the next generated vector is'
       PRINT *,'not equal to zero up to rounding and computational errors will be performed,'
       PRINT *,'thus, rci_request will not take the value 4'
    ELSE
       WRITE(*,'(a,i1,a)')'as ipar(12)=',ipar(12),', the automatic test for the norm of the next generated vector is'
       PRINT *,'not equal to zero up to rounding and computational errors will be skipped,'
       PRINT *,'thus, the user-defined test will be requested via rci_request=4'
    ENDIF
    PRINT *,'+++'
    !
    !
    main_cycle: DO 
       !
       ! Main iteration cycle of GMRES
       !
       CALL dfgmres(nn, xn, bb, rci_request, ipar, dpar, tmp)
       !
       SELECT CASE(rci_request)
       CASE(0)
          EXIT main_cycle
       CASE(1)
          !
          ! Only if we dont have preconditioning
          !
          CALL mkl_dcsrgemv('n',nn, A%aa, A%row, A%column, tmp(ipar(22)), tmp(ipar(23)))
          CYCLE main_cycle
          !
       CASE(2)
          !
          CYCLE main_cycle
          !
       CASE(3)
          !
          ! only if we use a preconditioner
          !
          CALL mkl_dcsrtrsv('L','N','U',nn,bilu0,A%row,A%column,tmp(ipar(22)),trvec)
          CALL mkl_dcsrtrsv('U','N','N',nn,bilu0,A%row,A%column,trvec,tmp(ipar(23)))
          CYCLE main_cycle
          !
       CASE(4)
          !
          IF(dpar(7) < 1.0e-12_dp) EXIT main_cycle
          CYCLE main_cycle
          !
       CASE default
          !
          PRINT*,'ERROR: main cycle, rci_request (S: gmres,  M: system) ', rci_request
          STOP
          !
       END SELECT
       !
    END DO main_cycle
    !
    ipar(13) = 0
    !
    CALL dfgmres_get(nn, xn, bb, rci_request, ipar, dpar, tmp, itercount)   
    !
    bb = xn
    !
    PRINT *, ''
    PRINT *,' number of iterations:      ',itercount
    PRINT *, ''
    !
    PRINT *, ''
    PRINT *,' norm of the last residual: ',dpar(5)
    PRINT *, ''
    !
    DEALLOCATE(tmp,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problems in  dealocation of tmp (S: gmres, M: system)!!'
       STOP
    END IF
    !
    IF(ASSOCIATED(bilu0)) THEN
       DEALLOCATE(bilu0,trvec,STAT=ierr)
       IF(ierr/=0) THEN
          PRINT*,'Problems in  dealocation of bilu0 and trvec (S: gmres, M: system)!!'
          STOP
       END IF
    END IF
    !
  END SUBROUTINE gmres
  !
  !
  SUBROUTINE mkl_pardiso(A,bb)
    !
    ! Solving the linear system using Pardiso solver (MKL library)
    !
    TYPE(t_sparse)               :: A
    REAL(kind=dp) ,INTENT(inout) :: bb(:)
    REAL(kind=dp),ALLOCATABLE    :: solution(:) 
    INTEGER, PARAMETER           :: bufLen = 20
    INTEGER                      :: error,nCols, nRows, nRhs,nn
    TYPE(MKL_DSS_HANDLE)         :: handle 
    REAL(kind=dp),ALLOCATABLE    :: statOUt( : )
    CHARACTER(len=15)            :: statIn
    INTEGER                      :: perm(1), ierr
    INTEGER                      :: buff(bufLen)
    !
    nn      = SIZE(bb)
    !
    nRows   = nn
    nCols   = nn
    nRhs    = 1
    perm(1) = 0
    ierr    = 0
    !
    ALLOCATE(solution(nn),STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problems in  alocation of vector solution (S: mkl_pardiso, M: system)!!'
       STOP
    END IF
    solution = 0.0_dp
    !
    ! Initialize the solver.
    error = DSS_CREATE( handle, MKL_DSS_DEFAULTS )
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
    !
    ! Define the non-zero structure of the matrix.
    error = DSS_DEFINE_STRUCTURE( handle, MKL_DSS_NON_SYMMETRIC, A%row, nRows,nCols, A%column, A%nzero)
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
    !
    ! Reorder the matrix.
    !    
    !error = DSS_REORDER( handle, MKL_DSS_METIS_OPENMP_ORDER, perm )
    error = DSS_REORDER( handle, MKL_DSS_AUTO_ORDER, perm )
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
    !
    ! Factor the matrix.
    !
    error = DSS_FACTOR_REAL( handle, MKL_DSS_INDEFINITE, A%aa )
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
    !
    ! Solving the linear system
    !
    error = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS, bb, nRhs, solution )    
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
    !
    bb = solution !the solution exits the soubroutine trough bb
    !
    DEALLOCATE(solution,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problems in  dealocation of vector solution (S: mkl_pardiso, M: system)!!'
       STOP
    END IF   
    !
    ! Print Out the quantity of memory used
    !
    ALLOCATE(statOut( 5 ) )
    statIn = 'Solvemem' !'determinant'
    CALL mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)    
    !
    error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut )
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
    !
    WRITE(*,'(a,1x,f8.2,a)')  'Memoria Total Utilizada = ',statOut(1)/1024.0_dp, ' (Mb)' 
    !    
    ! Deallocate solver storage and various local arrays.
    !
    error = DSS_DELETE( handle, MKL_DSS_DEFAULTS )
    IF (error /= MKL_DSS_SUCCESS ) GOTO 999    
    !
    GOTO 1000    !  en este punto guatea (ramiro)
    ! Print an error message and exit
999 WRITE(*,*) "Solver returned error code ", error    
1000 CONTINUE
    !
  END SUBROUTINE mkl_pardiso
  !
  !
  SUBROUTINE mkl_pardiso_by_block(A,b)
    !
    ! Misma idea de mkl_pardiso pero con lado derecho multiplce (una matriz)
    !
    INTEGER                          :: iparm(64),maxfct, mnum, mtype, phase, &
         n, nrhs, error, msglvl,pt(64),ierr
    INTEGER                          :: i,idum(1),ddum(1)
    REAL(kind=dp),     ALLOCATABLE   :: x(:,:)
    !REAL(kind=dp),    POINTER       :: x(:,:)
    REAL(kind=dp),     INTENT(inout) :: b(:,:)  
    TYPE(t_sparse)                   :: A
    !
    n      = SIZE(b,1)
    nrhs   = SIZE(b,2)
    maxfct = 1
    mnum   = 1
    !
    !NULLIFY(x)
    PRINT*,'allocate'
    ALLOCATE(x(n,nrhs))
    !
    PRINT*,'x'
    x = 0.0_dp
    PRINT*,'x (end)'
    !
!!$    DATA n /5/, nrhs /3/, maxfct /1/, mnum /1/
!!$    DATA A%row /1,4,6,9,12,14/
!!$    DATA A%column
!!$1   /   1,    2,          4,
!!$2     1 ,    2,
!!$3                 3 ,    4,    5,
!!$4     1 ,          3,    4,
!!$5           2 ,                5/
!!$    DATA a
!!$1   /1.d0,-1.d0,      -3.d0,
!!$2   -2.d0, 5.d0,
!!$3              4 .d0, 6.d0, 4.d0,
!!$4   -4.d0,       2.d0, 7.d0,
!!$5        8 .d0,            -5.d0/
    !
    !C..
    !C.. Set up PARDISO control parameter
    !C..
    DO i = 1, 64
       iparm(i) = 0
    END DO
    iparm(1) = 1 ! no solver default
    iparm(2) = 2 ! fill-in reordering from METIS
    iparm(3) = 1 ! numbers of processors
    iparm(4) = 0 ! no iterative-direct algorithm
    iparm(5) = 0 ! no user fill-in reducing permutation
    iparm(6) = 0 ! =0 solution on the first n compoments of x
    iparm(7) = 0 ! not in use
    iparm(8) = 9 ! numbers of iterative refinement steps
    iparm(9) = 0 ! not in use
    iparm(10) = 13 ! perturbe the pivot elements with 1E-13
    iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
    iparm(12) = 0 ! not in use
    iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
    iparm(14) = 0 ! Output: number of perturbed pivots
    iparm(15) = 0 ! not in use
    iparm(16) = 0 ! not in use
    iparm(17) = 0 ! not in use
    iparm(18) = -1 ! Output: number of nonzeros in the factor LU
    iparm(19) = -1 ! Output: Mflops for LU factorization
    iparm(20) = 0 ! Output: Numbers of CG Iterations
    error = 0 ! initialize error flag
    msglvl = 1 ! print statistical information
    mtype = 11 ! real unsymmetric
    !C.. Initiliaze the internal solver memory pointer. This is only
    !C necessary for the FIRST call of the PARDISO solver.
    DO i = 1, 64
       pt(i) = 0
    END DO
    !C.. Reordering and Symbolic Factorization, This step also allocates
    !C all memory that is necessary for the factorization
    phase = 11 ! only reordering and symbolic factorization
    PRINT*,'1'
    CALL pardiso_64(pt, maxfct, mnum, mtype, phase, n, A%aa, A%row, A%column,idum, nrhs, iparm, msglvl, ddum, ddum, error)
    PRINT*,'2'
    WRITE(*,*) 'Reordering completed ... '
    IF (error .NE. 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
       STOP 1
    END IF
    WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
    WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
    !C.. Factorization.
    phase = 22 ! only factorization
    PRINT*,'3'
    CALL pardiso_64(pt, maxfct, mnum, mtype, phase, n, A%aa, A%row, A%column,idum, nrhs, iparm, msglvl, ddum, ddum, error)
    PRINT*,'4'
    WRITE(*,*) 'Factorization completed ... '
    IF (error .NE. 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
       STOP 1
    ENDIF
    !C.. Back substitution and iterative refinement
    iparm(8) = 2 ! max numbers of iterative refinement steps
    phase = 33 ! only factorization
!!$    do i = 1, n
!!$       b(i,1) = 1.d0
!!$    end do
!!$    write(*,*) ' '
!!$    do i = 1, n
!!$       b(i,2) = 2.d0
!!$    end do
!!$    write(*,*) ' '
!!$    do i = 1, n
!!$       b(i,3) = 3.d0
!!$    end do
    PRINT*,'5'
    CALL pardiso_64(pt, maxfct, mnum, mtype, phase, n, A%aa, A%row, A%column,idum, nrhs, iparm, msglvl, b, x, error)
    PRINT*,'6'
    WRITE(*,*) 'Solve completed ... '
    WRITE(*,*) 'The solution of the system is '
!!$    DO i = 1, n
!!$       WRITE(*,*) ' x(',i,',1) = ', x(i,1)
!!$    END DO
!!$    DO i = 1, n
!!$       WRITE(*,*) ' x(',i,',2) = ', x(i,2)
!!$    END DO
!!$    DO i = 1, n
!!$       WRITE(*,*) ' x(',i,',3) = ', x(i,3)
!!$    END DO
    !C.. Termination and release of memory
    phase = -1 ! release internal memory
    PRINT*,'pardiso  ***'
    CALL pardiso_64(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,idum, nrhs, iparm, msglvl, ddum, ddum, error)
    PRINT*,'pardiso *** (end)'
    !
    b = 0.0_dp
    PRINT*,'size(b,1)=',SIZE(b,1),'size(b,2)=',SIZE(b,2)
    PRINT*,'size(x,1)=',SIZE(x,1),'size(x,2)=',SIZE(x,2)
    !print*,'x=',x
    !b = x
    DO i=1,SIZE(b,1)
       b(i,:) = x(i,:)
    END DO
    !
!    PRINT*,'deallocando A'
!    !
!    DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
!    IF(ierr/=0) THEN
!       PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!!'
!       STOP
!    END IF
!    NULLIFY(A%row,A%column,A%aa)
    !
    DEALLOCATE(x)
    ! NULLIFY(x)
    !
    PRINT*,' '
    PRINT*,' ------------------------------------------------------------------- '
    !
    !
  END SUBROUTINE mkl_pardiso_by_block
  !
  !
END MODULE system
