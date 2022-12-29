MODULE mumps_system
  !
  USE decimal
  USE types
  !
  IMPLICIT NONE
  INCLUDE 'dmumps_struc.h'
  !
  PRIVATE
  !
  PUBLIC system_mumps, system_mumps_by_blocks
  !
CONTAINS
  !
  SUBROUTINE system_mumps(A,bb)
    !
    ! Linear system solve driver
    !
    TYPE(t_sparse)               :: A
    REAL(kind=dp),INTENT(inout)  :: bb(:)
    INTEGER                      :: precon,ierr
    TYPE (DMUMPS_STRUC)          :: mumps_par
    INTEGER                      :: I,k,j,iunit
    !
    PRINT*,'  Initialize an instance of the package'
    !
    mumps_par%JOB = -1
    mumps_par%SYM = 0
    mumps_par%PAR = 1
    !
    CALL DMUMPS(mumps_par)
    !
    mumps_par%N  = SIZE(bb)
    mumps_par%NZ = A % nzero
    !
    !************************************
    ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
    ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
    ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )
    ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
    !
    k = 1
    !
    DO i=1,mumps_par%N
       !
       DO j=A%row(i),A%row(i+1)-1
          !
          mumps_par%IRN(k) = i
          !
          k= k + 1
       END DO
       !
    END DO
    !
    ! mumps_par%JCN = A%column ! falla esto si nzero es grande
    !
    DO i=1,A%nzero
       !
       mumps_par%JCN(i) = A%column(i)
       mumps_par%A(i)   = A%aa(i)
       !
    END DO
    !
    mumps_par%RHS = bb
    !
    !************************************
    !
    DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problems in deallocating matrix A (S: system_mumps M:system)!!'
       STOP
    END IF
    NULLIFY(A%row,A%column,A%aa)
    !
    !************************************
    !
    PRINT*,' '
    PRINT*,' ---------------------------------------------------------------------------- '
    PRINT*,'              Solving linear system using MUMPS 5.5.1            '
    PRINT*,' ---------------------------------------------------------------------------- '
    PRINT*,' '
    !
    mumps_par%ICNTL(28)  = 2 ! determines whether a sequential or a parallel analysis is performed ! 0: automatic choice
    !
    mumps_par%ICNTL(7)  = 7  ! ordenamiento ver manual
    !
    mumps_par% icntl(22) = 0 ! 0 in-core 1 out-of-core  !! Prueba 28 May 2015
    !
    PRINT*,'************ ANALISIS ***********'
    !
    mumps_par%JOB = 1 ! analisis
    !
    CALL DMUMPS(mumps_par)
    !
    PRINT*,'************ FACTORIZA ***********'
    !
    mumps_par%JOB = 2 ! factoriza
    !
    mumps_par%ICNTL(14) = 70 
    !
    CALL DMUMPS(mumps_par)
    !
    PRINT*,'************ RESUELVE ***********'
    !
    mumps_par%ICNTL(10) = -10 ! Refinamiento iterativo =0 por defecto, <0 no hace nada, mumps_par%ICNTL(10) = numero maximo de refinamiento iterativo
    !
    mumps_par%ICNTL(11) = 1 ! 2=compute main statistics =1 Full statistics
    !
    mumps_par%JOB = 3 ! resuelve
    !
    CALL DMUMPS(mumps_par)
    !
    bb = 0.0_dp
    !
    bb = mumps_par%RHS
    !
    !  Deallocate user data
    !
    IF ( mumps_par%MYID .EQ. 0 )THEN
       DEALLOCATE( mumps_par%IRN )
       DEALLOCATE( mumps_par%JCN )
       DEALLOCATE( mumps_par%A   )
       DEALLOCATE( mumps_par%RHS )
    END IF
    !
    iunit = 10
    !
    OPEN(iunit,file='info_esqueleto.dat',status='replace', action='write' )
    !
    !WRITE(iunit,*) "dimension esqueleto & Upper bound ERROR         CONDITION NUMBER (1)      CONDITION NUMBER (2) "
    WRITE(iunit,'(1x,i8,5x,e25.15,1x,e25.15,1x,e25.15)') mumps_par%N,mumps_par%RINFOG(9),mumps_par%RINFOG(10),mumps_par%RINFOG(11)
    !
    CLOSE(iunit)
    !
    ! Destroy the instance (deallocate internal data structures)
    !
    mumps_par%JOB = -2
    !
    CALL DMUMPS(mumps_par)
    !
  END SUBROUTINE system_mumps
  !
  !
  SUBROUTINE system_mumps_by_blocks(A,bb)
    !
    ! Linear system solve driver
    !
    IMPLICIT NONE
    !INCLUDE 'mpif.h'
    INCLUDE 'dmumps_struc.h'
    TYPE(t_sparse)               :: A
    REAL(kind=dp),INTENT(inout)  :: bb(:,:)
    INTEGER                      :: precon,ierr
    TYPE (DMUMPS_STRUC)          :: mumps_par
    INTEGER                      :: I,k,j
    !
    ierr = 0
    !
    PRINT*,'NZERO=',A%nzero
    !
    !CALL MPI_INIT(IERR)
    !
    PRINT*,'Define a communicator for the package.'
    ! Define a communicator for the package.
    !
    ! mumps_par%COMM = MPI_COMM_WORLD
    !
    !  Initialize an instance of the package
    !  for L U factorization (sym = 0, with working host)
    !
    PRINT*,'  Initialize an instance of the package'
    !
    mumps_par%JOB = -1
    mumps_par%SYM = 0
    mumps_par%PAR = 1
    !
    CALL DMUMPS(mumps_par)
    !
    PRINT*,'block JOB=-1'
    !pause
    !
    !  Define problem on the host (processor 0)
    !
    PRINT*,'  Define problem on the host (processor 0)'
    !
    mumps_par%N  = SIZE(bb,1)
    mumps_par%NZ = A % nzero
    !
    mumps_par%ICNTL(20) = 0 ! Para definir el lado derecho como una matriz
    mumps_par%ICNTL(21) = 0 ! Para definir el lado derecho como una matriz
    mumps_par%LRHS = SIZE(bb,1)
    mumps_par%NRHS = SIZE(bb,2)
    !
    PRINT*,'MUMPS_PAR%IRN'
    ! mumps_par%IRN
    !
    !************************************
    ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
    ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
    ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )
    NULLIFY(mumps_par%RHS)
    ALLOCATE( mumps_par%RHS (mumps_par%LRHS*mumps_par%NRHS))
    !
    !
    k = 1
    !
    DO i=1,mumps_par%N
       !
       DO j=A%row(i),A%row(i+1)-1
          !
          ! PRINT*,'k=',k
          !
          mumps_par%IRN(k) = i
          !
          k= k + 1
          !PRINT*,'k=',k
          !
       END DO
       !
    END DO
    !
    PRINT*,'listo 1'
    PRINT*,'k=',k
    PRINT*,'nzero=',A%nzero
    !stop 
    !
    !
    mumps_par%JCN = A%column
    !
    mumps_par%A = A%aa
    !
    !mumps_par%RHS = bb
    !
    i=0
    !
    DO j=1,mumps_par%NRHS
       !
       i = i + 1
       !
       mumps_par%RHS((j-1)*mumps_par%N+1:j*mumps_par%N) = bb(:,i)
       !
    END DO
    !
    PRINT*,'listo 2'
    !
    !************************************
    !
    DEALLOCATE(A%row,A%column,A%aa,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problems in deallocating matrix A (S: system_solution M:system)!!'
       STOP
    END IF
    NULLIFY(A%row,A%column,A%aa)
    !
    !
    PRINT*,' '
    PRINT*,' ------------------------------------------------------------------- '
    !
    !************************************
    !
    !
!!$      IF ( mumps_par%MYID .EQ. 0 ) THEN
!!$         READ(5,*)  mumps_par%NZ6!!$         READ(5,*) mumps_par%NZ
!!$         ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
!!$         ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
!!$         ALLOCATE( mumps_par%A( mumps_par%NZ ) )
!!$         ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
!!$         DO I = 1, mumps_par%NZ
!!$            READ(5,*) mumps_par%IRN(I),mumps_par%JCN(I),mumps_par%A(I)
!!$         END DO
!!$         DO I = 1, mumps_par%N
!!$            READ(5,*) mumps_par%RHS(I)
!!$         END DO
!!$      END IF
    !
    !  Call package for solution
    !
    !PRINT*,' Call package for solution'
    !
    !mumps_par%JOB = 6
    !
    !  PRINT*,' '
    !  PRINT*,' ---------------------------------------------------------------------------- '
    !  PRINT*,'              Solving linear system using MUMPS_5.5.1            '
    !  PRINT*,' ---------------------------------------------------------------------------- '
    !  PRINT*,' '
    !
    ! PRINT*,'ICNTL(2)',mumps_par%ICNTL(2)
    !
    !
    ! CALL DMUMPS(mumps_par)
    !
    !!PRINT*,'error??'
    !pause
    !
    ! Impongo la forma del pivoteo (porque no funciona con PORD [revisar instalacion])
    !
    PRINT*,'ICNTL(7)=',mumps_par%ICNTL(7)
    ! mumps_par%ICNTL(7)  = 6 ! Ordenamiento QAMD
    ! mumps_par%ICNTL(7) = 2 ! AMF
    !mumps_par%ICNTL(7)  = 5
    !
    PRINT*,'Actualizado ICNTL(7)=',mumps_par%ICNTL(7)
    !
    !
    PRINT*,'************ ANALISIS ***********'
    !
    mumps_par%JOB = 1 ! analisis
    !
    CALL DMUMPS(mumps_par)
    !
    !pause
    !
    PRINT*,'************ FACTORIZA ***********'
    !
    mumps_par%JOB = 2 ! factoriza
    !
    !PRINT*,'ICNTL(14)=',mumps_par%ICNTL(14)
    !mumps_par%ICNTL(14) = 50
    !PRINT*,'Actualizado ICNTL(14)=',mumps_par%ICNTL(14)
    !
    CALL DMUMPS(mumps_par)
    !
    !pause
    !
    PRINT*,'************ RESUELVE ***********'
    !
    !mumps_par%ICNTL(10) = 5 ! Refinamiento iterativo =0 por defecto, <0 no hace nada
    !
    mumps_par%JOB = 3 ! resuelve
    !
    !
    CALL DMUMPS(mumps_par)
    !
    !pause
    !  Solution has been assembled on the host
    !
    !IF ( mumps_par%MYID .EQ. 0 ) THEN
    !   WRITE( 6, * ) ' Solution is ',(mumps_par%RHS(I),I=1,mumps_par%N)
    !END IF
    !
    bb = 0.0_dp
    !
    !
    i=0
    !
    DO j=1,mumps_par%NRHS
       !
       i = i + 1
       !
       bb(:,i) = mumps_par%RHS((j-1)*mumps_par%N+1:j*mumps_par%N) 
       !
    END DO
    !
    !  Deallocate user data
    !
    IF ( mumps_par%MYID .EQ. 0 )THEN
       DEALLOCATE( mumps_par%IRN )
       DEALLOCATE( mumps_par%JCN )
       DEALLOCATE( mumps_par%A   )
       DEALLOCATE( mumps_par%RHS )
    END IF
    !
    !  Destroy the instance (deallocate internal data structures)
    !
    mumps_par%JOB = -2
    !
    CALL DMUMPS(mumps_par)
    !
    !CALL MPI_FINALIZE(IERR)
    !
    ! stop
    !
  END SUBROUTINE system_mumps_by_blocks
  !
  !
END MODULE mumps_system
