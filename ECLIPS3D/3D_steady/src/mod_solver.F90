MODULE mod_solver

  USE mod_data
  USE mod_3D_fill_matrix
  USE mod_writing

  IMPLICIT NONE
  INTEGER, DIMENSION(50) :: descSolution
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: solution


CONTAINS

  SUBROUTINE solve

  IMPLICIT NONE
  INTEGER*8 :: LWORK

  INTEGER :: i,j,k,t, sel_dum, mput

  INTEGER :: IA,JA

  INTEGER :: info

  INTEGER :: ILO, IHI
  INTEGER, DIMENSION(:), ALLOCATABLE :: scale

  INTEGER, DIMENSION(:), ALLOCATABLE :: pivot ! for linear algebra
  INTEGER :: pivot_size


  DOUBLE PRECISION :: ONE=1.d0
  DOUBLE PRECISION :: ZERO=0.d0


  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: WORK

  INTEGER :: numroc
  DOUBLE PRECISION, DIMENSION(ntot) :: heating  ! heating, solution to obtain

  EXTERNAL  BLACS_GRIDINFO, DESCINIT, PDELSET,PDGETRF, PDGETRDS &
  , BLACS_GRIDEXIT, BLACS_EXIT

  print *, 'Matrix ready'

  LWORK=10*ntot
  ALLOCATE(WORK(LWORK))

  ! First, initialize the solution vector
  CALL DESCINIT(descSolution, ntot,1,nb,nb,0,0,info_txt,nlld,info)

  ALLOCATE(solution(nlld))

  open(unit=1,file=TRIM(DIRDATA) // "heating.dat", access= &
  'sequential')

  read(1,*) heating ! the initialisation is done with the heating array
  close(1)


  DO t=1, ntot
    CALL PDELSET(solution,t,1,descSolution,heating(t))
  END DO


  ILO=1
  IHI=ntot

  IA=1
  JA=1

  pivot_size= NUMROC(ntot,desc_mat(5),myrow,desc_mat(7),nprow) ! cf pzgetrf documentation
  pivot_size = pivot_size + desc_mat(5)

  ALLOCATE(pivot(pivot_size))


  print *, '     PLU decomposition ...'


  CALL PDGETRF(ntot, ntot, mat_evol, IA, JA, desc_mat, pivot, info )

  IF (info==0) THEN
  ELSE
    print *, info
    print *, "     PLU decomposition failure"
  END IF

  print *, '     Linear equation solver ...'
  CALL PDGETRS('N',ntot,1,mat_evol,IA,JA,desc_mat,pivot,solution,1,1,descSolution,info)

  IF (info==0) THEN

  ELSE
    print *, info
    print *, "     Linear solver failure"
  END IF


  IF (myrow==0 .and. mycol==0) THEN
    open(unit=44,file=TRIM(DIRDATA) // 'size.dat', &
    access="SEQUENTIAL")

    WRITE(44,*) nlat,nz,ntot
    CLOSE(44)
  END IF



  CALL PDLAWRITE(TRIM(DIRDATA) // 'solution.dat',ntot, &
       1,solution, 1,1,descSolution,0,0,WORK)

  DEALLOCATE(WORK)


9999 FORMAT( D30.18 )

  END SUBROUTINE

END MODULE
