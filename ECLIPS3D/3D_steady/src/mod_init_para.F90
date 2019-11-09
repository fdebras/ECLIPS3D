MODULE mod_init_para

  USE MPI
  INTEGER, PUBLIC :: info_txt, nprocs, nprow, npcol
  INTEGER, PUBLIC :: ierr
  INTEGER, PUBLIC :: myrow, mycol

CONTAINS

  SUBROUTINE init_para


    IMPLICIT NONE


    EXTERNAL blacs_exit,blacs_get,blacs_gridexit,blacs_gridinfo,blacs_gridinit,&
    blacs_pinfo,blacs_setup

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)


    nprow = FLOOR(SQRT(REAL(nprocs)))
    npcol = nprocs/nprow
    DO WHILE (nprow*npcol .NE. nprocs)
      nprow = nprow - 1
      npcol = nprocs/nprow
    END DO

    CALL SL_INIT(info_txt,nprow,npcol)
    CALL BLACS_GRIDINFO(info_txt,nprow,npcol,myrow,mycol)


  END SUBROUTINE

      SUBROUTINE SL_INIT( ICTXT, NPROW, NPCOL )
! *
! *     .. Scalar Arguments ..
      INTEGER, INTENT(OUT) :: ICTXT
      INTEGER, INTENT(IN) :: NPCOL, NPROW
! *     ..

! *  Purpose
! *  =======
! *
! *  SL_INIT initializes an NPROW x NPCOL process grid using a row-major
! *  ordering  of  the  processes. This routine retrieves a default system
! *  context  which  will  include all available processes. In addition it
! *  spawns the processes if needed.
! *
! *  Arguments
! *  =========
! *
! *  ICTXT   (global output) INTEGER
! *          ICTXT specifies the BLACS context handle identifying the
! *          created process grid.  The context itself is global.
! *
! *  NPROW   (global input) INTEGER
! *          NPROW specifies the number of process rows in the grid
! *          to be created.
! *
! *  NPCOL   (global input) INTEGER
! *          NPCOL specifies the number of process columns in the grid
! *          to be created.
! *
! *  =====================================================================
! *
! *     .. Local Scalars ..
      INTEGER :: IAM, NPROCS
      EXTERNAL BLACS_PINFO, BLACS_SETUP, BLACS_GET, BLACS_GRIDINIT
! *     Get starting information
! *
      CALL BLACS_PINFO( IAM, NPROCS )
! *
! *     If machine needs additional set up, do it now
! *

      IF( NPROCS.LT.1 ) THEN
         IF( IAM.EQ.0 ) NPROCS = NPROW*NPCOL
         CALL BLACS_SETUP( IAM, NPROCS )
      END IF
! *
! *     Define process grid
! *
      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL )

      RETURN
    END SUBROUTINE SL_INIT


END MODULE
