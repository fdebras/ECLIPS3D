PROGRAM read_solution
  IMPLICIT NONE

  INTEGER, PARAMETER :: nlong = 30
  INTEGER, PARAMETER :: nlat = 25
  INTEGER,PARAMETER :: nz = 33
  INTEGER,PARAMETER :: ntot = nlong*(2*nlat*nz+(nlat+1)*nz+2*nlat*(nz))
  DOUBLE PRECISION, DIMENSION(ntot) :: solution


  CHARACTER(LEN=*), PARAMETER :: DIRDATA='../data/'
  OPEN(unit=1,file=DIRDATA // 'solution.dat', &
    access='SEQUENTIAL')

  OPEN(unit=3,file=DIRDATA // 'solution_to_read.dat', &
    access="SEQUENTIAL")

  read(1,*)

  read(1,9999) solution

  write(3,888) DBLE(solution)

  close(1)
  close(3)

9999 FORMAT( D30.18)
888 FORMAT( 10E24.16 )

END PROGRAM
