PROGRAM write_eigenvectors

  IMPLICIT NONE
  INTEGER :: nlong
  INTEGER :: nlat
  INTEGER :: nz
  INTEGER :: ntot

  COMPLEX*16 :: freq
  INTEGER :: i,j
  INTEGER :: num, num_prev

  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: eigenvec, tmp

  CHARACTER(LEN=*), PARAMETER :: DIRDATA ='../data/'

  OPEN(unit=4,file=DIRDATA // 'size.dat', &
    access="SEQUENTIAL")

  READ(4,*) nlong,nlat,nz,ntot
  CLOSE(4)


  ALLOCATE(eigenvec(ntot))
  ALLOCATE(tmp(ntot))


  OPEN(unit=1,file=DIRDATA // 'nummodes.dat', &
    access='SEQUENTIAL')
  OPEN(unit=2, file=DIRDATA // 'scalapack.dat', &
    access='SEQUENTIAL')
  OPEN(unit=3,file=DIRDATA // 'selected_modes.dat', &
    access="SEQUENTIAL")
   ! First line is the size, we dont care
  READ(2,*)
  READ(1,*) num, freq

  print *, num,freq

  num_prev=0
  ! We read the numbers of the selected eigenvectors
  DO WHILE (num>0) !num=0 is the end of the file
   IF (num==num_prev+1) THEN ! If the next eigenvectors is good,
    !we read and write it
      READ(2,9999) eigenvec

      WRITE(3,888) DBLE(freq), DIMAG(freq)
      WRITE(3,888) DBLE(eigenvec)
      WRITE(3,888) DIMAG(eigenvec)
    ELSE
    ! Else we read useless stuff until the good eigenvectors
      DO i=num_prev+1,num-1
        DO j=1,ntot
          READ(2,9999)
        END DO
      END DO

      READ(2,9999) eigenvec

      WRITE(3,888) DBLE(freq), DIMAG(freq)
      WRITE(3,888) DBLE(eigenvec)
      WRITE(3,888) DIMAG(eigenvec)
    END IF

    !Next iteration
    num_prev=num
    READ(1,*) num, freq
    print *, num
  END DO
  CLOSE(1)
  CLOSE(2)
  CLOSE(3)


  DEALLOCATE(eigenvec)
  DEALLOCATE(tmp)
9999 FORMAT( E15.8,E15.8 )
888 FORMAT( 10E14.6 )

END PROGRAM
