PROGRAM normal_modes_2D_para

  USE MPI
  USE mod_init_para
  USE mod_data
  USE mod_init_matrix
  USE mod_fill_matrix
  USE mod_eigenvalues

  IMPLICIT NONE

  INTEGER ::  ntab, ndtab
  INTEGER :: num_args, n_input_file,n_timescale_file
  CHARACTER(LEN=300) :: c_namelist_file, c_timescale_file


  CALL init_para

  IF (myrow==0 .AND. mycol==0) THEN
    num_args=command_argument_count()
    IF (num_args .lt. 1) THEN
      WRITE(6,*) "Please supply a namelist file as a command line argument."
      STOP
    END IF
    CALL get_command_argument(1,c_namelist_file)
    CALL get_command_argument(2,c_timescale_file)
    WRITE(6,*)"parameter file name=",trim(c_namelist_file)
  END IF

  CALL mpi_bcast(c_namelist_file,300,mpi_character,0,MPI_COMM_WORLD,ierr)
  CALL mpi_bcast(c_timescale_file,300,mpi_character,0,MPI_COMM_WORLD,ierr)
  n_input_file = 65
  n_timescale_file = 75

  OPEN(UNIT=n_input_file,FILE=trim(c_namelist_file),STATUS="OLD")

  READ(n_input_file, resolution)
  READ(n_input_file,folder)
  READ(n_input_file,planet)

  CLOSE(n_input_file)

  OPEN(UNIT=n_timescale_file,FILE=trim(c_timescale_file),STATUS="OLD")
  READ(n_timescale_file,timescales)
  CLOSE(n_timescale_file)

  ntot=2*nlat*nz+(nlat+1)*nz+2*nlat*(nz-1)
  ntab = 6*(nlat*nz+(nlat+1)*nz+nlat*(nz+1))
  ndtab = 8*nlat*nz+7*nlat*(nz+1)+5*(nlat+1)*nz

  kappa = gascons/cp

  CALL init_matrix
  CALL fill_matrix(ntab,ndtab)
  CALL evalues



  CALL BLACS_GRIDEXIT( info_txt )
  CALL BLACS_EXIT(1)

!

  print *, 'lol'
  CALL MPI_FINALIZE(ierr)


END PROGRAM
