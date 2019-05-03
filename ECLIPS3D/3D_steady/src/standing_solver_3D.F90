PROGRAM standing_solver_3D

  USE MPI
  USE mod_init_para
  USE mod_data
  USE mod_init_matrix
  USE mod_3D_fill_matrix
  USE mod_solver
  
  

  
  IMPLICIT NONE

  INTEGER ::  ntab, ndutab, ndvtab, ndwtab, ndptab, nqtab
  INTEGER :: num_args, n_input_file
  CHARACTER(LEN=300) :: c_namelist_file


  CALL init_para

  IF (myrow==0 .AND. mycol==0) THEN
    num_args=command_argument_count()
    IF (num_args .lt. 1) THEN
      WRITE(6,*) "Please supply a namelist file as a command line argument."
      STOP
    END IF
    CALL get_command_argument(1,c_namelist_file)
    WRITE(6,*)"parameter file name=",trim(c_namelist_file)
  END IF

  CALL mpi_bcast(c_namelist_file,300,mpi_character,0,MPI_COMM_WORLD,ierr)
  n_input_file = 65

  OPEN(UNIT=n_input_file,FILE=trim(c_namelist_file),STATUS="OLD")

  READ(n_input_file, resolution)
  READ(n_input_file,folder)
  READ(n_input_file,planet)
  READ(n_input_file,timescales)
  CLOSE(n_input_file)

  ntot=nlong*(2*nlat*nz+(nlat+1)*nz+2*nlat*(nz))
  kappa = gascons/cp
    
  ntab=6*nlong*(2*nlat*nz+(nlat+1)*nz+nlat*(nz+1))
  ndutab=nz*nlat*nlong*7
  ndvtab=nz*(nlat+1)*nlong*7
  ndptab=nz*nlat*nlong*8
  ndwtab=(nz+1)*nlat*nlong*10
  nqtab=2*nz*nlat*nlong
   
   
  IF (myrow==0 .AND. mycol==0) THEN 
    print *, 'Job starting'
  END IF 
  CALL init_matrix
  CALL fill_matrix(ntab, ndutab, ndvtab, ndptab, ndwtab, nqtab)
  CALL solve
  
  IF (myrow==0 .AND. mycol==0) THEN
    print *, 'All done'
  END IF


  CALL BLACS_GRIDEXIT( info_txt )
  CALL BLACS_EXIT(1)


  CALL MPI_FINALIZE(ierr)

END PROGRAM
  
