PROGRAM normal_modes_shallow
  
  USE MPI
  USE mod_init_para
  USE mod_data
  USE mod_init_matrix
  USE mod_fill_matrix
  USE mod_eigenvalues
  
  
  IMPLICIT NONE
  
  INTEGER :: num_args, n_input_file 
  CHARACTER(LEN=300) :: c_namelist_file

  INTEGER :: nu_tab, nv_tab, nh_tab
  
  
  
  CALL init_para
  
  n_input_file = 65
 
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

  OPEN(UNIT=n_input_file,FILE=trim(c_namelist_file),STATUS="OLD")

  READ(n_input_file, resolution)
  READ(n_input_file,folder)
  READ(n_input_file,planet)
  READ(n_input_file,timescales)
  CLOSE(n_input_file)

  ntot=2*nlong*nlat + nlong*(nlat+1) ! u and h same size, v from equator to pole
  kappa = gascons/cp
  nu_tab=4*nlong*nlat
  nv_tab=4*nlong*(nlat+1)
  nh_tab = 7*nlong*nlat
  
  
  CALL init_matrix
  CALL fill_matrix(nu_tab,nv_tab,nh_tab)
  CALL evalues
  
  ! Wait call for all processors, I broadcast some useless stuff
  CALL mpi_bcast(c_namelist_file,300,mpi_character,0,MPI_COMM_WORLD,ierr)

  CALL BLACS_GRIDEXIT( info_txt )
  CALL BLACS_EXIT(1)

!   

  print *, 'all done'
  CALL MPI_FINALIZE(ierr)

  
END PROGRAM
