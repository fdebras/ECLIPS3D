MODULE mod_data

  IMPLICIT NONE

  !Parameters
  DOUBLE PRECISION,PUBLIC, PARAMETER :: pi=3.141592654d0 ! Well ... pi
  ! Resolution
  INTEGER, PUBLIC :: nlong
  INTEGER, PUBLIC :: nlat
  INTEGER, PUBLIC :: ntot
  INTEGER, PUBLIC :: nb
  ! data folder
  CHARACTER(LEN=100), PUBLIC :: DIRDATA

  !Planet inputs
  DOUBLE PRECISION, PUBLIC :: height ! Typical height considered
  DOUBLE PRECISION, PUBLIC :: phi  ! latitude of the beta plane
  DOUBLE PRECISION, PUBLIC :: g  !Gravity acceleration
  DOUBLE PRECISION, PUBLIC :: omega !Rotation rate
  DOUBLE PRECISION, PUBLIC :: gascons ! gas constant
  DOUBLE PRECISION, PUBLIC :: cp !specific heat capacity of air
  DOUBLE PRECISION, PUBLIC :: kappa ! gascons/cp
  DOUBLE PRECISION, PUBLIC :: rtot! planet radius
  DOUBLE PRECISION, PUBLIC :: p0 ! Reference pressure

  ! Drag and radiative timescales
  DOUBLE PRECISION, PUBLIC :: trad_invert
  DOUBLE PRECISION, PUBLIC :: tdrag_invert


  NAMELIST /resolution/ nlong,nlat,nb
  NAMELIST /folder/ DIRDATA
  NAMELIST /planet/ phi,height, &
  g, omega, gascons, cp, rtot, p0
  NAMELIST /timescales/ trad_invert,tdrag_invert

END MODULE
