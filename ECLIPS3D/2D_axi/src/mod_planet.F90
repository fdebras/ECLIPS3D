MODULE mod_planet

! Uncomment the planet you wish to use
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: m=1 ! Longitidunal wave number
  LOGICAL, PUBLIC, PARAMETER :: deep=.TRUE. ! true if deep atmosphere
  LOGICAL, PUBLIC, PARAMETER :: g_var=.FALSE. ! True if g varies in atmosphere
  LOGICAL, PUBLIC, PARAMETER :: eq_sym=.TRUE. ! True if symmetric about equator
  DOUBLE PRECISION, PUBLIC, PARAMETER :: ymax=89.9999
  
  ! Earth
  DOUBLE PRECISION, PUBLIC, PARAMETER :: height_min=0.d0
  DOUBLE PRECISION, PUBLIC, PARAMETER :: height_max=80000.d0
  DOUBLE PRECISION,PUBLIC, PARAMETER :: g=9.80616   !Gravity acceleration 
  DOUBLE PRECISION,PUBLIC, PARAMETER :: omega=7.292E-5 !Rotation rate
!  DOUBLE PRECISION,PUBLIC, PARAMETER :: omega=2.0E-4 !Rotation rate
  DOUBLE PRECISION,PUBLIC, PARAMETER :: gascons=287.05 ! gas constant
  DOUBLE PRECISION,PUBLIC, PARAMETER :: cp=1005.0 !specific heat capacity of air
  DOUBLE PRECISION,PUBLIC, PARAMETER :: kappa= gascons/cp
  DOUBLE PRECISION,PUBLIC, PARAMETER :: rtot= 6371220.0 ! planet radius
  DOUBLE PRECISION,PUBLIC, PARAMETER :: depth=80000.0 ! Depth of atmosphere 
  DOUBLE PRECISION, PUBLIC, PARAMETER :: p0=1.0E5 

  ! HD209
!   DOUBLE PRECISION, PUBLIC, PARAMETER :: height_min=1.0E6
!   DOUBLE PRECISION, PUBLIC, PARAMETER :: height_max=1.09E7
!   DOUBLE PRECISION,PUBLIC, PARAMETER :: g=9.42  !Gravity acceleration 
!   DOUBLE PRECISION,PUBLIC, PARAMETER :: omega=2.06E-5 !Rotation rate
!   DOUBLE PRECISION,PUBLIC, PARAMETER :: gascons=4593.0! gas constant
!   DOUBLE PRECISION,PUBLIC, PARAMETER :: cp=14308.4 !specific heat capacity of air
!   DOUBLE PRECISION,PUBLIC, PARAMETER :: kappa= gascons/cp
!   DOUBLE PRECISION,PUBLIC, PARAMETER :: rtot= 9.44E7 ! planet radius
!   DOUBLE PRECISION,PUBLIC, PARAMETER :: depth=1.1E7 ! Depth of atmosphere 

END MODULE
