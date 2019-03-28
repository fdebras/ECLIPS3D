MODULE mod_planet


  IMPLICIT NONE
  
  DOUBLE PRECISION, PARAMETER :: pi=3.141592654 ! Well ... pi
  
  LOGICAL, PARAMETER :: deep=.TRUE. ! true if deep atmosphere
  LOGICAL, PARAMETER :: g_var=.TRUE. ! True if g varies in atmosphere
  LOGICAL, PARAMETER :: eq_sym=.TRUE. ! True if symmetric about equator


  !DOUBLE PRECISION, PARAMETER :: tdrag=1.0E07
  !DOUBLE PRECISION, PARAMETER :: trad = 1.0E-5  
  ! Uncomment the planet you wish to use
  
  ! Earth
!  DOUBLE PRECISION, PUBLIC, PARAMETER :: ymax=88.9999 ! Maximum latitude
!  DOUBLE PRECISION, PUBLIC, PARAMETER :: height_min=501 ! Minimium height
!  DOUBLE PRECISION, PUBLIC, PARAMETER :: height_max=31499
!  DOUBLE PRECISION, PUBLIC, PARAMETER :: lambda_min=1.25001
!  DOUBLE PRECISION, PUBLIC, PARAMETER :: lambda_max=361.25001


 ! DOUBLE PRECISION,PUBLIC, PARAMETER :: g=9.80616   !Gravity acceleration 
 ! DOUBLE PRECISION,PUBLIC, PARAMETER :: omega=7.292E-5 !Rotation rate
  !REAL,PUBLIC, PARAMETER :: omega=2.0E-4 !Rotation rate
 ! DOUBLE PRECISION,PUBLIC, PARAMETER :: gascons=287.05 ! gas constant
 ! DOUBLE PRECISION,PUBLIC, PARAMETER :: cp=1005.0 !specific heat capacity of air
 ! DOUBLE PRECISION,PUBLIC, PARAMETER :: kappa= gascons/cp
 ! DOUBLE PRECISION,PUBLIC, PARAMETER :: rtot= 6371220.0 ! planet radius
 ! DOUBLE PRECISION,PUBLIC, PARAMETER :: depth=32000.0 ! Depth of atmosphere 
 ! DOUBLE PRECISION, PUBLIC, PARAMETER :: p0=1.0E5  

  ! HD209
  DOUBLE PRECISION, PUBLIC, PARAMETER :: ymax=88.9999 ! Maximum latitude
  DOUBLE PRECISION, PUBLIC, PARAMETER :: height_min=1.0E6 ! Minimium height
  DOUBLE PRECISION, PUBLIC, PARAMETER :: height_max=1.09E7
  DOUBLE PRECISION, PUBLIC, PARAMETER :: lambda_min=1.25001
  DOUBLE PRECISION, PUBLIC, PARAMETER :: lambda_max=361.25001


!  Newtonian
  DOUBLE PRECISION,PUBLIC, PARAMETER :: g=9.42  !Gravity acceleration 
  DOUBLE PRECISION,PUBLIC, PARAMETER :: omega=2.06E-5 !Rotation rat 
  DOUBLE PRECISION,PUBLIC, PARAMETER :: gascons=4593.0! gas constant
  DOUBLE PRECISION,PUBLIC, PARAMETER :: cp=14308.4 !specific heat capacity of air
  DOUBLE PRECISION,PUBLIC, PARAMETER :: kappa= gascons/cp
  DOUBLE PRECISION,PUBLIC, PARAMETER :: rtot= 9.44E7 ! planet radius
  DOUBLE PRECISION,PUBLIC, PARAMETER :: depth=1.1E7 ! Depth of atmosphere 
  DOUBLE PRECISION, PUBLIC, PARAMETER :: p0=2.2E7


END MODULE
