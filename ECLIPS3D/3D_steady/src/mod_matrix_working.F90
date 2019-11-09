! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE mod_3D_fill_matrix

  USE mod_init_para
  USE mod_data
  USE mod_init_matrix

  IMPLICIT NONE
  ! Description:
  !  A routine to calculate the normal modes from linear perturbations
  !  around a steady state where U,P and T depend on R and PHI only

  ! Method:
  !  Takes a steady state from the UM, then
  !  solve an eigenvalue problem to find frequencies and modes. Special
  !  attention must be taken on the grid for calculating the values of
  !  the perturbed variables.
  CONTAINS

  SUBROUTINE fill_matrix(ntab,ndutab,ndvtab,ndptab,ndwtab,nqtab)

    IMPLICIT NONE
    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Data arrays
    !----------------------------------------------------------
    !----------------------------------------------------------
    INTEGER, INTENT(IN) ::  ntab, ndutab, ndvtab, ndwtab, ndptab, &
    nqtab

    DOUBLE PRECISION, ALLOCATABLE :: dataarray(:),duarray(:), &
    dvarray(:), dwarray(:), dparray(:), qarray(:)

    !----------------------------------------------------------
    !----------------------------------------------------------
    !Local Variables
    !----------------------------------------------------------
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! Steady state variables (see grid)
    !----------------------------------------------------------

    DOUBLE PRECISION, DIMENSION(nlong,nlat,nz) :: Us_u, Vs_u,Ps_u, Ws_u,Thetas_u, &
    Rho_u ,C_u_square
            ! Steady state on u levels

    DOUBLE PRECISION, DIMENSION(nlong,0:nlat,nz) :: Us_v, Vs_v,Ps_v,Ws_v,Thetas_v, &
    Rho_v , C_v_square
    ! Steady state on v levels

    DOUBLE PRECISION, DIMENSION(nlong,nlat,nz) :: Us_p, Vs_p,Ps_p, Ws_p, Thetas_p, &
    Rho_p, C_p_square, N_p_square, Q_p
    ! Steady state on p levels

    DOUBLE PRECISION, DIMENSION(nlong,nlat,0:nz) :: Us_w, Vs_w,Ps_w, Ws_w,Thetas_w, &
    Rho_w ,C_w_square, N_w_square
    ! Steady state on w levels

    DOUBLE PRECISION :: Q_theta(nlong,nlat,nz-1)
    ! exception for Q_theta

    DOUBLE PRECISION, DIMENSION(nlong,nlat,nz) :: dUs_dl_u, dUs_dphi_u, &
    dUs_dr_u, dPs_dl_u, dRho_dl_u,dRho_dphi_u,dRho_dr_u
    ! Derivatives on u levels

    DOUBLE PRECISION, DIMENSION(nlong,0:nlat,nz) :: dVs_dl_v, dVs_dphi_v, &
    dVs_dr_v, dPs_dphi_v, dRho_dl_v, dRho_dphi_v, dRho_dr_v
    ! Derivatives on v levels

    DOUBLE PRECISION, DIMENSION(nlong,nlat,nz) :: dUs_dl_p, dVs_dphi_p, &
    dWs_dr_p, dPs_dl_p, dPs_dphi_p, dRho_dl_p, dRho_dphi_p, dThetas_dr_p
    ! derivaties on p levels

    DOUBLE PRECISION, DIMENSION(nlong,nlat,0:nz) :: dWs_dl_w, dWs_dphi_w, &
    dWs_dr_W, dPs_dr_w, dThetas_dl_w,dThetas_dphi_w, dThetas_dr_w, &
    dRho_dl_w, dRho_dphi_w, dRho_dr_w
     ! Derivatives on w levels


    DOUBLE PRECISION :: trad_theta(nlong,nlat,0:nz),trad_p(nlong,nlat,nz), &
    tdrag_u(nlong,nlat,nz), tdrag_v(nlong,0:nlat,nz),tdrag_w(nlong,nlat,0:nz)

!
!     !----------------------------------------------------------
!     ! Other variables
!     !----------------------------------------------------------

    INTEGER :: i,j,k,t,eig ! Loops integers

    DOUBLE PRECISION :: lambda_u(nlong), lambda_v(nlong), phi_u(nlat), &
    phi_v(0:nlat) ,z_u(nz), z_w(0:nz), gz_u(nz), gz_w(0:nz)
    ! Longitude, Latitude, Height,g on u, v and w levels
    DOUBLE PRECISION :: dlambda,dphi,dz ! longitude, latitude and height steps
    DOUBLE PRECISION :: S(ntot) ! Just used for computing the matrix

    INTEGER :: info
    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Subroutine
    !----------------------------------------------------------
    !----------------------------------------------------------

    open(unit=1,file=TRIM(DIRDATA) // "data.dat", access= &
    'sequential')

    open(unit=2,file=TRIM(DIRDATA) // "dudata.dat", access= &
    'sequential')

    open(unit=3,file=TRIM(DIRDATA) // "dvdata.dat", access= &
    'sequential')

    open(unit=4,file=TRIM(DIRDATA) // "dpdata.dat", access= &
    'sequential')

    open(unit=10,file=TRIM(DIRDATA) // "dwdata.dat", access= &
    'sequential')

    open(unit=11,file=TRIM(DIRDATA) // "qdata.dat", access= &
    'sequential')



    ALLOCATE(dataarray(ntab))
    ALLOCATE(duarray(ndutab))
    ALLOCATE(dvarray(ndvtab))
    ALLOCATE(dparray(ndptab))
    ALLOCATE(dwarray(ndwtab))
    ALLOCATE(qarray(nqtab))


    read(1,*) dataarray
    read(2,*) duarray
    read(3,*) dvarray
    read(4,*) dparray
    read(10,*) dwarray
    read(11,*) qarray


    close(1)
    close(2)
    close(3)
    close(4)
    close(10)
    close(11)

    !----------------------------------------------------------
    ! Initialisation of the variables
    !------------------------------------------------ß----------

    !-----------------------------------------------------------
    ! Reading the input array according to grid
    t=1
    !------------
    !U
    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Us_u(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          Us_v(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Us_p(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          Us_w(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO


    !------------
    !V
    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Vs_u(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          Vs_v(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Vs_p(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          Vs_w(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO


    !------------
    !P
    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Ps_u(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          Ps_v(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Ps_p(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          Ps_w(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO


    !------------
    !W
    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Ws_u(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          Ws_v(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Ws_p(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          Ws_w(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO



    !------------
    !Theta
    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Thetas_u(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          Thetas_v(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Thetas_p(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          Thetas_w(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO


    !------------
    !Rho
    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Rho_u(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          Rho_v(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Rho_p(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          Rho_w(i,j,k)=dataarray(t)
          t=t+1
        END DO
      END DO
    END DO



    !------------------------------------------------------------------------
    ! Reading the u derivative array
    t=1
    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dUs_dl_u(i,j,k)=duarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dUs_dphi_u(i,j,k)=duarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dUs_dr_u(i,j,k)=duarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dPs_dl_u(i,j,k)=duarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dRho_dl_u(i,j,k)=duarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dRho_dphi_u(i,j,k)=duarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dRho_dr_u(i,j,k)=duarray(t)
          t=t+1
        END DO
      END DO
    END DO

    !------------------------------------------------------------------------
    ! Reading the v derivative array

    t=1
    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          dVs_dl_v(i,j,k)=dvarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          dVs_dphi_v(i,j,k)=dvarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          dVs_dr_v(i,j,k)=dvarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          dPs_dphi_v(i,j,k)=dvarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          dRho_dl_v(i,j,k)=dvarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          dRho_dphi_v(i,j,k)=dvarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          dRho_dr_v(i,j,k)=dvarray(t)
          t=t+1
        END DO
      END DO
    END DO


    !------------------------------------------------------------------------
    ! Reading the p derivative array

    t=1
    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dUs_dl_p(i,j,k)=dparray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dVs_dphi_p(i,j,k)=dparray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dPs_dl_p(i,j,k)=dparray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dPs_dphi_p(i,j,k)=dparray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dWs_dr_p(i,j,k)=dparray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dThetas_dr_p(i,j,k)=dparray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dRho_dl_p(i,j,k)=dparray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          dRho_dphi_p(i,j,k)=dparray(t)
          t=t+1
        END DO
      END DO
    END DO

    !------------------------------------------------------------------------
    ! Reading the w derivative array


    t=1

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dPs_dr_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dWs_dl_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dWs_dphi_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dWs_dr_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dThetas_dl_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dThetas_dphi_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dThetas_dr_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO


    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dRho_dl_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dRho_dphi_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          dRho_dr_w(i,j,k)=dwarray(t)
          t=t+1
        END DO
      END DO
    END DO

   print *, t
    !------------------------------------------------------------------------
    ! Reading the heating rate array
    t=1

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          Q_p(i,j,k)=qarray(t)
          t=t+1
        END DO
      END DO
    END DO


    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz-1
          Q_theta(i,j,k)=qarray(t)
          t=t+1
        END DO
      END DO
    END DO



    ! Creating the radiative timescale array
    CALL set_trad(trad_theta, trad_p, Ps_w, Ps_p)
    CALL set_tdrag(tdrag_u, tdrag_v, tdrag_w, Ps_u, Ps_v, Ps_w)

    DO i=1,nlong
      DO j=1,nlat
        DO k=0,nz
          tdrag_w(i,j,k)=0.0
        END DO
      END DO
   END DO


    IF (myrow==0 .and. mycol==0) THEN
      print *, 'trad set'
    END IF
   !print*, Theta_s

    !Defining height and latitude steps

    dlambda=2*pi/nlong
    dphi=ymax/180.*pi/(nlat)
    dz=(height_max-height_min)/nz


    DO j=1,nlat
      phi_u(j)=(j-0.5)*dphi
      phi_v(j)=(j)*dphi
    END DO
    phi_v(0)=0.0


    ! Defining height and gravity
    DO k=0,nz-1
      IF (deep) THEN
        z_w(k)=rtot+height_min+k*dz
        z_u(k+1)=rtot+height_min+(k+0.5)*dz
      ELSE
        z_w(k)=rtot
        z_u(k+1)=rtot
      END IF

      IF (g_var) THEN
      	gz_w(k)=g*rtot*rtot/(z_w(k)*z_w(k))
      	gz_u(k+1)=g*rtot*rtot/(z_u(k+1)*z_u(k+1))
      ELSE
        gz_w(k)=g
        gz_u(k+1)=g
      END IF
    END DO

    IF (deep) THEN
      z_w(nz)=rtot+height_min+nz*dz
      ELSE
      z_w(nz)=rtot
    END IF

    IF (g_var) THEN
      gz_w(nz)=g*rtot*rtot/(z_w(nz)*z_w(nz))
    ELSE
      gz_w(nz)=g
    END IF

   !Defining latitude full and half LATERAL levels
   if (myrow==0 .and. mycol==0) THEN
     print *, z_u
     print*, z_w
   END IF

	DO i=1,nlong
	  DO j=1,nlat
	    DO k=1,nz

          N_w_square(i,j,k)=gz_w(k)*(dThetas_dr_w(i,j,k))/(Thetas_w(i,j,k))

          N_p_square(i,j,k)=gz_u(k)*(dThetas_dr_p(i,j,k))/(Thetas_p(i,j,k))

          C_u_square(i,j,k)=gascons/(1.0-kappa)*(Ps_u(i,j,k)/p0)**(kappa)* &
          Thetas_u(i,j,k)

          C_v_square(i,j,k)=gascons/(1.0-kappa)*(Ps_v(i,j,k)/p0)**(kappa)* &
          Thetas_v(i,j,k)

          C_w_square(i,j,k)=gascons/(1.0-kappa)*(Ps_w(i,j,k)/p0)**(kappa)* &
          Thetas_w(i,j,k)

          C_p_square(i,j,k)=gascons/(1.0-kappa)*(Ps_p(i,j,k)/p0)**(kappa)* &
          Thetas_p(i,j,k)
        END DO
      END DO
    END DO

    j=0
    DO i=1,nlong
      DO k=1,nz
        C_v_square(i,j,k)=gascons/(1.0-kappa)*(Ps_v(i,j,k)/p0)**(kappa)* &
        Thetas_v(i,j,k)
      END DO
    END DO

    k=0
    DO i=1,nlong
      DO j=1,nlat
        N_w_square(i,j,k)=gz_w(k)*(dThetas_dr_w(i,j,k))/(Thetas_w(i,j,k))

        C_w_square(i,j,k)=gascons/(1.0-kappa)*(Ps_w(i,j,k)/p0)**(kappa)* &
        Thetas_w(i,j,k)
      END DO
    END DO
    ! Sound Speed and Brunt Vaisala frequency with c=gamma*R*T

    !----------------------------------------------------------
    ! Seek for eigenvectors
    !----------------------------------------------------------

    ! To create the matrix, we will just set all the variables to
    ! 0 at each point of the grid and allow one to be non zero at one point
    ! of the grid. Then we will calculate the evolution coefficients from the
    ! perturbed equations due to this value and store it in the matrix from
    ! which we will get the eigenvectors and values.



    DO i=1,ntot
      S=0.0
      S(i)=1.0

      CALL evol_coeff(i,dlambda,dphi, &
      dz,S,lambda_u,lambda_v,phi_u,phi_v,z_u,z_w,gz_u,gz_w, &
      N_w_square, N_p_square, C_u_square,C_v_square, C_w_square, C_p_square, &
      Us_u, dUs_dl_u, dUs_dphi_u, dUs_dr_u, Vs_u,dPs_dl_u, Ws_u, Rho_u, dRho_dl_u, &
      dRho_dphi_u, dRho_dr_u, Us_v,Vs_v, dVs_dl_v, dVs_dphi_v, dVs_dr_v, dPs_dphi_v, &
      Ws_v, Rho_v, dRho_dl_v, dRho_dphi_v, dRho_dr_v, &
      Us_w, Vs_w,Ps_w,dPs_dr_w, Ws_w, dWs_dl_w, dWs_dphi_w, dWs_dr_w,Thetas_w,dThetas_dl_w, &
      dThetas_dphi_w, Rho_w, dRho_dl_w, dRho_dphi_w, dRho_dr_w, Us_p, dUs_dl_p, Vs_p, &
      dVs_dphi_p,Ps_p,dPs_dl_p, dPs_dphi_p, Ws_p, dWs_dr_p, Rho_p, dRho_dl_p, &
      dRho_dphi_p,Q_p,Q_theta, tdrag_u, tdrag_v,trad_p,tdrag_w, trad_theta )


    END DO

    IF (myrow==0 .AND. mycol==0) THEN
      open(unit=44,file=TRIM(DIRDATA) // 'initial_state.dat', &
      access='SEQUENTIAL')

       WRITE(44,*) nlong,nlat,nz,height_max
!
       WRITE(44,888) Us_u
       WRITE(44,888) Vs_v
       WRITE(44,888) Ws_w
       WRITE(44,888) Ps_p
       WRITE(44,888) Thetas_w


       WRITE(44,888) Rho_u
       WRITE(44,888) Rho_v
       WRITE(44,888) Rho_w
       WRITE(44,888) Rho_p
       WRITE(44,888) C_p_square
       WRITE(44,888) N_w_square
      CLOSE(44)
    END IF


  888 FORMAT( 10E24.16 )

  END SUBROUTINE

  SUBROUTINE evol_coeff(pos,dlambda,dphi, &
      dz,S,lambda_u,lambda_v,phi_u,phi_v,z_u,z_w,gz_u,gz_w, &
      N_w_square, N_p_square, C_u_square,C_v_square, C_w_square, C_p_square, &
      Us_u, dUs_dl_u, dUs_dphi_u, dUs_dr_u, Vs_u,dPs_dl_u, Ws_u, Rho_u, dRho_dl_u, &
      dRho_dphi_u, dRho_dr_u, Us_v,Vs_v, dVs_dl_v, dVs_dphi_v, dVs_dr_v, dPs_dphi_v, &
      Ws_v, Rho_v, dRho_dl_v, dRho_dphi_v, dRho_dr_v, &
      Us_w, Vs_w,Ps_w, dPs_dr_w, Ws_w, dWs_dl_w, dWs_dphi_w, dWs_dr_w,Thetas_w,dThetas_dl_w, &
      dThetas_dphi_w, Rho_w, dRho_dl_w, dRho_dphi_w, dRho_dr_w, Us_p, dUs_dl_p, Vs_p, &
      dVs_dphi_p,Ps_p,dPs_dl_p, dPs_dphi_p, Ws_p, dWs_dr_p, Rho_p, dRho_dl_p, &
      dRho_dphi_p, Q_p,Q_theta,tdrag_u, tdrag_v, trad_p,tdrag_w, trad_theta )


    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Function to calculate the evolution of perturbed variables
    !----------------------------------------------------------
    !----------------------------------------------------------


    IMPLICIT NONE

    ! Declaring the inputs already known by sedulous reader
    INTEGER, INTENT(IN) :: pos
    DOUBLE PRECISION, INTENT(IN) :: dlambda,dphi,dz
    DOUBLE PRECISION, INTENT(IN) :: S(ntot)

    DOUBLE PRECISION, INTENT(IN) :: lambda_u(nlong), lambda_v(nlong), &
    phi_u(nlat),phi_v(0:nlat), z_u(nz),z_w(0:nz), gz_u(nz), gz_w(0:nz)

    ! Brunt vaisala and sound speed
    DOUBLE PRECISION, INTENT(IN) :: C_u_square(nlong,nlat,nz), &
    C_v_square(nlong,0:nlat,nz),C_w_square(nlong,nlat,0:nz), &
    C_p_square(nlong,nlat,nz), N_w_square(nlong,nlat,0:nz), &
    N_p_square(nlong,nlat,nz)

    !Steady on u
    DOUBLE PRECISION, DIMENSION(nlong,nlat,nz), INTENT(IN) :: Us_u, dUs_dl_u, dUs_dphi_u,&
    dUs_dr_u, Vs_u,dPs_dl_u, Ws_u, Rho_u, dRho_dl_u, dRho_dphi_u, dRho_dr_u


    !Steady on v
    DOUBLE PRECISION, DIMENSION(nlong,0:nlat,nz), INTENT(IN) :: Us_v, Vs_v, dVs_dl_v, &
    dVs_dphi_v, dVs_dr_v, dPs_dphi_v, Ws_v, Rho_v, dRho_dl_v, dRho_dphi_v, dRho_dr_v


    !Steady on w
    DOUBLE PRECISION, DIMENSION(nlong,nlat,0:nz), INTENT(IN) :: Us_w, Vs_w,Ps_w,dPs_dr_w, &
    Ws_w, dWs_dl_w, dWs_dphi_w, dWs_dr_w, Thetas_w, dThetas_dl_w, dThetas_dphi_w, Rho_w,&
    dRho_dl_w,dRho_dphi_w, dRho_dr_w

    !Steady on p
    DOUBLE PRECISION, DIMENSION(nlong,nlat,nz), INTENT(IN) :: Us_p, dUs_dl_p, Vs_p, &
    dVs_dphi_p,Ps_p, dPs_dl_p, dPs_dphi_p, Ws_p, dWs_dr_p,Rho_p, dRho_dl_p, dRho_dphi_p, &
    Q_p

    DOUBLE PRECISION, INTENT(IN) :: Q_theta(nlong,nlat,nz-1)

    DOUBLE PRECISION, INTENT(IN) :: trad_theta(nlong,nlat,0:nz), trad_p(nlong,nlat,nz), &
    tdrag_u(nlong,nlat,nz), tdrag_v(nlong,0:nlat,nz), tdrag_w(nlong,nlat,0:nz)

    DOUBLE PRECISION :: u(nlong,nlat,nz), v(nlong,0:nlat,nz), &
    p(nlong,nlat,nz), w(nlong,nlat,0:nz),theta(nlong,nlat,0:nz)
    ! perturbed quantities, all zero but one

    DOUBLE PRECISION :: alpha ! to fill the matrix


!----------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
! diffusion attempt
    DOUBLE PRECISION :: diff_u(nlong,nlat,nz), diff_v(nlong,0:nlat,nz)



    INTEGER :: i,j,k,t ! just for convenience and loops
  	!----------------------------------------------------------
    ! Initialisation of the non zero variable
    !----------------------------------------------------------

  	! Just in case
  	u=0.0
  	v=0.0
  	p=0.0
  	w=0.0
  	theta=0.0

  	! Reading the stored array


  	t=1
  	DO i=1,nlong
  	  DO j=1,nlat
  	    DO k=1,nz
  	      u(i,j,k)=S(t)
  	      t=t+1
  	    END DO
  	  END DO
    END DO

  	DO i=1,nlong
  	  DO j=0,nlat
  	    DO k=1,nz
  	      v(i,j,k)=S(t)
  	      t=t+1
  	    END DO
  	  END DO
    END DO

  	DO i=1,nlong
  	  DO j=1,nlat
  	    DO k=1,nz
  	      p(i,j,k)=S(t)
  	      t=t+1
  	    END DO
  	  END DO
    END DO

 	DO i=1,nlong
  	  DO j=1,nlat
  	    DO k=1,nz-1
  	      w(i,j,k)=S(t)
  	      t=t+1
  	    END DO
  	  END DO
    END DO

  	DO i=1,nlong
  	  DO j=1,nlat
  	    DO k=1,nz-1
  	      theta(i,j,k)=S(t)
  	      t=t+1
  	    END DO
  	  END DO
    END DO


  	!Boundary conditions : always zero

    !----------------------------------------------------------
    ! Diffusion
    !----------------------------------------------------------

    !diff_u=diffusion_u(u,ddu_ddl_on_u(u,dlambda),du_dphi_u, &
    !ddu_ddphi_on_u(u,dphi),du_dr_u,ddu_ddr_on_u(u,dz),&
    !Rho_u,dRho_dr_u,ddu_ddr_on_u(Rho_u,dz),phi_u,z_u)

    !diff_v = diffusion_v(v,ddv_ddl_on_v(v,dlambda),dv_dphi_v, &
    !ddv_ddphi_on_v(v,dphi),dv_dr_v,ddv_ddr_on_v(v,dz), &
    !Rho_v, dRho_dr_v, ddv_ddr_on_v(Rho_v,dz),phi_v,z_u)



  !----------------------------------------------------------
    ! FILL THE MATRIX
    !----------------------------------------------------------

  	t=1
  	! --------------------- u evolution ------------------------------------
  	DO i=1, nlong
  	  DO j=1,nlat
  	    DO k=1,nz
  	      alpha=u(i,j,k)*((dUs_dl_u(i,j,k)-Us_u(i,j,k)*dRho_dl_u(i,j,k)/Rho_u(i,j,k)) &
  	      /(z_u(k)*COS(phi_u(j))) - Vs_u(i,j,k)*dRho_dphi_u(i,j,k)/(z_u(k)*Rho_u(i,j,k)) - &
  	      Ws_u(i,j,k)*dRho_dr_u(i,j,k)/Rho_u(i,j,k) + Ws_u(i,j,k)/z_u(k) - &
  	      TAN(phi_u(j))*Vs_u(i,j,k)/z_u(k)) + &
  	      du_dl_on_u(u,i,j,k,dlambda)*Us_u(i,j,k)/(z_u(k)*COS(phi_u(j))) + &
  	      du_dphi_on_u(u,i,j,k,dphi)*Vs_u(i,j,k)/z_u(k) + &
  	      du_dr_on_u(u,i,j,k,dz)*Ws_u(i,j,k) + &
  	      v_on_u(v,i,j,k)*(dUs_dphi_u(i,j,k)/z_u(k)-2.0*omega*SIN(phi_u(j))- &
  	      Us_u(i,j,k)*TAN(phi_u(j))/z_u(k) ) + &
  	      w_on_u(w,i,j,k)*(dUs_dr_u(i,j,k)+2.0*omega*COS(phi_u(j))+Us_u(i,j,k)/z_u(k)) + &
  	      p_on_u(p,i,j,k)*(-dPs_dl_u(i,j,k)/(C_u_square(i,j,k)*Rho_u(i,j,k)*z_u(k)*COS(phi_u(j)))) + &
  	      dp_dl_on_u(p,i,j,k,dlambda)*(1.0/(z_u(k)*COS(phi_u(j)))) + &
  	      w_on_u(theta,i,j,k)*(dPs_dl_u(i,j,k)/(gz_u(k)*z_u(k)*Rho_u(i,j,k)*COS(phi_u(j)))) + &
              !(-1)*diff_u(i,j,k)*tdrag_u(i,j,k)
  	      u(i,j,k)*tdrag_u(i,j,k)
  	      
              CALL PDELSET(mat_evol,t,pos,desc_mat,alpha)
  	      t=t+1
  	    END DO
  	  END DO
  	END DO


  	! --------------------- v evolution ------------------------------------
  	DO i=1, nlong
  	  DO j=0,nlat
  	    DO k=1,nz
  	      IF (j==0) THEN
  	        alpha=v(i,j,k)*tdrag
  	      ELSE IF (j==nlat) THEN
  	        alpha=v(i,j,k)*tdrag
  	      ELSE
  	        alpha=u_on_v(u,i,j,k)*(dVs_dl_v(i,j,k)/(z_u(k)*COS(phi_v(j)))+ &
  	        2.0*omega*SIN(phi_v(j)) + &
  	        2.0*TAN(phi_v(j))*Us_v(i,j,k)/z_u(k)) + &
  	        v(i,j,k)*(-Us_v(i,j,k)*dRho_dl_v(i,j,k)/(z_u(k)*Rho_v(i,j,k)*COS(phi_v(j))) + &
  	        dVs_dphi_v(i,j,k)/z_u(k) - Vs_v(i,j,k)*dRho_dphi_v(i,j,k)/ &
  	        (z_u(k)*Rho_v(i,j,k)) - Ws_v(i,j,k)*dRho_dr_v(i,j,k)/Rho_v(i,j,k) + &
  	        Ws_v(i,j,k)/z_u(k)) + &
  	        dv_dl_on_v(v,i,j,k,dlambda)*(Us_v(i,j,k)/(z_u(k)*COS(phi_v(j)))) + &
  	        dv_dphi_on_v(v,i,j,k,dphi)*(Vs_v(i,j,k)/z_u(k)) + &
  	        dv_dr_on_v(v,i,j,k,dz)*Ws_v(i,j,k) + &
  	        w_on_v(w,i,j,k)*(dVs_dr_v(i,j,k)+Vs_v(i,j,k)/z_u(k)) + &
  	        p_on_v(p,i,j,k)*(-dPs_dphi_v(i,j,k)/(C_v_square(i,j,k)*Rho_v(i,j,k)*z_u(k))) + &
  	        dp_dphi_on_v(p,i,j,k,dphi)*(1.0/z_u(k)) + &
  	        w_on_v(theta,i,j,k)*(dPs_dphi_v(i,j,k)/(gz_u(k)*Rho_v(i,j,k)*z_u(k))) + &
  	        !(-1)*diff_v(i,j,k)*tdrag_v(i,j,k)
  	        v(i,j,k)*tdrag_v(i,j,k)
  	      END IF

  	      CALL PDELSET(mat_evol,t,pos,desc_mat,alpha)
  	      t=t+1
  	    END DO
  	  END DO
  	END DO

  	! --------------------- p evolution ------------------------------------
  	DO i=1, nlong
  	  DO j=1,nlat
  	    DO k=1,nz
  	      alpha =u_on_p(u,i,j,k) * ((dPs_dl_p(i,j,k) - C_p_square(i,j,k)*dRho_dl_p(i,j,k))/ &
          (Rho_p(i,j,k)*z_u(k)*COS(phi_u(j))) ) + &
          du_dl_on_p(u,i,j,k,dlambda)*C_p_square(i,j,k)/(z_u(k)*COS(phi_u(j))) + &
          v_on_p(v,i,j,k) * ((dPs_dphi_p(i,j,k) - C_p_square(i,j,k)*dRho_dphi_p(i,j,k))/ &
          (Rho_p(i,j,k)*z_u(k)) - C_p_square(i,j,k)*TAN(phi_u(j))/z_u(k) ) + &
          dv_dphi_on_p(v,i,j,k,dphi)*C_p_square(i,j,k)/z_u(k) + &
          w_on_p(w,i,j,k)*(C_p_square(i,j,k) * (2.0/z_u(k) + N_p_square(i,j,k)/gz_u(k)) ) + &
          dw_dr_on_p(w,i,j,k,dz)*C_p_square(i,j,k) + &
          p(i,j,k)* (1.0+kappa/(1.0-kappa)) * (dUs_dl_p(i,j,k)/(z_u(k)*COS(phi_u(j))) + &
          dVs_dphi_p(i,j,k)/z_u(k) - Vs_p(i,j,k)*TAN(phi_u(j))/z_u(k) + &
          dWs_dr_p(i,j,k) + 2.d0*Ws_p(i,j,k)/z_u(k) ) + &
          dp_dl_on_p(p,i,j,k,dlambda)*(Us_p(i,j,k)/ (z_u(k)*COS(phi_u(j))) ) + &
          dp_dphi_on_p(p,i,j,k,dphi)*Vs_p(i,j,k)/z_u(k) + &
          dp_dr_on_p(p,i,j,k,dz)*Ws_p(i,j,k) + & ! Additional terms
          !Q_p(i,j,k) / Rho_p(i,j,k)* ( &
          !theta_p(i,j,k) / (gz_u(k)) - p(i,j,k)/C_p_square(i,j,k) ) +&
          w_on_p(theta,i,j,k) * trad_p(i,j,k) * C_p_square(i,j,k)/gz_u(k) + &
          p(i,j,k) * trad_p(i,j,k) * C_p_square(i,j,k) * Rho_p(i,j,k) * &
          (kappa) / Ps_p(i,j,k)

              CALL PDELSET(mat_evol,t,pos,desc_mat,alpha)
  	      t=t+1
  	    END DO
  	  END DO
  	END DO

  	! --------------------- w evolution ------------------------------------
  	DO i=1, nlong
  	  DO j=1,nlat
  	    DO k=1,nz-1
  	      alpha = u_on_w(u,i,j,k)*(-2.0*omega*COS(phi_u(j))-2.0*Us_w(i,j,k)/z_w(k) + &
  	          dWs_dl_w(i,j,k)/(z_w(k)*COS(phi_u(j))) ) + &
              v_on_w(v,i,j,k)*(-2.0*Vs_w(i,j,k)/z_w(k) + dWs_dphi_w(i,j,k)/z_w(k) ) + &
              w(i,j,k)*(-Us_w(i,j,k)*dRho_dl_w(i,j,k)/(Rho_w(i,j,k)*z_w(k)*COS(phi_u(j))) - &
              Vs_w(i,j,k)*dRho_dphi_w(i,j,k)/(z_w(k)*Rho_w(i,j,k)) + &
              dWs_dr_w(i,j,k) - Ws_w(i,j,k)*dRho_dr_w(i,j,k)/Rho_w(i,j,k) ) + &
              dw_dl_on_w(w,i,j,k,dlambda)*Us_w(i,j,k)/(z_w(k)*COS(phi_u(j))) + &
              dw_dphi_on_w(w,i,j,k,dphi)*Vs_w(i,j,k)/z_w(k) + &
              dw_dr_on_w(w,i,j,k,dz)*Ws_w(i,j,k) + &
              p_on_w(p,i,j,k)*(-dPs_dr_w(i,j,k)/(C_w_square(i,j,k)*Rho_w(i,j,k))) + &
              dp_dr_on_w(p,i,j,k,dz) + &
              theta(i,j,k)*(dPs_dr_w(i,j,k)/(Rho_w(i,j,k)*gz_w(k))) + &

              w(i,j,k)*tdrag_w(i,j,k)

  	      CALL PDELSET(mat_evol,t,pos,desc_mat,alpha)
  	      t=t+1
  	    END DO
  	  END DO
  	END DO


  	! --------------------- theta evolution ------------------------------------
  	DO i=1, nlong
  	  DO j=1,nlat
  	    DO k=1,nz-1
  	      alpha = u_on_w(u,i,j,k)*(gz_w(k)*dThetas_dl_w(i,j,k)/(z_w(k)*COS(phi_u(j))* &
              Thetas_w(i,j,k))) + &
              v_on_w(v,i,j,k)*(gz_w(k)*dThetas_dphi_w(i,j,k)/(z_w(k)*Thetas_w(i,j,k))) + &
              w(i,j,k)*N_w_square(i,j,k) + &
              theta(i,j,k)*(Us_w(i,j,k)/(z_u(k)*COS(phi_u(j))) *  &
              (dThetas_dl_w(i,j,k)/Thetas_w(i,j,k)-dRho_dl_w(i,j,k)/Rho_w(i,j,k)) + &
              Vs_w(i,j,k)/(z_u(k)) *  &
              (dThetas_dphi_w(i,j,k)/Thetas_w(i,j,k)-dRho_dphi_w(i,j,k)/Rho_w(i,j,k)) + &
              Ws_w(i,j,k)* (N_w_square(i,j,k)/gz_w(k) + 2.d0/z_w(k) - & ! derivative of g
              dRho_dr_w(i,j,k)/Rho_w(i,j,k) ) ) + &
              dw_dl_on_w(theta,i,j,k,dlambda)*Us_w(i,j,k)/(z_w(k)*COS(phi_u(j))) + &
              dw_dphi_on_w(theta,i,j,k,dphi)*Vs_w(i,j,k)/z_w(k) + &
              dw_dr_on_w(theta,i,j,k,dz)*Ws_w(i,j,k) + & ! additional terms
              theta(i,j,k)*trad_theta(i,j,k) + &
              p_on_w(p,i,j,k) * gz_w(k) *Rho_w(i,j,k)* trad_theta(i,j,k) * &
              kappa / Ps_w(i,j,k) !+ &
              !Q_theta(i,j,k) * kappa*p_w(i,j,k) / Ps_w(i,j,k)



  	      CALL PDELSET(mat_evol,t,pos,desc_mat,alpha)
  	      t=t+1
  	    END DO
  	  END DO
  	END DO

  	RETURN

  END SUBROUTINE


  !---------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------
  ! Subroutine to set the radiative timescale
  !---------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------

  SUBROUTINE set_trad(trad_theta, trad_p, Ps_w, Ps_p)

    DOUBLE PRECISION, INTENT(INOUT) :: trad_theta(nlong,nlat,0:nz), &
    trad_p(nlong,nlat,nz)

    DOUBLE PRECISION, INTENT(IN) :: Ps_w(nlong,nlat,0:nz), &
    Ps_p(nlong,nlat,nz)

    DOUBLE PRECISION ::logP, P_high, P_low, trad_bottom


    INTEGER :: i,j,k

    IF (trad_type=='iroetal') THEN
      P_high=1.0E6
      P_low=10.0
! First for theta, then for p
      DO i=1,nlong
        DO j=1,nlat
          DO k=0,nz
            IF (Ps_w(i,j,k)<P_high) THEN
              IF (Ps_w(i,j,k)<P_low) THEN
                logP=LOG10(P_low/1.0E5)
              ELSE
                logP=LOG10(Ps_w(i,j,k)/1.0E5)
              END IF
              trad_theta(i,j,k)=1.0/(10**(5.4659686+1.4940124*logP+ &
              0.66079196*(logP**2)+0.16475329*(logP**3)+ &
              0.014241552*(logP**4)) )
            ELSE
              trad_theta(i,j,k)=0.0
            END IF
          END DO

          DO k=1,nz
            IF (Ps_p(i,j,k)<P_high) THEN
              IF (Ps_p(i,j,k)<P_low) THEN
                logP=LOG10(P_low/1.0E5)
              ELSE
                logP=LOG10(Ps_p(i,j,k)/1.0E5)
              END IF
              trad_p(i,j,k)=1.0/(10**(5.4659686+1.4940124*logP+ &
              0.66079196*(logP**2)+0.16475329*(logP**3)+ &
              0.014241552*(logP**4)) )
            ELSE
              trad_p(i,j,k)=0.0
            END IF
          END DO

        END DO
      END DO

    ELSE IF (trad_type == 'komacek') THEN
    ! In that case, be careful that trad is actually in s-1 and
    ! trad bottom in seconds.
      P_high=1.0E6
      P_low=1000.0
      trad_bottom = 1.0E7
      DO i=1,nlong
        DO j=1,nlat
          DO k=0,nz
           IF (Ps_w(i,j,k)<P_low) THEN
             trad_theta(i,j,k) = trad
           ELSE IF (Ps_w(i,j,k)<P_high) THEN
             trad_theta(i,j,k) =1.0/ ( trad_bottom * (Ps_w(i,j,k) / P_high) ** &
             (LOG((1.0/trad)/trad_bottom)/LOG(P_low/P_high)) )
           ELSE
             trad_theta(i,j,k) =1.0 / trad_bottom
           END IF
         END DO

         DO k=1,nz
           IF (Ps_p(i,j,k)<P_low) THEN
             trad_p(i,j,k) = trad
           ELSE IF (Ps_p(i,j,k)<P_high) THEN
             trad_p(i,j,k) = 1.0/ ( trad_bottom * (Ps_p(i,j,k) / P_high) ** &
             (LOG((1.0/trad)/trad_bottom)/LOG(P_low/P_high)) )
           ELSE
             trad_p(i,j,k) = 1.0/ trad_bottom
           END IF
         END DO

       END DO
     END DO

    ELSE IF(trad_type == 'constnt') THEN
      DO i=1,nlong
        DO j=1,nlat
          DO k=0,nz
            trad_theta(i,j,k) = trad
          END DO

          DO k=1,nz
            trad_p(i,j,k) = trad
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE






  SUBROUTINE set_tdrag(tdrag_u, tdrag_v,tdrag_w,Ps_u,Ps_v, Ps_w)

    DOUBLE PRECISION, INTENT(INOUT) :: tdrag_u(nlong,nlat,nz), &
    tdrag_v(nlong,0:nlat,nz),tdrag_w(nlong,nlat,0:nz)


    DOUBLE PRECISION, INTENT(IN) :: Ps_u(nlong,nlat,nz), &
    Ps_v(nlong,0:nlat,nz) , Ps_w(nlong,nlat,0:nz)


    DOUBLE PRECISION ::logP, P_high, P_low, tdrag_bottom


    INTEGER :: i,j,k


    IF (tdrag_type == 'komacek') THEN
      P_low = 1.0E6
      P_high = p0
      tdrag_bottom = 1.0E-6
      DO i=1,nlong
        DO j=1,nlat
          DO k=1,nz
            tdrag_u(i,j,k) = MAX(tdrag,tdrag_bottom * (Ps_u(i,j,k)-P_low)/(P_high-P_low))
            tdrag_v(i,j,k) = MAX(tdrag,tdrag_bottom * (Ps_v(i,j,k)-P_low)/(P_high-P_low))
            tdrag_w(i,j,k) = MAX(tdrag,tdrag_bottom * (Ps_w(i,j,k)-P_low)/(P_high-P_low))
          END DO
          tdrag_w(i,j,0) = MAX(tdrag,tdrag_bottom * (Ps_w(i,j,0)-P_low)/(P_high-P_low))
        END DO
        DO k=1,nz
          tdrag_v(i,0,k) = MAX(tdrag,tdrag_bottom * (Ps_v(i,0,k)-P_low)/(P_high-P_low))
        END DO
      END DO

    ELSE IF (tdrag_type == 'constnt' .OR. tdrag_type == 'iro' ) THEN
      DO i=1,nlong
        DO j=1,nlat
          DO k=1,nz
            tdrag_u(i,j,k) = tdrag
            tdrag_v(i,j,k) = tdrag
            tdrag_w(i,j,k) = tdrag
          END DO
          tdrag_w(i,j,0) =   tdrag
        END DO
        DO k=1,nz
          tdrag_v(i,0,k) = tdrag
        END DO
      END DO

    END IF

  END SUBROUTINE
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! INTERPOLATION FUNCTIONS
!--------------------------------------------------------------------
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! For u
!--------------------------------------------------------------------

  FUNCTION u_on_v(u,i,j,k)
    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k

    INTEGER :: t

    DOUBLE PRECISION :: u_on_v

    u_on_v=0.0
    ! At the pole, we average taking into account the direction thus cosinus
    IF (j==nlat) THEN
      DO t=1,nlong
        u_on_v=u_on_v+u(t,j,k)*COS((lambda_min+(t-1)*360/nlong)*pi/180.0)/nlong
      END DO

    ELSE IF (j==0) THEN
      IF (i==nlong) THEN
        u_on_v=0.25*(2.5*u(i,1,k)-0.5*u(i,2,k)+ &
        2.5*u(1,1,k)-0.5*u(1,2,k))
      ELSE
        u_on_v=0.25*(2.5*u(i,1,k)-0.5*u(i,2,k)+ &
        2.5*u(i+1,1,k)-0.5*u(i+1,2,k))
      END IF

    ELSE IF (i==nlong) THEN
      u_on_v=0.25*(u(i,j,k)+u(i,j+1,k)+u(1,j,k)+u(1,j+1,k))

    ELSE
      u_on_v=0.25*(u(i,j,k)+u(i,j+1,k)+u(i+1,j,k)+u(i+1,j+1,k))
    END IF

  RETURN
  END FUNCTION

  FUNCTION u_on_w(u,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k

    DOUBLE PRECISION :: u_on_w



    IF (k==0) THEN
      IF (i==nlong) THEN
        u_on_w=0.25*(2.5*u(i,j,1)-0.5*u(i,j,2)+ &
        2.5*u(1,j,1)-0.5*u(1,j,2))
      ELSE
        u_on_w=0.25*(2.5*u(i,j,1)-0.5*u(i,j,2)+ &
        2.5*u(i+1,j,1)-0.5*u(i+1,j,2))
      END IF

    ELSE IF (k==nz) THEN
      IF (i==nlong) THEN
        u_on_w=0.25*(2.5*u(i,j,nz)-0.5*u(i,j,nz-1)+ &
        2.5*u(1,j,nz)-0.5*u(i,j,nz-1))
      ELSE
        u_on_w=0.25*(2.5*u(i,j,nz)-0.5*u(i,j,nz-1)+ &
        2.5*u(i+1,j,nz)-0.5*u(i+1,j,nz-1))
      END IF
    ELSE
      IF (i==nlong) THEN
        u_on_w=0.25*(u(i,j,k)+u(i,j,k+1)+u(1,j,k)+u(1,j,k+1))
      ELSE
        u_on_w=0.25*(u(i,j,k)+u(i,j,k+1)+u(i+1,j,k)+u(i+1,j,k+1))
      END IF
    END IF

   RETURN
  END FUNCTION

  FUNCTION u_on_p(u,i,j,k)
    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    INTEGER, INTENT(IN) ::i,j,k
    DOUBLE PRECISION :: u_on_p

    IF (i==nlong) THEN
      u_on_p=0.5*(u(i,j,k)+u(1,j,k))
    ELSE
      u_on_p=0.5*(u(i,j,k)+u(i+1,j,k))
    END IF

    RETURN

  END FUNCTION

  FUNCTION du_dl_on_u(u,i,j,k,dlambda)

    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dlambda
    DOUBLE PRECISION :: du_dl_on_u

    IF (i==1) THEN
      du_dl_on_u=0.5*(1.5*u(i+1,j,k)-1.5*u(i,j,k))/dlambda
    ELSE IF (i==nlong) THEN
      du_dl_on_u=0.5*(1.5*u(i,j,k)-1.5*u(i-1,j,k))/dlambda
    ELSE
      du_dl_on_u=0.5*(u(i+1,j,k)-u(i-1,j,k))/dlambda
    END IF


    RETURN
  END FUNCTION

  FUNCTION du_dl_on_w(u,i,j,k,dlambda)

    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dlambda
    DOUBLE PRECISION :: du_dl_on_w


    IF (k==0) THEN
      IF (i==nlong) THEN
        du_dl_on_w=0.5*(2.5*u(1,j,1)-2.5*u(i,j,1) + &
        0.5*u(i,j,2)-0.5*u(1,j,2))/dlambda
      ELSE
        du_dl_on_w=0.5*(2.5*u(i+1,j,1)-2.5*u(i,j,1) + &
        0.5*u(i,j,2)-0.5*u(i+1,j,2))/dlambda
      END IF

    ELSE IF (k==nz) THEN
      IF (i==nlong) THEN
        du_dl_on_w=0.5*(2.5*u(1,j,nz)-2.5*u(i,j,nz) + &
        0.5*u(i,j,nz-1)-0.5*u(1,j,nz-1))/dlambda
      ELSE
        du_dl_on_w=0.5*(2.5*u(i+1,j,nz)-2.5*u(i,j,nz) + &
        0.5*u(i,j,nz-1)-0.5*u(i+1,j,nz-1))/dlambda
      END IF
    ELSE
      IF (i==nlong) THEN
        du_dl_on_w=0.5*(u(1,j,k+1)-u(i,j,k+1) + &
        u(1,j,k)-u(i,j,k))/dlambda
      ELSE
        du_dl_on_w=0.5*(u(i+1,j,k+1)-u(i,j,k+1) + &
        u(i+1,j,k)-u(i,j,k))/dlambda
      END IF
    END IF

    RETURN
  END FUNCTION

  FUNCTION du_dl_on_p(u,i,j,k,dlambda)
    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dlambda

    DOUBLE PRECISION :: du_dl_on_p


    IF (i==nlong) THEN
      du_dl_on_p=(u(1,j,k)-u(i,j,k))/dlambda
    ELSE
      du_dl_on_p=(u(i+1,j,k)-u(i,j,k))/dlambda
    END IF

    RETURN

  END FUNCTION

  FUNCTION du_dphi_on_u(u,i,j,k,dphi)

    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    INTEGER, INTENT(IN) ::i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dphi

    DOUBLE PRECISION :: du_dphi_on_u

    IF (j==1) THEN
      du_dphi_on_u=0.5*(1.5*u(i,j+1,k)-1.5*u(i,j,k))/dphi
    ELSE IF (j==nlat) THEN
      du_dphi_on_u=0.5*(1.5*u(i,j,k)-1.5*u(i,j-1,k))/dphi
    ELSE
      du_dphi_on_u=0.5*(u(i,j+1,k)-u(i,j-1,k))/dphi
    END IF
   RETURN
  END FUNCTION

  FUNCTION du_dr_on_u(u,i,j,k,dr)

    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dr

    DOUBLE PRECISION :: du_dr_on_u

    IF (k==1) THEN
      du_dr_on_u=0.5*(1.5*u(i,j,k+1)-1.5*u(i,j,k))/dr
    ELSE IF (k==nz) THEN
      du_dr_on_u=0.5*(1.5*u(i,j,k)-1.5*u(i,j,k-1))/dr
    ELSE
      du_dr_on_u=0.5*(u(i,j,k+1)-u(i,j,k-1))/dr
    END IF

   RETURN
  END FUNCTION



!--------------------------------------------------------------------
! For v
!--------------------------------------------------------------------

  FUNCTION v_on_u(v,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k

    DOUBLE PRECISION :: v_on_u

    IF (i==1) THEN
      v_on_u=0.25*(v(nlong,j-1,k)+v(nlong,j,k)+v(1,j-1,k)+v(1,j,k))
    ELSE
      v_on_u=0.25*(v(i-1,j-1,k)+v(i-1,j,k)+v(i,j-1,k)+v(i,j,k))
    END IF

    RETURN
  END FUNCTION

  FUNCTION v_on_w(v,i,j,k)
    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k

    DOUBLE PRECISION :: v_on_w

    IF (k==0) THEN
      v_on_w=0.25*(2.5*v(i,j-1,1)-0.5*v(i,j-1,2)+ &
            2.5*v(i,j,1)-0.5*v(i,j,2))
    ELSE IF (k==nz) THEN
      v_on_w=0.25*(2.5*v(i,j-1,nz)-0.5*v(i,j-1,nz-1)+2.5*v(i,j,nz)-0.5*v(i,j,nz-1))
    ELSE
      v_on_w=0.25*(v(i,j-1,k)+v(i,j,k)+v(i,j-1,k+1)+v(i,j,k+1))
    END IF

    RETURN
  END FUNCTION

  FUNCTION v_on_p(v,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k

    DOUBLE PRECISION :: v_on_p

    v_on_p=0.5*(v(i,j-1,k)+v(i,j,k))
    RETURN
  END FUNCTION

  FUNCTION dv_dl_on_v(v,i,j,k,dlambda)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dlambda

    DOUBLE PRECISION :: dv_dl_on_v

    IF (i==1) THEN
      dv_dl_on_v=0.5*(1.5*v(i+1,j,k)-1.5*v(i,j,k))/dlambda
    ELSE IF(i==nlong) THEN
      dv_dl_on_v=0.5*(1.5*v(i,j,k)-1.5*v(i-1,j,k))/dlambda
    ELSE
      dv_dl_on_v=0.5*(v(i+1,j,k)-v(i-1,j,k))/dlambda
    END IF
    RETURN
  END FUNCTION

  FUNCTION dv_dphi_on_v(v,i,j,k,dphi)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dphi

    DOUBLE PRECISION :: dv_dphi_on_v

    IF (j==0) THEN
      dv_dphi_on_v=0.5*(1.5*v(i,j+1,k)-1.5*v(i,j,k))/dphi
    ELSE IF (j==nlat) THEN
      dv_dphi_on_v=0.5*(1.5*v(i,j,k)-1.5*v(i,j-1,k))/dphi
    ELSE
      dv_dphi_on_v=0.5*(v(i,j+1,k)-v(i,j-1,k))/dphi
    END IF

    RETURN
  END FUNCTION

  FUNCTION dv_dphi_on_w(v,i,j,k,dphi)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dphi

    DOUBLE PRECISION :: dv_dphi_on_w


    IF (k==0) THEN
      dv_dphi_on_w=0.5*(2.5*v(i,j,1)-2.5*v(i,j-1,1) + &
        0.5*v(i,j-1,2)-0.5*v(i,j,2))/dphi
    ELSE IF (k==nz) THEN
      dv_dphi_on_w=0.5*(2.5*v(i,j,nz)-2.5*v(i,j-1,nz) + &
        0.5*v(i,j-1,nz-1)-0.5*v(i,j,nz-1))/dphi
    ELSE
      dv_dphi_on_w=0.5*(v(i,j,k+1)-v(i,j-1,k+1) + &
          v(i,j,k)-v(i,j-1,k))/dphi
    END IF

    RETURN
  END FUNCTION

  FUNCTION dv_dphi_on_p(v,i,j,k,dphi)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dphi

    DOUBLE PRECISION :: dv_dphi_on_p

    dv_dphi_on_p=(v(i,j,k)-v(i,j-1,k))/dphi

    RETURN
  END FUNCTION

  FUNCTION dv_dr_on_v(v,i,j,k,dr)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dr

    DOUBLE PRECISION :: dv_dr_on_v

    IF (k==1) THEN
      dv_dr_on_v=0.5*(1.5*v(i,j,k+1)-1.5*v(i,j,k))/dr
    ELSE IF (k==nz) THEN
      dv_dr_on_v=0.5*(1.5*v(i,j,k)-1.5*v(i,j,k-1))/dr
    ELSE
      dv_dr_on_v=0.5*(v(i,j,k+1)-v(i,j,k-1))/dr
    END IF

    RETURN
  END FUNCTION

!--------------------------------------------------------------------
! For w
!--------------------------------------------------------------------

  FUNCTION w_on_u(w,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: w(nlong,nlat,0:nz)
    INTEGER, INTENT(IN) :: i,j,k

    DOUBLE PRECISION :: w_on_u

    IF (i==1) THEN
      w_on_u=0.25*(w(nlong,j,k-1)+w(nlong,j,k)+w(1,j,k-1)+w(1,j,k))
    ELSE
      w_on_u=0.25*(w(i-1,j,k-1)+w(i-1,j,k)+w(i,j,k-1)+w(i,j,k))
    END IF

    RETURN
  END FUNCTION


  FUNCTION w_on_v(w,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: w(nlong,nlat,0:nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION :: w_on_v

    IF (j==0) THEN
      w_on_v=0.25*(2.5*w(i,1,k-1)-0.5*w(i,2,k-1)+ &
        2.5*w(i,1,k)-0.5*w(i,2,k))
    ELSE IF (j==nlat) THEN
      w_on_v=0.25*(2.5*w(i,nlat,k-1)-0.5*w(i,nlat-1,k-1)+ &
        2.5*w(i,nlat,k)-0.5*w(i,nlat-1,k))
    ELSE
      w_on_v=0.25*(w(i,j,k-1)+w(i,j,k)+w(i,j+1,k-1)+w(i,j+1,k))
    END IF
   RETURN
  END FUNCTION

  FUNCTION w_on_p(w,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: w(nlong,nlat,0:nz)
    INTEGER, INTENT(IN) :: i,j,k

    DOUBLE PRECISION :: w_on_p

    w_on_p=0.5*(w(i,j,k-1)+w(i,j,k))
    RETURN

  END FUNCTION


  FUNCTION dw_dl_on_w(w,i,j,k,dlambda)

    DOUBLE PRECISION, INTENT(IN) :: w(nlong,nlat,0:nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dlambda

    DOUBLE PRECISION :: dw_dl_on_w

    IF (i==1) THEN
      dw_dl_on_w=0.5*(1.5*w(i+1,j,k)-1.5*w(i,j,k))/dlambda
    ELSE IF (i==nlong) THEN
      dw_dl_on_w=0.5*(1.5*w(i,j,k)-1.5*w(i-1,j,k))/dlambda
    ELSE
      dw_dl_on_w=0.5*(w(i+1,j,k)-w(i-1,j,k))/dlambda
    END IF

    RETURN
  END FUNCTION

  FUNCTION dw_dphi_on_w(w,i,j,k,dphi)

    DOUBLE PRECISION, INTENT(IN) :: w(nlong,nlat,0:nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dphi

    DOUBLE PRECISION :: dw_dphi_on_w

    IF (j==1) THEN
      dw_dphi_on_w=0.5*(1.5*w(i,j+1,k)-1.5*w(i,j,k))/dphi
    ELSE IF (j==nlat) THEN
      dw_dphi_on_w=0.5*(1.5*w(i,j,k)-1.5*w(i,j-1,k))/dphi
    ELSE
      dw_dphi_on_w=0.5*(w(i,j+1,k)-w(i,j-1,k))/dphi
    END IF

    RETURN
  END FUNCTION

  FUNCTION dw_dr_on_w(w,i,j,k,dz)

    DOUBLE PRECISION, INTENT(IN) :: w(nlong,nlat,0:nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dz

    DOUBLE PRECISION :: dw_dr_on_w

    ! Unfortunately need to take an ugly derivative
    IF (k==0) THEN
      dw_dr_on_w=0.5*(1.5*w(i,j,1)-1.5*w(i,j,0))/dz
    ELSE IF (k==nz) THEN
      dw_dr_on_w=0.5*(1.5*w(i,j,nz)-1.5*w(i,j,nz-1))/dz
    ELSE
      dw_dr_on_w=0.5*(w(i,j,k+1)-w(i,j,k-1))/dz
    END IF

    RETURN
  END FUNCTION

  FUNCTION dw_dr_on_p(w,i,j,k,dz)

    DOUBLE PRECISION, INTENT(IN) :: w(nlong,nlat,0:nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dz

    DOUBLE PRECISION :: dw_dr_on_p

         dw_dr_on_p = (w(i,j,k)-w(i,j,k-1))/dz
    RETURN
  END FUNCTION


  !--------------------------------------------------------------------
! For p
!--------------------------------------------------------------------


  FUNCTION p_on_u(p,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) ::i,j,k

    DOUBLE PRECISION :: p_on_u

    IF (i==1) THEN
      p_on_u=0.5*(p(nlong,j,k)+p(i,j,k))
    ELSE
      p_on_u=0.5*(p(i-1,j,k)+p(i,j,k))
    END IF

  RETURN

  END FUNCTION

 FUNCTION p_on_v(p,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k

    DOUBLE PRECISION :: p_on_v


    IF (j==0) THEN
      p_on_v=0.5*(2.5*p(i,j+1,k)-0.5*p(i,j+2,k))
    ELSE IF (j==nlat) THEN
      p_on_v=0.5*(2.5*p(i,j,k)-0.5*p(i,j-1,k))
    ELSE
      p_on_v=0.5*(p(i,j,k)+p(i,j+1,k))
    END IF

    RETURN
  END FUNCTION

  FUNCTION p_on_w(p,i,j,k)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k

    DOUBLE PRECISION :: p_on_w

    IF (k==0) THEN
      p_on_w=0.5*(2.5*p(i,j,k+1)-0.5*p(i,j,k+2))
    ELSE IF (k==nz) THEN
      p_on_w=0.5*(2.5*p(i,j,k)-0.5*p(i,j,k-1))
    ELSE
      p_on_w=0.5*(p(i,j,k+1)+p(i,j,k))
    END IF

    RETURN
  END FUNCTION

  FUNCTION dp_dl_on_u(p,i,j,k,dlambda)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dlambda


    DOUBLE PRECISION :: dp_dl_on_u

    IF (i==1) THEN
      dp_dl_on_u=(p(1,j,k)-p(nlong,j,k))/dlambda
    ELSE
      dp_dl_on_u=(p(i,j,k)-p(i-1,j,k))/dlambda
    END IF

    RETURN
  END FUNCTION

  FUNCTION dp_dl_on_p(p,i,j,k,dlambda)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dlambda

    DOUBLE PRECISION :: dp_dl_on_p

    IF (i==1) THEN
      dp_dl_on_p=0.5*(1.5*p(i+1,j,k)-1.5*p(i,j,k))/dlambda
    ELSE IF (i==nlong) THEN
      dp_dl_on_p=0.5*(1.5*p(i,j,k)-1.5*p(i-1,j,k))/dlambda
    ELSE
      dp_dl_on_p=0.5*(p(i+1,j,k)-p(i-1,j,k))/dlambda
    END IF

    RETURN
  END FUNCTION

  FUNCTION dp_dphi_on_v(p,i,j,k,dphi)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) ::i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dphi

    DOUBLE PRECISION :: dp_dphi_on_v

    IF (j==0) THEN
      dp_dphi_on_v=0.5*(p(i,2,k)-p(i,1,k))/dphi
    ELSE IF (j==nlat) THEN
      dp_dphi_on_v=0.5*(p(i,nlat,k)-p(i,nlat-1,k))/dphi
    ELSE
      dp_dphi_on_v=(p(i,j+1,k)-p(i,j,k))/dphi
    END IF

    RETURN
  END FUNCTION

  FUNCTION dp_dphi_on_p(p,i,j,k,dphi)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dphi

    DOUBLE PRECISION :: dp_dphi_on_p

    IF (j==1) THEN
      dp_dphi_on_p=0.5*(1.5*p(i,j+1,k)-1.5*p(i,j,k))/dphi
    ELSE IF (j==nlat) THEN
      dp_dphi_on_p=0.5*(1.5*p(i,nlat,k)-1.5*p(i,nlat-1,k))/dphi
    ELSE
      dp_dphi_on_p=0.5*(p(i,j+1,k)-p(i,j-1,k))/dphi
    END IF
    RETURN
  END FUNCTION

  FUNCTION dp_dr_on_p(p,i,j,k,dr)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dr

    DOUBLE PRECISION :: dp_dr_on_p

    IF (k==1) THEN
      dp_dr_on_p=0.5*(1.5*p(i,j,k+1)-1.5*p(i,j,k))/dr
    ELSE IF (k==nz) THEN
      dp_dr_on_p=0.5*(1.5*p(i,j,k)-1.5*p(i,j,k-1))/dr
    ELSE
      dp_dr_on_p=0.5*(p(i,j,k+1)-p(i,j,k-1))/dr
    END IF
    RETURN
  END FUNCTION


  FUNCTION dp_dr_on_w(p,i,j,k,dz)

    DOUBLE PRECISION, INTENT(IN) :: p(nlong,nlat,nz)
    INTEGER, INTENT(IN) :: i,j,k
    DOUBLE PRECISION, INTENT(IN) :: dz

    DOUBLE PRECISION :: dp_dr_on_w

    IF (k==0) THEN
      dp_dr_on_w=0.5*(p(i,j,2)-0.5*p(i,j,1))/dz
    ELSE IF (k==nz) THEN
      dp_dr_on_w=0.5*(p(i,j,nz)-p(i,j,nz-1))/dz
    ELSE
      dp_dr_on_w=(p(i,j,k+1)-p(i,j,k))/dz
    END IF

    RETURN
  END FUNCTION

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! DIFFUSION AND SECOND DERIVATIVE
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  FUNCTION ddu_ddl_on_u(u,dlambda)

    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dlambda

    INTEGER :: i,j,k

    DOUBLE PRECISION :: ddu_ddl_on_u(nlong,nlat,nz)

    DO j=1,nlat
      DO k=1,nz
        DO i=1,nlong
          IF(i==1) THEN
            ddu_ddl_on_u(i,j,k) = (0.5*u(i+1,j,k)-0.5*u(i,j,k))/ &
            (dlambda*dlambda)
          ELSE IF (i==nlong) THEN
            ddu_ddl_on_u(i,j,k) = (0.5*u(i-1,j,k)-0.5*u(i,j,k))/ &
            (dlambda*dlambda)
          ELSE
            ddu_ddl_on_u(i,j,k) = (u(i+1,j,k)+u(i-1,j,k)-2*u(i,j,k))/ &
            (dlambda*dlambda)
          END IF
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION


  FUNCTION ddu_ddphi_on_u(u,dphi)

    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi

    INTEGER :: i,j,k

    DOUBLE PRECISION :: ddu_ddphi_on_u(nlong,nlat,nz)

    DO i=1,nlong
      DO k=1,nz
        DO j=1,nlat
          IF(j==1) THEN
            ddu_ddphi_on_u(i,j,k) = (0.5*u(i,j+1,k)-0.5*u(i,j,k))/ &
            (dphi*dphi)
          ELSE IF (j==nlat) THEN
            ddu_ddphi_on_u(i,j,k) = (0.5*u(i,j-1,k)-0.5*u(i,j,k))/ &
            (dphi*dphi)
          ELSE
            ddu_ddphi_on_u(i,j,k) = (u(i,j+1,k)+u(i,j-1,k)-2*u(i,j,k))/ &
            (dphi*dphi)
          END IF
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION


  FUNCTION ddu_ddr_on_u(u,dr)

    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dr

    INTEGER :: i,j,k

    DOUBLE PRECISION :: ddu_ddr_on_u(nlong,nlat,nz)

    DO i=1,nlong
      DO j=1,nlat
        DO k=1,nz
          IF(k==1) THEN
            ddu_ddr_on_u(i,j,k) = (0.5*u(i,j,k+1)-0.5*u(i,j,k))/ &
            (dr*dr)
          ELSE IF (k==nz) THEN
            ddu_ddr_on_u(i,j,k) = (0.5*u(i,j,k-1)-0.5*u(i,j,k))/ &
            (dr*dr)
          ELSE
            ddu_ddr_on_u(i,j,k) = (u(i,j,k+1)+u(i,j,k-1)-2*u(i,j,k))/ &
            (dr*dr)
          END IF
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION


  FUNCTION ddv_ddl_on_v(v,dlambda)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dlambda

    INTEGER :: i,j,k

    DOUBLE PRECISION :: ddv_ddl_on_v(nlong,0:nlat,nz)

    DO j=0,nlat
      DO k=1,nz
        DO i=1,nlong
          IF(i==1) THEN
            ddv_ddl_on_v(i,j,k) = (0.5*v(i+1,j,k)-0.5*v(i,j,k))/ &
            (dlambda*dlambda)
          ELSE IF (i==nlong) THEN
            ddv_ddl_on_v(i,j,k) = (0.5*v(i-1,j,k)-0.5*v(i,j,k))/ &
            (dlambda*dlambda)
          ELSE
            ddv_ddl_on_v(i,j,k) = (v(i+1,j,k)+v(i-1,j,k)-2*v(i,j,k))/ &
            (dlambda*dlambda)
          END IF
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION


  FUNCTION ddv_ddphi_on_v(v,dphi)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi

    INTEGER :: i,j,k

    DOUBLE PRECISION :: ddv_ddphi_on_v(nlong,0:nlat,nz)

    DO i=1,nlong
      DO k=1,nz
        DO j=0,nlat
          IF(j==0) THEN
            ddv_ddphi_on_v(i,j,k) = (0.5*v(i,j+1,k)-0.5*v(i,j,k))/ &
            (dphi*dphi)
          ELSE IF (j==nlat) THEN
            ddv_ddphi_on_v(i,j,k) = (0.5*v(i,j-1,k)-0.5*v(i,j,k))/ &
            (dphi*dphi)
          ELSE
            ddv_ddphi_on_v(i,j,k) = (v(i,j+1,k)+v(i,j-1,k)-2*v(i,j,k))/ &
            (dphi*dphi)
          END IF
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION


  FUNCTION ddv_ddr_on_v(v,dr)

    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dr

    INTEGER :: i,j,k

    DOUBLE PRECISION :: ddv_ddr_on_v(nlong,0:nlat,nz)

    DO i=1,nlong
      DO j=0,nlat
        DO k=1,nz
          IF(k==1) THEN
            ddv_ddr_on_v(i,j,k) = (0.5*v(i,j,k+1)-0.5*v(i,j,k))/ &
            (dr*dr)
          ELSE IF (k==nz) THEN
            ddv_ddr_on_v(i,j,k) = (0.5*v(i,j,k-1)-0.5*v(i,j,k))/ &
            (dr*dr)
          ELSE
            ddv_ddr_on_v(i,j,k) = (v(i,j,k+1)+v(i,j,k-1)-2*v(i,j,k))/ &
            (dr*dr)
          END IF
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION





  FUNCTION diffusion_u(u,ddu_ddl,du_dphi,ddu_ddphi,du_dr,ddu_ddr,&
           Rho_u, dRho_dr_u, ddRho_ddr_u,phi_u,z_u)

    DOUBLE PRECISION, DIMENSION(nlong,nlat,nz), INTENT(IN) :: u, ddu_ddl, &
    du_dphi,ddu_ddphi,du_dr,ddu_ddr, Rho_u, dRho_dr_u, ddRho_ddr_u
    DOUBLE PRECISION, INTENT(IN) :: phi_u(nlat), z_u(nz)

    INTEGER :: i,j,k

    DOUBLE PRECISION :: diffusion_u(nlong,nlat,nz)


    DO i=1, nlong
      DO j=1,nlat
        DO k=1,nz
          diffusion_u(i,j,k) = 1./COS(phi_u(j)) * (ddu_ddl(i,j,k) - &
          SIN(phi_u(j))*du_dphi(i,j,k) + COS(phi_u(j))*ddu_ddphi(i,j,k) + &
          u(i,j,k)* (-2.* (rtot+z_u(k))/rho_u(i,j,k) * dRho_dr_u(i,j,k) - 2 + &
          2.* ((rtot+z_u(k))/rho_u(i,j,k) * dRho_dr_u(i,j,k))**2 - &
          (rtot+z_u(k))**2./rho_u(i,j,k) * ddRho_ddr_u(i,j,k) ) + &
          du_dr(i,j,k)*(2.*(rtot+z_u(k))- &
          2* (rtot+z_u(k))**2/rho_u(i,j,k) *dRho_dr_u(i,j,k)) + &
          (rtot+z_u(k))**2*ddu_ddr(i,j,k) )
        END DO
      END DO
    END DO

    RETURN

  END FUNCTION


  FUNCTION diffusion_v(v,ddv_ddl,dv_dphi,ddv_ddphi,dv_dr,ddv_ddr,&
           Rho_v, dRho_dr_v,ddRho_ddr_v,phi_v,z_u)

    DOUBLE PRECISION, DIMENSION(nlong,0:nlat,nz), INTENT(IN) :: v, ddv_ddl, &
    dv_dphi,ddv_ddphi,dv_dr,ddv_ddr, Rho_v, dRho_dr_v, ddRho_ddr_v
    DOUBLE PRECISION, INTENT(IN) :: phi_v(0:nlat), z_u(nz)

    INTEGER :: i,j,k

    DOUBLE PRECISION :: diffusion_v(nlong,0:nlat,nz)


    DO i=1, nlong
      DO j=0,nlat
        DO k=1,nz
          diffusion_v(i,j,k) = 1/COS(phi_v(j)) * (ddv_ddl(i,j,k) - &
          SIN(phi_v(j))*dv_dphi(i,j,k) + COS(phi_v(j))*ddv_ddphi(i,j,k) + &
          v(i,j,k)* (-2.* (rtot+z_u(k))/rho_v(i,j,k) * dRho_dr_v(i,j,k) - 2 + &
          2.* ((rtot+z_u(k))/rho_v(i,j,k) * dRho_dr_v(i,j,k))**2 - &
          (rtot+z_u(k))**2./rho_v(i,j,k) * ddRho_ddr_v(i,j,k) ) + &
          dv_dr(i,j,k)*(2.*(rtot+z_u(k))- &
          2* (rtot+z_u(k))**2/rho_v(i,j,k) *dRho_dr_v(i,j,k)) + &
          (rtot+z_u(k))**2*ddv_ddr(i,j,k) )
       END DO
      END DO
    END DO

    RETURN

  END FUNCTION




END MODULE
