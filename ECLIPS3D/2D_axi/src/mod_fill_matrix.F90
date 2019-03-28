! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE mod_fill_matrix


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
  SUBROUTINE fill_matrix(ntab,ndtab)
    
    IMPLICIT NONE
    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Subroutine arguments
    !----------------------------------------------------------
    !----------------------------------------------------------

    !points for matrix inversion

    INTEGER, INTENT(IN) :: ntab
    INTEGER, INTENT(IN) :: ndtab

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dataarray ! name of steady state
    !data array (U,V,P,W,Theta). U, V and P having nz vertical points 
    !and W and Theta having nz+1 thus ths size (see Thuburn et al. 2002)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: derivativearray ! Steady state derivative array used here
    ! name of output array for perturbed quantities (u,v,p,w,theta) for each
    ! eigenfrequency. Because of boundary condition w prime and theta must 
    ! be zero at top and botttom so no need to include them in the inversion
    ! : thus the  minus 2. The plus 2 in the second dimension is for im 
    ! and real part of frequency.
                          
    !----------------------------------------------------------
    !----------------------------------------------------------                     
    !Local Variables
    !----------------------------------------------------------
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! Steady state variables (see grid)
    !----------------------------------------------------------
    DOUBLE PRECISION :: T_s(nlat,nz), P_s(nlat,nz), Theta_s(nlat,nz), &
            U_s(nlat,nz), V_s(nlat,nz), W_s(nlat,nz), C_s_square(nlat,nz), &
            N_s_square(nlat,nz), Exner_s(nlat,nz), Rho_s(nlat,nz) 
            ! Full z level steady state plus sound speed and buoyancy freq
            ! Full level correspond to U,p coordinates
         
    DOUBLE PRECISION :: T_hs(nlat,0:nz), P_hs(nlat,0:nz), Theta_hs(nlat,0:nz), &
            U_hs(nlat,0:nz), V_hs(nlat,0:nz), W_hs(nlat,0:nz), &
            C_hs_square(nlat,0:nz), N_hs_square(nlat,0:nz), &
            Exner_hs(nlat,0:nz), Rho_hs(nlat,0:nz) 
            ! half z level steady state
            ! correspond to theta,W coordinate
            
    DOUBLE PRECISION :: T_ls(0:nlat,nz), P_ls(0:nlat,nz), Theta_ls(0:nlat,nz), &
            U_ls(0:nlat,nz), V_ls(0:nlat,nz),W_ls(0:nlat,nz), &
            Rho_ls(0:nlat,nz) 
            ! half latitude steady state
            ! Correspond to V coordinate, no need of N and C 
            
    DOUBLE PRECISION :: dU_dphi_s(nlat,nz),dU_dr_s(nlat,nz), dV_dphi_s(nlat,nz), &
            dV_dphi_hs(nlat,0:nz), dV_dphi_ls(0:nlat,nz),dV_dr_ls(0:nlat,nz),dP_dphi_s(nlat,nz), &
            dP_dphi_hs(nlat,0:nz), dP_dphi_ls(0:nlat,nz),dP_dr_hs(nlat,0:nz), &
            dTheta_dr_s(nlat,nz),dTheta_dphi_hs(nlat,0:nz),dTheta_dr_hs(nlat,0:nz), &
            dRho_dphi_s(nlat,nz),dRho_dphi_hs(nlat,0:nz),dRho_dphi_ls(0:nlat,nz), &
            dRho_dr_s(nlat,nz), dRho_dr_hs(nlat,0:nz), dRho_dr_ls(0:nlat,nz), &
            dW_dr_s(nlat,nz), dW_dphi_hs(nlat,0:nz), dW_dr_hs(nlat,0:nz)
    
    DOUBLE PRECISION :: P_high, P_low, logP,trad_invert(nlat,0:nz), &
    tdrag_u(nlat,nz), tdrag_v(0:nlat,nz), tdrag_w(nlat,0:nz)           
               
    !----------------------------------------------------------
    ! Results variables
    !----------------------------------------------------------
                               
    DOUBLE PRECISION :: freq_r(ntot), freq_i(ntot), eigenvec_r(ntot,ntot), &
                           eigenvec_i(ntot,ntot) ! Frequency and 
    !eigenvector, real and imaginary part
    DOUBLE PRECISION:: RWORK(2*ntot) ! Used for inversion routine
    
    COMPLEX*16 :: mat_solve(ntot,ntot) , eigen_left(ntot,ntot), &
        eigen_right(ntot,ntot), frequencies(ntot) ! Matrix and eigenvectors
        ! Eigenvectors are same size than matrix because one for each frequency
    INTEGER :: LWORK ! Useful integer for matrix inversion
    COMPLEX*16, ALLOCATABLE :: WORK(:) ! Used for the inversion     
    
    
    COMPLEX*16 :: u(nlat,nz), v(0:nlat,nz), p(nlat,nz), theta(nlat,1:nz-1), &
            w(nlat,1:nz-1) ! Perturbed variables
            
    DOUBLE PRECISION :: pp(nlat,nz),uv(nlat,nz), timescale
            
    !----------------------------------------------------------
    ! Other variables
    !----------------------------------------------------------
    
    INTEGER :: zero_lat, zero_z ! Number of zeros in height and latitude
    INTEGER :: i,j,k,t ! Loops integers 
    
    DOUBLE PRECISION :: z(nz), z_h(0:nz), gz(nz), gz_h(0:nz), phi(nlat), phi_hl(0:nlat) 
    ! Height,g and latitude full and half
    DOUBLE PRECISION :: dz, dphi ! height and latitude steps
    COMPLEX*16 :: S(ntot) ! Just used for computing the matrix 
        

    !----------------------------------------------------------
    !----------------------------------------------------------                     
    ! Subroutine
    !----------------------------------------------------------
    !----------------------------------------------------------

    
    !----------------------------------------------------------
    ! Initialisation of the variables
    !----------------------------------------------------------
    

   open(unit=1,file=TRIM(DIRDATA) // "data.dat", access= &
    'sequential')

    open(unit=2,file=TRIM(DIRDATA) // "ddata.dat", access= &
    'sequential')

    ALLOCATE(dataarray(ntab))
    ALLOCATE(derivativearray(ndtab))

    READ(1,*) dataarray
    READ(2,*) derivativearray


    !-----------------------------------------------------------
    ! Reading the input array according to grid
    k=1
    
    !------------
    !U
    DO i=1,nlat
      DO j=1,nz
        U_s(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        U_hs(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        U_ls(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    !------------
    !V
    DO i=1,nlat
      DO j=1,nz
        V_s(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        V_hs(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        V_ls(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    !------------
    !P
    DO i=1,nlat
      DO j=1,nz
        P_s(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        P_hs(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        P_ls(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    !------------
    !Theta
    DO i=1,nlat
      DO j=1,nz
        Theta_s(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        Theta_hs(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        Theta_ls(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    !------------
    !Rho
    DO i=1,nlat
      DO j=1,nz
        Rho_s(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        Rho_hs(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        Rho_ls(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    
    !------------
    !W
    
    DO i=1,nlat
      DO j=1,nz
        W_s(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        W_hs(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        W_ls(i,j)=dataarray(k)
        k=k+1
      END DO
    END DO
    
   
    !------------------------------------------------------------------------
    ! Reading the derivative array 
    k=1
    
    !------------
    !dU
    DO i=1,nlat
      DO j=1,nz
        dU_dphi_s(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=1,nz
        dU_dr_s(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    !------------
    !dV
    DO i=1,nlat
      DO j=1,nz
        dV_dphi_s(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO


    DO i=0,nlat
      DO j=1,nz
        dV_dphi_ls(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        dV_dr_ls(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    !------------
    !dP
    
    DO i=1,nlat
      DO j=1,nz
        dP_dphi_s(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    
    DO i=0,nlat
      DO j=1,nz
        dP_dphi_ls(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        dP_dr_hs(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO

    
    !------------
    !dTheta
    
    DO i=1,nlat
      DO j=1,nz
        dTheta_dr_s(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        dTheta_dphi_hs(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        dTheta_dr_hs(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    
    !------------
    !dRho
    DO i=1,nlat
      DO j=1,nz
        dRho_dphi_s(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    
    DO i=1,nlat
      DO j=0,nz
        dRho_dphi_hs(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        dRho_dphi_ls(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=1,nz
        dRho_dr_s(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        dRho_dr_hs(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        dRho_dr_ls(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    !------------
    !dW
    
    DO i=1,nlat
      DO j=0,nz
        dW_dphi_hs(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=1,nz
        dW_dr_s(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
        
    
    DO i=1,nlat
      DO j=0,nz
        dW_dr_hs(i,j)=derivativearray(k)
        k=k+1
      END DO
    END DO
     
        !Defining height and latitude steps
    dz=(height_max-height_min)/nz
    dphi=ymax*pi/180/(nlat) 
    
    ! Defining height and gravity
    DO j=0,nz-1 
      IF (deep) THEN
        z_h(j)=rtot+height_min+j*dz
        z(j+1)=rtot+height_min+(j+0.5)*dz
      ELSE 
        z_h(j)=rtot
        z(j+1)=rtot
      END IF
      
      IF (g_var) THEN
      	gz_h(j)=g*rtot*rtot/(z_h(j)*z_h(j))
      	gz(j+1)=g*rtot*rtot/(z(j+1)*z(j+1))
      ELSE
        gz_h(j)=g
        gz(j+1)=g
      END IF
    END DO
    
    IF (deep) THEN
      z_h(nz)=rtot+height_min+nz*dz
      ELSE 
      z_h(nz)=rtot
    END IF  
    
    IF (g_var) THEN
      gz_h(nz)=g*rtot*rtot/(z_h(nz)*z_h(nz))
    ELSE
      gz_h(nz)=g   
    END IF
    
    !Defining latitude full and half LATERAL levels
    DO i=1,nlat
      phi(i)=(i-0.5)*dphi
      phi_hl(i)=(i)*dphi
    END DO
    phi_hl(0)=0.0
   ! Sound Speed and Brunt Vaisala frequency with c=gamma*R*T
    DO i=1,nlat
      DO j=1,nz
        N_s_square(i,j)=gz(j)*(dTheta_dr_s(i,j))/ &
                     (Theta_s(i,j))
        C_s_square(i,j)=gascons/(1.0-kappa)*(P_s(i,j)/p0)**(kappa)*Theta_s(i,j)
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
       
        N_hs_square(i,j)=gz_h(j)*(dTheta_dr_hs(i,j))/ &
                     (Theta_hs(i,j))
        C_hs_square(i,j)=gascons/(1.0-kappa)*(P_hs(i,j)/p0)**(kappa) &
                      *Theta_hs(i,j)
      END DO
    END DO

    P_high=1.0E6
    P_low=1.0
    timescale=0.0
    DO j=1,nlat
      DO k=0,nz
        IF (P_hs(j,k)<P_high) THEN
          IF (P_hs(j,k)<P_low) THEN
            logP=LOG10(P_low/1.0E5)
          ELSE
            logP=LOG10(P_hs(j,k)/1.0E5)
          END IF
          trad_invert(j,k)=1.0/(10**(5.4659686+1.4940124*logP+ &
          0.66079196*(logP**2)+0.16475329*(logP**3)+ &
          0.014241552*(logP**4)) )
          trad_invert(j,k)=timescale
        ELSE
          trad_invert(j,k)=timescale
        END IF
      END DO
    END DO
    

    DO j=1,nlat
      DO k=1,nz
        IF (P_s(j,k)<P_high) THEN
          IF (P_s(j,k)<P_low) THEN
            logP=LOG10(P_low/1.0E5)
          ELSE
            logP=LOG10(P_s(j,k)/1.0E5)
          END IF
          tdrag_u(j,k)=1.0/(10**(5.4659686+1.4940124*logP+ &
          0.66079196*(logP**2)+0.16475329*(logP**3)+ &
          0.014241552*(logP**4)) )
          tdrag_u(j,k)=timescale/5.0
        ELSE
          tdrag_u(j,k)=timescale/5.0
        END IF
      END DO
    END DO

    DO j=0,nlat
      DO k=1,nz
        IF (P_ls(j,k)<P_high) THEN
          IF (P_ls(j,k)<P_low) THEN
            logP=LOG10(P_low/1.0E5)
          ELSE
            logP=LOG10(P_ls(j,k)/1.0E5)
          END IF
          tdrag_v(j,k)=1.0/(10**(5.4659686+1.4940124*logP+ &
          0.66079196*(logP**2)+0.16475329*(logP**3)+ &
          0.014241552*(logP**4)) )
          tdrag_v(j,k)=timescale/5.0
        ELSE
          tdrag_v(j,k)=timescale/5.0
        END IF
      END DO
    END DO

    DO j=1,nlat
      DO k=0,nz
        IF (P_hs(j,k)<P_high) THEN
          IF (P_hs(j,k)<P_low) THEN
            logP=LOG10(P_low/1.0E5)
          ELSE
            logP=LOG10(P_hs(j,k)/1.0E5)
          END IF
          tdrag_w(j,k)=1.0/(10**(5.4659686+1.4940124*logP+ &
          0.66079196*(logP**2)+0.16475329*(logP**3)+ &
          0.014241552*(logP**4)) )
          tdrag_w(j,k)=timescale/5.0
        ELSE
          tdrag_w(j,k)=timescale/5.0
        END IF
      END DO
    END DO


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
      S(i)=(1.0,0.0)

!      CALL evol_coeff_rest_UMgrid(i,dz, &
!      dphi,S,&
!      z,z_h,N_s_square,N_hs_square,C_s_square,C_hs_square,phi,phi_hl, &
!      gz,gz_h,U_s,U_hs,V_s,V_hs,P_s,P_hs, Theta_s, Theta_hs,Rho_s,Rho_hs, &
 !     U_ls, V_ls, P_ls, Rho_ls, Theta_ls)
   



      

      CALL evol_coeff_UMgrid(i,dz,dphi,S,&
      z,z_h,N_s_square,N_hs_square,C_s_square,C_hs_square,phi,phi_hl, &
       gz,gz_h,U_s,U_hs,V_s,V_hs,P_s,P_hs, Theta_s, Theta_hs,Rho_s,Rho_hs, &
       U_ls, V_ls, P_ls, Rho_ls, Theta_ls, &
       dU_dphi_s,dU_dr_s, dV_dphi_s,dV_dphi_ls,dV_dr_ls,dP_dphi_s, &
       dP_dphi_ls,dP_dr_hs, dTheta_dr_s,dTheta_dphi_hs,dTheta_dr_hs, &
       dRho_dphi_s,dRho_dphi_hs,dRho_dphi_ls, dRho_dr_s, dRho_dr_hs, dRho_dr_ls, &
       dW_dr_s, dW_dphi_hs, dW_dr_hs, W_s, W_hs, W_ls, trad_invert,&
       tdrag_u, tdrag_v, tdrag_w)



    END DO

    IF (myrow==0 .AND. mycol==0) THEN 
      open(unit=44,file= TRIM(DIRDATA) // 'rho_cs_ns.dat', &
      access='SEQUENTIAL') 
     
       WRITE(44,*) nlat,nz,height_max
!     
       WRITE(44,888) Rho_s
       WRITE(44,888) Rho_ls
       WRITE(44,888) Rho_hs
       WRITE(44,888) C_s_square
       WRITE(44,888) N_hs_square
    CLOSE(44)
  END IF


  CLOSE(44)

  888 FORMAT(10E14.6 )

  END SUBROUTINE
    
  SUBROUTINE evol_coeff_UMgrid(pos,dz,dphi,S, &
      z,z_h,N_s_square,N_hs_square,C_s_square,C_hs_square,phi,phi_hl, &
      gz,gz_h,U_s,U_hs,V_s,V_hs,P_s,P_hs, Theta_s, Theta_hs,Rho_s,Rho_hs, &
      U_ls, V_ls, P_ls, Rho_ls, Theta_ls, &
      dU_dphi_s,dU_dr_s, dV_dphi_s,dV_dphi_ls,dV_dr_ls,dP_dphi_s, &
      dP_dphi_ls,dP_dr_hs, dTheta_dr_s,dTheta_dphi_hs,dTheta_dr_hs, &
      dRho_dphi_s,dRho_dphi_hs,dRho_dphi_ls, dRho_dr_s, dRho_dr_hs, dRho_dr_ls,&
      dW_dr_s, dW_dphi_hs, dW_dr_hs, W_s, W_hs, W_ls, trad_invert,&
      tdrag_u, tdrag_v, tdrag_w)    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Function to calculate the evolution of perturbed variables
    !----------------------------------------------------------
    !----------------------------------------------------------
    
    
    IMPLICIT NONE
    
    ! Declaring the inputs already known by sedulous reader
    INTEGER, INTENT(IN) :: pos
    DOUBLE PRECISION, INTENT(IN) :: z(nz), z_h(0:nz), gz(nz), gz_h(0:nz), phi(nlat), &
    phi_hl(0:nlat),P_s(nlat,nz), Theta_s(nlat,nz), &
    U_s(nlat,nz), V_s(nlat,nz), C_s_square(nlat,nz), &
    N_s_square(nlat,nz), Rho_s(nlat,nz),&
    P_hs(nlat,0:nz), Theta_hs(nlat,0:nz),U_hs(nlat,0:nz), V_hs(nlat,0:nz), &
    C_hs_square(nlat,0:nz), N_hs_square(nlat,0:nz), &
    Rho_hs(nlat,0:nz), W_s(nlat,nz), W_hs(nlat,0:nz) 
    
    DOUBLE PRECISION, INTENT(IN):: dU_dphi_s(nlat,nz),dU_dr_s(nlat,nz),dV_dphi_s(nlat,nz), &
    dV_dphi_ls(0:nlat,nz),dV_dr_ls(0:nlat,nz),dP_dphi_s(nlat,nz), &
    dP_dphi_ls(0:nlat,nz),dP_dr_hs(nlat,0:nz), &
    dTheta_dr_s(nlat,nz), dTheta_dphi_hs(nlat,0:nz),dTheta_dr_hs(nlat,0:nz), &
    dRho_dphi_s(nlat,nz),dRho_dphi_hs(nlat,0:nz),dRho_dphi_ls(0:nlat,nz), &
    dRho_dr_s(nlat,nz), dRho_dr_hs(nlat,0:nz), dRho_dr_ls(0:nlat,nz), &
    dW_dr_s(nlat,nz), dW_dphi_hs(nlat,0:nz), &
    dW_dr_hs(nlat,0:nz), trad_invert(nlat,0:nz), tdrag_u(nlat,nz), &
    tdrag_v(0:nlat,nz), tdrag_w(nlat,0:nz)
        
            
    DOUBLE PRECISION, INTENT(IN) :: P_ls(0:nlat,nz), Theta_ls(0:nlat,nz), &
    U_ls(0:nlat,nz), V_ls(0:nlat,nz), Rho_ls(0:nlat,nz), W_ls(0:nlat,nz)
    
            
    COMPLEX*16, INTENT(IN) :: S(ntot)
    DOUBLE PRECISION, INTENT(IN) :: dz,dphi    
        
    COMPLEX*16 :: u(nlat,nz), v(0:nlat,nz), p(nlat,nz), w(nlat,0:nz), &
                     theta(nlat,0:nz) ! perturbed quantities, all zero but one
    
    COMPLEX*16, DIMENSION(nlat,nz) :: du_dphi_u, du_dr_u, v_u, dv_dphi_u,  w_u, &
    dw_dr_u, theta_u, dp_dphi_u, dp_dr_u

    COMPLEX*16, DIMENSION(0:nlat,nz) :: u_v,dv_dphi_v, dv_dr_v, w_v, p_v, &
    dp_dphi_v, theta_v
    
    COMPLEX*16, DIMENSION(nlat,1:nz-1) :: u_w,v_w, dw_dphi_w, dw_dr_w, p_w, &
    dp_dr_w, dtheta_dphi_w, dtheta_dr_w
                     
    COMPLEX*16 :: u_t(nlat,nz), v_t(0:nlat,nz), p_t(nlat,nz), &
                     w_t(nlat,1:nz-1),theta_t(nlat,1:nz-1) ! evolution rate
                     ! of perturbed quantities
                     
    INTEGER :: i,j,k ! just for convenience and loops          
    
    COMPLEX*16 :: ic = (0.0, 1.0) ! complex i      
  	
  	COMPLEX*16 :: temp
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
  	
  	
  	k=1
  	DO i=1,nlat
  	  DO j=1,nz
  	    u(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=0,nlat
  	  DO j=1,nz
  	    v(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz
  	    p(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    w(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    theta(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	!Boundary conditions : always zero
  	w(:,0)=0.0
  	w(:,nz)=0.0
  	theta(:,0)=0.0
  	theta(:,nz)=0.0
  	
   	!----------------------------------------------------------
    ! Assigning the grid variable
    !----------------------------------------------------------
 	
  	du_dphi_u=du_dphi_on_u(u,nlat,nz,dphi)
  	du_dr_u=du_dr_on_u(u,nlat,nz,dz)
  	v_u=v_on_u(v,nlat,nz)
  	dv_dphi_u=dv_dphi_on_u(v,nlat,nz,dphi)
  	w_u=w_on_u(w,nlat,nz)
  	dw_dr_u=dw_dr_on_u(w,nlat,nz,dz)
  	theta_u=w_on_u(theta,nlat,nz)
  	dp_dphi_u=du_dphi_on_u(p,nlat,nz,dphi)
  	dp_dr_u=du_dr_on_u(p,nlat,nz,dz)
  	

  	u_v=u_on_v(u,nlat,nz)
  	dv_dphi_v=dv_dphi_on_v(v,nlat,nz,dphi)
  	dv_dr_v=dv_dr_on_v(v,nlat,nz,dz)
  	w_v=w_on_v(w,nlat,nz)
  	p_v=u_on_v(p,nlat,nz)
  	dp_dphi_v=du_dphi_on_v(p,nlat,nz,dphi)
  	theta_v=w_on_v(theta,nlat,nz)

    u_w=u_on_w(u,nlat,nz)
    v_w=v_on_w(v,nlat,nz)
    dw_dphi_w=dw_dphi_on_w(w,nlat,nz,dphi)
    dw_dr_w=dw_dr_on_w(w,nlat,nz,dz)
    p_w=u_on_w(p,nlat,nz)
    dp_dr_w=du_dr_on_w(p,nlat,nz,dz)
    dtheta_dphi_w=dw_dphi_on_w(theta,nlat,nz,dphi)
    dtheta_dr_w=dw_dr_on_w(theta,nlat,nz,dz)

  	!----------------------------------------------------------
    ! Evolution for each variable
    !----------------------------------------------------------


    !----------------------------------------------------------
    ! u
    !----------------------------------------------------------
    
    DO i=1,nlat
      DO j=1,nz
        u_t(i,j)= -ic * u(i,j) * ( -V_s(i,j)*dRho_dphi_s(i,j)/(z(j)*Rho_s(i,j)) - &
        W_s(i,j)*dRho_dr_s(i,j)/(Rho_s(i,j)) + W_s(i,j)/z(j) - &
        TAN(phi(i))*V_s(i,j)/z(j) ) + &
        m*u(i,j) * U_s(i,j)/(z(j)*COS(phi(i)) ) - &
        ic * du_dphi_u(i,j)*V_s(i,j)/z(j) - &
        ic * du_dr_u(i,j)*W_s(i,j) - &
        ic * v_u(i,j)* ( dU_dphi_s(i,j)/z(j) - 2.0*omega*SIN(phi(i)) - &
        U_s(i,j)*TAN(phi(i))/z(j) ) - &
        ic*w_u(i,j) * (dU_dr_s(i,j) + 2.0*omega*COS(phi(i)) + U_s(i,j)/z(j) ) + &
        m*p(i,j) / (z(j)*COS(phi(i)))! - &
       ! ic*u(i,j)*tdrag_u(i,j)
      END DO
    END DO
  	
    !----------------------------------------------------------
    ! v
    DO i=1,nlat-1
      DO j=1,nz
        v_t(i,j)=-ic * u_v(i,j) * ( 2.0*omega*SIN(phi_hl(i)) + &
        2.0*TAN(phi_hl(i))*U_ls(i,j)/z(j) ) - &
        ic*v(i,j)*(dV_dphi_ls(i,j)/z(j) - &
        V_ls(i,j)*dRho_dphi_ls(i,j)/ (z(j)*Rho_ls(i,j)) - &
        W_ls(i,j)*dRho_dr_ls(i,j)/Rho_ls(i,j) + W_ls(i,j)/z(j) ) + &
        m*v(i,j) * U_ls(i,j)/(z(j)*COS(phi_hl(i))) - &
        ic*dv_dphi_v(i,j) * V_ls(i,j)/z(j) - &
        ic*dv_dr_v(i,j) * W_ls(i,j) - &
        ic*w_v(i,j) * (dV_dr_ls(i,j) + V_ls(i,j)/z(j) ) - &
        ic*p_v(i,j) * (-dP_dphi_ls(i,j) * (1.0-kappa) / (z(j) * P_ls(i,j)) ) - &
        ic*dp_dphi_v(i,j) / z(j) - &
        ic*theta_v(i,j) * (dP_dphi_ls(i,j) / (gz(j) * z(j) * Rho_ls(i,j)) ) !- &
       ! ic*v(i,j)*tdrag_v(i,j)
      END DO
    END DO
  
  
    i=0
    DO j=1,nz
      IF (eq_sym) THEN
        v_t(0,j)=0.0
      ELSE
        v_t(0,j)= - ic*v(i,j)*(dV_dphi_ls(i,j)/z(j) - &
        V_ls(i,j)*dRho_dphi_ls(i,j)/ (z(j)*Rho_ls(i,j)) - &
        W_ls(i,j)*dRho_dr_ls(i,j)/Rho_ls(i,j) + W_ls(i,j)/z(j) ) + &
        m*v(i,j) * U_ls(i,j)/z(j) - &
        ic*dv_dr_v(i,j) * W_ls(i,j) - &
        ic*2.0*p(1,j)/(dphi* z(j)) - &
        ic*v(i,j)*tdrag_v(i,j)
      END IF
    END DO
    
    i=nlat
    
    DO j=1,nz
      IF (m==1) THEN 
      !  v_t(nlat,j)=-ic * u_v(i,j) * ( 2.0*U_ls(i,j)/z(j) ) + &
      !  m*v(i,j) * U_ls(i,j)/(z(j))
         v_t(nlat,j)=0.0
        ! v_t(nlat,j)= 1.5*u(nlat,j)-0.5*u(nlat-1,j)
      ELSE
        v_t(nlat,j)=0.0
      END IF
    END DO
    ! p
    !----------------------------------------------------------    
    DO i=1,nlat
      DO j=1,nz
        p_t(i,j)=m*u(i,j) * C_s_square(i,j)/(z(j)*COS(phi(i))) - &
        ic * v_u(i,j) * (dP_dphi_s(i,j)/(z(j)*Rho_s(i,j)) - &
        C_s_square(i,j)/z(j) * ( dRho_dphi_s(i,j)/Rho_s(i,j) ) ) - &
        ic * C_s_square(i,j)/(z(j)*COS(phi(i))) * ( &
        (COS(phi_hl(i))*v(i,j)-COS(phi_hl(i-1))*v(i-1,j))/dphi) - &
        ic * w_u(i,j) * C_s_square(i,j) * (N_s_square(i,j)/gz(j)) - &
        ic * C_s_square(i,j)*((z_h(j)*z_h(j)*w(i,j)-z_h(j-1)*z_h(j-1)*w(i,j-1)) &
        /(z(j)*z(j)*dz)) - &
        ic * p(i,j) / (1.0-kappa) * ( dV_dphi_s(i,j)/z(j) - &
        V_s(i,j)*TAN(phi(i))/z(j) + dW_dr_s(i,j) + 2.0*W_s(i,j)/z(j) ) + &
        m* p(i,j) * U_s(i,j)/(z(j)* COS(phi(i)))  - &
        ic * dp_dphi_u(i,j) * V_s(i,j)/z(j) - & 
        ic * dp_dr_u(i,j) * W_s(i,j)
      END DO
    END DO
    
    !----------------------------------------------------------
    ! w
    !----------------------------------------------------------    
    DO i=1,nlat
      DO j=1,nz-1
        w_t(i,j)= -ic*u_w(i,j) * (-2.0*omega*COS(phi(i))-2.0*U_hs(i,j)/z_h(j)) - &
        ic * v_w(i,j) * (dW_dphi_hs(i,j)/z_h(j) - 2.0*V_hs(i,j)/z_h(j) ) - &
        ic * w(i,j) * (-V_hs(i,j)*dRho_dphi_hs(i,j)/(z_h(j)*Rho_hs(i,j)) + &
        dW_dr_hs(i,j) - W_hs(i,j)*dRho_dr_hs(i,j)/Rho_hs(i,j)) + &
        m*w(i,j) * U_hs(i,j)/ (z_h(j) * COS(phi(i))) - &
        ic * dw_dphi_w(i,j) * V_hs(i,j)/z_h(j) - &
        ic * dw_dr_w(i,j) * W_hs(i,j) - &
        ic * p_w(i,j) * &! (gz_h(j) / C_hs_square(i,j)) - &
         ( -dP_dr_hs(i,j)/(C_hs_square(i,j)*Rho_hs(i,j)) ) - & ! if non
        ! hydrostatic
        ic * dp_dr_w(i,j) - &
        ic*theta(i,j)* & !(-1)
        (dP_dr_hs(i,j)/(gz_h(j) * Rho_hs(i,j)) ) !- &
        !ic * w(i,j)*tdrag_w(i,j)

      END DO
    END DO
  	!----------------------------------------------------------
    ! theta
    !---------------------------------------------------------- 
    DO i=1,nlat
      DO j=1,nz-1
        theta_t(i,j)= -ic * v_w(i,j) *gz_h(j)*dTheta_dphi_hs(i,j)/(z_h(j)*Theta_hs(i,j)) - &
        ic * w(i,j)*N_hs_square(i,j) - &
        ic * theta(i,j) * ( V_hs(i,j)/z_h(j)*( &
        dTheta_dphi_hs(i,j)/Theta_hs(i,j) - dRho_dphi_hs(i,j)/Rho_hs(i,j) ) + &
        W_hs(i,j) * (N_hs_square(i,j)/gz_h(j) + 2.0/z_h(j) -  &
        dRho_dr_hs(i,j)/Rho_hs(i,j) ) ) + &
        m*theta(i,j)*U_hs(i,j)/(z_h(j)*COS(phi(i))) - &
        ic * dtheta_dphi_w(i,j)*V_hs(i,j)/z_h(j) - &
        ic * dtheta_dr_w(i,j)* W_hs(i,j) !- &  
        !ic * theta(i,j) * trad_invert(i,j)
      END DO
    END DO
  	!----------------------------------------------------------
    ! Output
    !----------------------------------------------------------
  	
  	k=1
  	DO i=1,nlat
  	  DO j=1,nz
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,u_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=0,nlat
  	  DO j=1,nz
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,v_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,p_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,w_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,theta_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  
  END SUBROUTINE



 SUBROUTINE evol_coeff_rest_UMgrid(pos,dz, &
      dphi,S, &
      z,z_h,N_s_square,N_hs_square,C_s_square,C_hs_square,phi,phi_hl, &
      gz,gz_h,U_s,U_hs,V_s,V_hs,P_s,P_hs, Theta_s, Theta_hs,Rho_s,Rho_hs, &
      U_ls, V_ls, P_ls, Rho_ls, Theta_ls) 
    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Function to calculate the evolution of perturbed variables
    !----------------------------------------------------------
    !----------------------------------------------------------
    
    
    IMPLICIT NONE
    
    ! Declaring the inputs already known by sedulous reader
    INTEGER, INTENT(IN) :: pos
    DOUBLE PRECISION, INTENT(IN) :: z(nz), z_h(0:nz), gz(nz), gz_h(0:nz), phi(nlat), &
    phi_hl(0:nlat),P_s(nlat,nz), Theta_s(nlat,nz), &
    U_s(nlat,nz), V_s(nlat,nz), C_s_square(nlat,nz), &
    N_s_square(nlat,nz), Rho_s(nlat,nz),&
    P_hs(nlat,0:nz), Theta_hs(nlat,0:nz),U_hs(nlat,0:nz), V_hs(nlat,0:nz), &
    C_hs_square(nlat,0:nz), N_hs_square(nlat,0:nz), &
    Rho_hs(nlat,0:nz) 
           
    DOUBLE PRECISION, INTENT(IN) :: P_ls(0:nlat,nz), Theta_ls(0:nlat,nz), &
    U_ls(0:nlat,nz), V_ls(0:nlat,nz), Rho_ls(0:nlat,nz) 
    
            
    COMPLEX*16, INTENT(IN) :: S(ntot)
    DOUBLE PRECISION, INTENT(IN) :: dz,dphi
   
    COMPLEX*16 :: u(nlat,nz), v(0:nlat,nz), p(nlat,nz), w(nlat,0:nz), &
                     theta(nlat,0:nz) ! perturbed quantities, all zero but one
    COMPLEX*16 :: u_t(nlat,nz), v_t(0:nlat,nz), p_t(nlat,nz), &
                     w_t(nlat,nz-1),theta_t(nlat,nz-1) ! evolution rate
                     ! of perturbed quantities
                     
    INTEGER :: i,j,k ! just for convenience and loops          
    
    COMPLEX*16 :: ic = (0.0, 1.0) ! complex i      

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
  	k=1
  	DO i=1,nlat
  	  DO j=1,nz
  	    u(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=0,nlat
  	  DO j=1,nz
  	    v(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz
  	    p(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    w(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    theta(i,j)=S(k)
  	    k=k+1
  	  END DO
  	END DO
  	
  	!Boundary conditions : always zero
  	!v(nlat,:)=0.0
  	w(:,0)=0.0
  	w(:,nz)=0.0
  	theta(:,0)=0.0
  	theta(:,nz)=0.0
  	
    DO i=1,nlat
  	  DO j=1,nz
  	    u_t(i,j)=-0.5*(v(i,j)+v(i-1,j))*(2.0*omega*SIN(phi(i))) + &
  	    p(i,j)*m/(z(j)*COS(phi(i))) + &
  	    0.5*(w(i,j-1)+w(i,j))*(2.0*omega*COS(phi(i)))
  	  END DO
  	END DO
  	
  	DO i=1,nlat-1
  	  DO j=1,nz    	    
  	    v_t(i,j)= -0.5*(u(i,j)+u(i+1,j))*( 2.0*omega*SIN(phi_hl(i))) - &
  	      (p(i+1,j)-p(i,j))/(dphi*z(j)) 
  	  END DO
  	END DO
  	
  	i=0
  	DO j=1,nz    	
  	  IF (eq_sym) THEN 
  	    v_t(0,j)=0.0    
  	  ELSE
  	    v_t(0,j)= - 2.0*p(1,j)/(dphi*z(j))
  	  END IF
  	END DO
  	
  	i=nlat 
  	DO j=1,nz    	    
  	  v_t(nlat,j)= 1.5*u(nlat,j)-0.5*u(nlat-1,j)
  	  v_t(nlat,j) = 0.0  
        END DO
  	
  	
  	
  	
  	
  	DO i=1,nlat
  	  DO j=1,nz
  	    p_t(i,j)= C_s_square(i,j)*( 1.0/(z(j)*COS(phi(i)))*( &
  	    m*u(i,j) + &
  	    (COS(phi_hl(i))*v(i,j)-COS(phi_hl(i-1))*v(i-1,j))/dphi  ) + &
  	    (z_h(j)*z_h(j)*w(i,j)-z_h(j-1)*z_h(j-1)*w(i,j-1))/(z(j)*z(j)*dz) + &
  	    0.5*(w(i,j)+w(i,j-1))*N_s_square(i,j)/gz(j) )
  	  END DO
  	END DO

  	  	
  	! Same for w and theta, with height boundary conditions being just zero 
  	! top and bottom which are not even included here
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    w_t(i,j)=omega*(u(i,j)+u(i,j+1))*COS(phi(i))- &
  	    (p(i,j+1)-p(i,j))/(dz) + &
  	    theta(i,j) - &
  	    gz_h(j)*0.5*(p(i,j)+p(i,j+1))/C_hs_square(i,j)
  	    
  	  END DO
  	END DO 
  	
  	
  	DO i=1,nlat
  	  DO j=1,nz-1  	
  	    theta_t(i,j)=N_hs_square(i,j)*w(i,j)
  	  END DO
  	END DO

   
k=1
  	DO i=1,nlat
  	  DO j=1,nz
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,u_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=0,nlat
  	  DO j=1,nz
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,v_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,p_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,w_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    CALL PZELSET(mat_evol,k,pos,desc_mat,theta_t(i,j))
  	    k=k+1
  	  END DO
  	END DO
  
  END SUBROUTINE

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Interpolation from u 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  FUNCTION u_on_v(u,nlat,nz)    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: u(nlat,nz)
    
    COMPLEX*16 :: u_on_v(0:nlat,nz)
    INTEGER :: i,j
    
    DO i=0,nlat
      DO j=1,nz
        IF (i==0) THEN
          u_on_v(0,j)=0.5*(2.5*u(1,j)-0.5*u(2,j))
        ELSE IF (i==nlat) THEN
          u_on_v(nlat,j)=0.5*(2.5*u(nlat,j)-0.5*u(nlat-1,j))
        ELSE 
          u_on_v(i,j)=0.5*(u(i,j)+u(i+1,j))
        END IF
      END DO
    END DO  
    
    RETURN

  END FUNCTION    
 
    
  FUNCTION u_on_w(u,nlat,nz)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: u(nlat,nz)
    
    COMPLEX*16 :: u_on_w(nlat,1:nz-1)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz-1
        u_on_w(i,j)=0.5*(u(i,j)+u(i,j+1))
      END DO
    END DO
  
    RETURN
  END FUNCTION    
    
    
  FUNCTION du_dphi_on_u(u,nlat,nz,dphi)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: u(nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    
    COMPLEX*16 :: du_dphi_on_u(nlat,nz)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz
        IF (i==1) THEN
          du_dphi_on_u(1,j)=0.5*(1.5*u(2,j)-1.5*u(1,j))/dphi
        ELSE IF (i==nlat) THEN
          du_dphi_on_u(nlat,j)=0.5*(1.5*u(nlat,j)-1.5*u(nlat-1,j))/dphi
        ELSE 
          du_dphi_on_u(i,j)=0.5*(u(i+1,j)-u(i-1,j))/dphi
        END IF
      END DO
    END DO
  
    RETURN
  END FUNCTION    
  
  
  FUNCTION du_dr_on_u(u,nlat,nz,dr)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: u(nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    
    COMPLEX*16 :: du_dr_on_u(nlat,nz)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz
        IF (j==1) THEN
          du_dr_on_u(i,1)=0.5*(1.5*u(i,2)-1.5*u(i,1))/dr
        ELSE IF (j==nz) THEN
          du_dr_on_u(i,nz)=0.5*(1.5*u(i,nz)-1.5*u(i,nz-1))/dr
        ELSE 
          du_dr_on_u(i,j)=0.5*(u(i,j+1)-u(i,j-1))/dr
        END IF
      END DO
    END DO
  
    RETURN
  END FUNCTION      
  
  
  FUNCTION du_dphi_on_v(u,nlat,nz,dphi)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: u(nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    
    COMPLEX*16 :: du_dphi_on_v(0:nlat,nz)
    INTEGER :: i,j
    
    DO i=0,nlat
      DO j=1,nz
        IF (i==0) THEN
          du_dphi_on_v(0,j)=0.5*(u(2,j)-u(1,j))/dphi
        ELSE IF (i==nlat) THEN
          du_dphi_on_v(nlat,j)=0.5*(u(nlat,j)-u(nlat-1,j))/dphi
        ELSE 
          du_dphi_on_v(i,j)=(u(i+1,j)-u(i,j))/dphi
        END IF
      END DO
    END DO
  
    RETURN
  END FUNCTION   
  
  FUNCTION du_dr_on_w(u,nlat,nz,dr)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: u(nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    
    COMPLEX*16 :: du_dr_on_w(nlat,1:nz-1)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz-1
        du_dr_on_w(i,j)=(u(i,j+1)-u(i,j))/dr
      END DO
    END DO
  
    RETURN
  END FUNCTION      
   
   
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Interpolation from v 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  FUNCTION v_on_u(v,nlat,nz)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: v(0:nlat,nz)
    
    COMPLEX*16 :: v_on_u(nlat,nz)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz
        v_on_u(i,j)=0.5*(v(i-1,j)+v(i,j))
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
  FUNCTION v_on_w(v,nlat,nz)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: v(0:nlat,nz)
    
    COMPLEX*16 :: v_on_w(nlat,1:nz-1)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz-1
        v_on_w(i,j)=0.25*(v(i-1,j)+v(i,j)+v(i-1,j+1)+v(i,j+1))
      END DO
    END DO
    RETURN
    
  END FUNCTION  

  FUNCTION dv_dphi_on_u(v,nlat,nz,dphi)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: v(0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    
    COMPLEX*16 :: dv_dphi_on_u(nlat,nz)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz
        dv_dphi_on_u(i,j)=(v(i,j)-v(i-1,j))/dphi
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
  
  
  FUNCTION dv_dphi_on_v(v,nlat,nz,dphi)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: v(0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    COMPLEX*16 :: dv_dphi_on_v(0:nlat,nz)
    INTEGER :: i,j
    
    DO i=0,nlat
      DO j=1,nz
        IF (i==0) THEN
          dv_dphi_on_v(0,j)=0.5*(1.5*v(1,j)-1.5*v(0,j))/dphi
        ELSE IF (i==nlat) THEN
          dv_dphi_on_v(nlat,j)=0.5*(1.5*v(nlat,j)-1.5*v(nlat-1,j))/dphi
        ELSE
          dv_dphi_on_v(i,j)=0.5*(v(i+1,j)-v(i-1,j))/dphi
        END IF
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
  FUNCTION dv_dr_on_v(v,nlat,nz,dr)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: v(0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    COMPLEX*16 :: dv_dr_on_v(0:nlat,nz)
    INTEGER :: i,j
    
    DO i=0,nlat
      DO j=1,nz
        IF (j==1) THEN
          dv_dr_on_v(i,j)=0.5*(1.5*v(i,2)-1.5*v(i,1))/dr
        ELSE IF (j==nz) THEN
          dv_dr_on_v(i,j)=0.5*(1.5*v(i,nz)-1.5*v(i,nz-1))/dr
        ELSE
          dv_dr_on_v(i,j)=0.5*(v(i,j+1)-v(i,j-1))/dr
        END IF
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Interpolation from w 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  FUNCTION w_on_u(w,nlat,nz)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: w(nlat,0:nz)
    
    COMPLEX*16 :: w_on_u(nlat,nz)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz
        w_on_u(i,j)=0.5*(w(i,j)+w(i,j-1))
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
  FUNCTION w_on_v(w,nlat,nz)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: w(nlat,0:nz)
    
    COMPLEX*16 :: w_on_v(0:nlat,nz)
    INTEGER :: i,j
    
    DO i=0,nlat
      DO j=1,nz
        IF (i==0) THEN
          w_on_v(0,j)=0.25*(2.5*w(1,j-1)-0.5*w(2,j-1)+2.5*w(1,j)-0.5*w(2,j))
        ELSE IF (i==nlat) THEN 
          w_on_v(nlat,j)=0.25*(2.5*w(nlat,j-1)-0.5*w(nlat-1,j-1)+2.5*w(nlat,j)-0.5*w(nlat-1,j))
        ELSE
          w_on_v(i,j)=0.25*(w(i,j-1)+w(i+1,j-1)+w(i,j)+w(i+1,j))
        END IF
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
  FUNCTION dw_dphi_on_w(w,nlat,nz,dphi)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: w(nlat,0:nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    COMPLEX*16 :: dw_dphi_on_w(nlat,1:nz-1)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz-1
        IF (i==1) THEN
          dw_dphi_on_w(i,j)=0.5*(1.5*w(2,j)-1.5*w(1,j))/dphi
        ELSE IF (i==nlat) THEN
          dw_dphi_on_w(i,j)=0.5*(1.5*w(nlat,j)-1.5*w(nlat-1,j))/dphi
        ELSE
          dw_dphi_on_w(i,j)=0.5*(w(i+1,j)-w(i-1,j))/dphi
        END IF
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
  FUNCTION dw_dr_on_u(w,nlat,nz,dr)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: w(nlat,0:nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    COMPLEX*16 :: dw_dr_on_u(nlat,nz)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz
        dw_dr_on_u(i,j)=(w(i,j)-w(i,j-1))/dr
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
  FUNCTION dw_dr_on_w(w,nlat,nz,dr)
    
    INTEGER, INTENT(IN) :: nlat,nz
    COMPLEX*16, INTENT(IN) :: w(nlat,0:nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    COMPLEX*16 :: dw_dr_on_w(nlat,1:nz-1)
    INTEGER :: i,j
    
    DO i=1,nlat
      DO j=1,nz-1
        dw_dr_on_w(i,j)=0.5*(w(i,j+1)-w(i,j-1))/dr
      END DO
    END DO
  
    RETURN
  END FUNCTION
  
END MODULE    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
