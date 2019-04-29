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
            Rho_ls(0:nlat,nz), C_ls_square(0:nlat,nz)
            ! half latitude steady state
            ! Correspond to V coordinate, no need of N and C 
            
    DOUBLE PRECISION :: dU_dphi_s(nlat,nz),dU_dr_s(nlat,nz), dV_dphi_s(nlat,nz), &
            dV_dphi_hs(nlat,0:nz), dV_dphi_ls(0:nlat,nz),dV_dr_ls(0:nlat,nz),dP_dphi_s(nlat,nz), &
            dP_dphi_hs(nlat,0:nz), dP_dphi_ls(0:nlat,nz),dP_dr_hs(nlat,0:nz), &
            dTheta_dr_s(nlat,nz),dTheta_dphi_hs(nlat,0:nz),dTheta_dr_hs(nlat,0:nz), &
            dRho_dphi_s(nlat,nz),dRho_dphi_hs(nlat,0:nz),dRho_dphi_ls(0:nlat,nz), &
            dRho_dr_s(nlat,nz), dRho_dr_hs(nlat,0:nz), dRho_dr_ls(0:nlat,nz), &
            dW_dr_s(nlat,nz), dW_dphi_hs(nlat,0:nz), dW_dr_hs(nlat,0:nz)
    
    DOUBLE PRECISION :: P_high, P_low, logP,trad_theta(nlat,0:nz),trad_p(nlat,nz), &
    tdrag_u(nlat,nz), tdrag_v(0:nlat,nz), tdrag_w(nlat,0:nz)           
               
    !----------------------------------------------------------
    ! Results variables
    !----------------------------------------------------------

    
    DOUBLE PRECISION :: u(nlat,nz), v(0:nlat,nz), p(nlat,nz), theta(nlat,0:nz), &
            w(nlat,0:nz) ! Perturbed variables

            
    !----------------------------------------------------------
    ! Other variables
    !----------------------------------------------------------
    
    INTEGER :: i,j,pos,t ! Loops integers 
    
    DOUBLE PRECISION :: z(nz), z_h(0:nz), gz(nz), gz_h(0:nz), phi(nlat), phi_hl(0:nlat) 
    ! Height,g and latitude full and half
    DOUBLE PRECISION :: dz, dphi ! height and latitude steps
    DOUBLE PRECISION :: S(ntot) ! Just used for computing the matrix 
    
    COMPLEX*16 :: alpha, ic
        
    ic = (0.d0,1.d0)
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
    t=1
    
    !------------
    !U
    DO i=1,nlat
      DO j=1,nz
        U_s(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        U_hs(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        U_ls(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    !------------
    !V
    DO i=1,nlat
      DO j=1,nz
        V_s(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        V_hs(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        V_ls(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    !------------
    !P
    DO i=1,nlat
      DO j=1,nz
        P_s(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        P_hs(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        P_ls(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    !------------
    !Theta
    DO i=1,nlat
      DO j=1,nz
        Theta_s(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        Theta_hs(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        Theta_ls(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    !------------
    !Rho
    DO i=1,nlat
      DO j=1,nz
        Rho_s(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        Rho_hs(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        Rho_ls(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    
    !------------
    !W
    
    DO i=1,nlat
      DO j=1,nz
        W_s(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        W_hs(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        W_ls(i,j)=dataarray(t)
        t=t+1
      END DO
    END DO
    
   
    !------------------------------------------------------------------------
    ! Reading the derivative array 
    t=1
    
    !------------
    !dU
    DO i=1,nlat
      DO j=1,nz
        dU_dphi_s(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=1,nz
        dU_dr_s(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    !------------
    !dV
    DO i=1,nlat
      DO j=1,nz
        dV_dphi_s(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO


    DO i=0,nlat
      DO j=1,nz
        dV_dphi_ls(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        dV_dr_ls(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    !------------
    !dP
    
    DO i=1,nlat
      DO j=1,nz
        dP_dphi_s(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    
    DO i=0,nlat
      DO j=1,nz
        dP_dphi_ls(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        dP_dr_hs(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO

    
    !------------
    !dTheta
    
    DO i=1,nlat
      DO j=1,nz
        dTheta_dr_s(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        dTheta_dphi_hs(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        dTheta_dr_hs(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    
    !------------
    !dRho
    DO i=1,nlat
      DO j=1,nz
        dRho_dphi_s(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    
    DO i=1,nlat
      DO j=0,nz
        dRho_dphi_hs(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        dRho_dphi_ls(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=1,nz
        dRho_dr_s(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=0,nz
        dRho_dr_hs(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=0,nlat
      DO j=1,nz
        dRho_dr_ls(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    !------------
    !dW
    
    DO i=1,nlat
      DO j=0,nz
        dW_dphi_hs(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    DO i=1,nlat
      DO j=1,nz
        dW_dr_s(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
        
    
    DO i=1,nlat
      DO j=0,nz
        dW_dr_hs(i,j)=derivativearray(t)
        t=t+1
      END DO
    END DO
    
    
!!!!!!!!!!!!! all read


    CALL set_trad(trad_theta, trad_p, P_hs, P_s)
    CALL set_tdrag(tdrag_u, tdrag_v, tdrag_w, P_s, P_ls, P_hs)    
    
     
        !Defining height and latitude steps
    dz=(height_max-height_min)/DBLE(nz)
    dphi=ymax*pi/180.D0/DBLE(nlat) 
    
    ! Defining height and gravity
    DO j=0,nz-1 
      IF (deep) THEN
        z_h(j)=rtot+height_min+DBLE(j)*dz
        z(j+1)=rtot+height_min+(DBLE(j)+0.5d0)*dz
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
      z_h(nz)=rtot+height_min+DBLE(nz)*dz
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
      phi(i)=(DBLE(i)-0.5d0)*dphi
      phi_hl(i)=DBLE(i)*dphi
    END DO
    phi_hl(0)=0.0
   ! Sound Speed and Brunt Vaisala frequency with c=gamma*R*T
    DO i=1,nlat
      DO j=1,nz
        N_s_square(i,j)=gz(j)*(dTheta_dr_s(i,j))/ &
                     (Theta_s(i,j))
        C_s_square(i,j)=gascons/(1.d0-kappa)*(P_s(i,j)/p0)**(kappa)*Theta_s(i,j)
        
        C_ls_square(i,j)=gascons/(1.d0-kappa)*(P_ls(i,j)/p0)**(kappa)*Theta_ls(i,j)
      END DO
    END DO
    
    DO j=1,nz
      C_ls_square(0,j)=gascons/(1.d0-kappa)*(P_ls(0,j)/p0)**(kappa)*Theta_ls(0,j)
    END DO
    
    DO i=1,nlat
      DO j=0,nz
       
        N_hs_square(i,j)=gz_h(j)*(dTheta_dr_hs(i,j))/ &
                     (Theta_hs(i,j))
        C_hs_square(i,j)=gascons/(1.d0-kappa)*(P_hs(i,j)/p0)**(kappa) &
                      *Theta_hs(i,j)
      END DO
    END DO
    
    
    print *, N_hs_square(0,:)
    stop



    !----------------------------------------------------------
    ! Seek for eigenvectors
    !----------------------------------------------------------
                 
    ! To create the matrix, we will just set all the variables to 
    ! 0 at each point of the grid and allow one to be non zero at one point 
    ! of the grid. Then we will calculate the evolution coefficients from the 
    ! perturbed equations due to this value and store it in the matrix from 
    ! which we will get the eigenvectors and values.


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
    

  888 FORMAT(10E14.6 )
  
  END IF

    DO pos=1,ntot
      S=0.d0
      S(pos)= 1.D0
  	!----------------------------------------------------------
    ! Initialisation of the non zero variable
    !----------------------------------------------------------
  	
  	! Just in case
  	u=0.d0
  	v=0.d0
  	p=0.d0
  	w=0.d0
  	theta=0.d0
  	
  	
  	! Reading the stored array
  	
  	
  	t=1
  	DO i=1,nlat
  	  DO j=1,nz
  	    u(i,j)=S(t)
  	    t=t+1
  	  END DO
  	END DO
  	
  	DO i=0,nlat
  	  DO j=1,nz
  	    v(i,j)=S(t)
  	    t=t+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz
  	    p(i,j)=S(t)
  	    t=t+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    w(i,j)=S(t)
  	    t=t+1
  	  END DO
  	END DO
  	
  	DO i=1,nlat
  	  DO j=1,nz-1
  	    theta(i,j)=S(t)
  	    t=t+1
  	  END DO
  	END DO
  	
  	!Boundary conditions : always zero
  	IF (eq_sym) THEN
  	  DO j=1,nz
  	    v(0,j) = 0.d0
  	  END DO
  	END IF
  	
  	DO i=1,nlat
  	  w(:,0)=0.d0
  	  w(:,nz)=0.d0
  	  theta(:,0)=0.d0
  	  theta(:,nz)=0.d0
  	END DO
  	


  	!----------------------------------------------------------
    ! Evolution for each variable
    !----------------------------------------------------------


    !----------------------------------------------------------
    ! u
    !----------------------------------------------------------
 ! --------------------- u evolution ------------------------------------

    t=1

        DO i=1,nlat
          DO j=1,nz
              alpha=u(i,j)*(-1.d0*V_s(i,j)*dRho_dphi_s(i,j)/(z(j)*Rho_s(i,j)) - &
              W_s(i,j)*dRho_dr_s(i,j)/Rho_s(i,j) + W_s(i,j)/z(j) - &
              DTAN(phi(i))*V_s(i,j)/z(j)) + &
              ic*m*u(i,j)*U_s(i,j)/(z(j)*DCOS(phi(i))) + &
              du_dphi_on_u(u,i,j,dphi)*V_s(i,j)/z(j) + &
              du_dr_on_u(u,i,j,dz)*W_s(i,j) + &
              v_on_u(v,i,j)*(dU_dphi_s(i,j)/z(j)-2.d0*omega*DSIN(phi(i))- &
              U_s(i,j)*DTAN(phi(i))/z(j) ) + &
              w_on_u(w,i,j)*(dU_dr_s(i,j)+2.d0*omega*DCOS(phi(i))+U_s(i,j)/z(j)) + &
              ic*m*p(i,j)*(1.d0/(z(j)*DCOS(phi(i)))) + &
              u(i,j)*tdrag_u(i,j)

            CALL PZELSET(mat_evol,t,pos,desc_mat,-ic*alpha)  
            t=t+1
          END DO
        END DO
      
      
      ! --------------------- v evolution ------------------------------------

        DO i=0,nlat
          DO j=1,nz
            IF (i==0) THEN
              alpha=0.d0
            ELSE IF (i==nlat) THEN
              alpha=u_on_v(u,i,j)*(2.d0*omega*DSIN(phi_hl(i))) + &
              v(i,j)*(dV_dphi_ls(i,j)/z(j) - V_ls(i,j)*dRho_dphi_ls(i,j)/ &
              (z(j)*Rho_ls(i,j)) - W_ls(i,j)*dRho_dr_ls(i,j)/Rho_ls(i,j) + &
              W_ls(i,j)/z(j)) + &
              dv_dphi_on_v(v,i,j,dphi)*(V_ls(i,j)/z(j)) + &
              dv_dr_on_v(v,i,j,dz)*W_ls(i,j) + &
              w_on_v(w,i,j)*(dV_dr_ls(i,j)+V_ls(i,j)/z(j)) + &
              u_on_v(p,i,j)*(-dP_dphi_ls(i,j)/(C_ls_square(i,j)*Rho_ls(i,j)*z(j))) + &
              du_dphi_on_v(p,i,j,dphi)*(1.d0/z(j)) + &
              w_on_v(theta,i,j)*(dP_dphi_ls(i,j)/(gz(j)*Rho_ls(i,j)*z(j)))+ &
              v(i,j)*tdrag_v(i,j)
            ELSE
              alpha=u_on_v(u,i,j)*(2.d0*omega*DSIN(phi_hl(i)) + &
              2.d0*DTAN(phi_hl(i))*U_ls(i,j)/z(j)) + &
              v(i,j)*(dV_dphi_ls(i,j)/z(j) - V_ls(i,j)*dRho_dphi_ls(i,j)/ &
              (z(j)*Rho_ls(i,j)) - W_ls(i,j)*dRho_dr_ls(i,j)/Rho_ls(i,j) + &
              W_ls(i,j)/z(j)) + &
              ic*m*v(i,j)*(U_ls(i,j)/(z(j)*DCOS(phi_hl(i)))) + &
              dv_dphi_on_v(v,i,j,dphi)*(V_ls(i,j)/z(j)) + &
              dv_dr_on_v(v,i,j,dz)*W_ls(i,j) + &
              w_on_v(w,i,j)*(dV_dr_ls(i,j)+V_ls(i,j)/z(j)) + &
              u_on_v(p,i,j)*(-dP_dphi_ls(i,j)/(C_ls_square(i,j)*Rho_ls(i,j)*z(j))) + &
              du_dphi_on_v(p,i,j,dphi)*(1.d0/z(j)) + &
              w_on_v(theta,i,j)*(dP_dphi_ls(i,j)/(gz(j)*Rho_ls(i,j)*z(j)))+ &
              v(i,j)*tdrag_v(i,j)
            END IF
            CALL PZELSET(mat_evol,t,pos,desc_mat,-ic*alpha)
            t=t+1
          END DO
        END DO
      
      ! --------------------- p evolution ------------------------------------
        DO i=1,nlat
          DO j=1,nz
            alpha = ic*m*u(i,j)*C_s_square(i,j)/(z(j)*DCOS(phi(i))) + &
            v_on_u(v,i,j) * ((dP_dphi_s(i,j) - C_s_square(i,j)*dRho_dphi_s(i,j))/ &
            (Rho_s(i,j)*z(j)) - C_s_square(i,j)*DTAN(phi(i))/z(j) ) + &
            dv_dphi_on_u(v,i,j,dphi)*C_s_square(i,j)/z(j) + &
            w_on_u(w,i,j)*(C_s_square(i,j) * (2.d0/z(j) + N_s_square(i,j)/gz(j)) ) + &
            dw_dr_on_u(w,i,j,dz)*C_s_square(i,j) + &
            p(i,j)* (1.d0+kappa/(1.d0-kappa)) * (dV_dphi_s(i,j)/z(j) - &
            V_s(i,j)*DTAN(phi(i))/z(j) + &
            dW_dr_s(i,j) + 2.d0*W_s(i,j)/z(j) ) + &
            ic*m*p(i,j)*(U_s(i,j)/ (z(j)*DCOS(phi(i))) ) + &
            du_dphi_on_u(p,i,j,dphi)*V_s(i,j)/z(j) + &
            du_dr_on_u(p,i,j,dz)*W_s(i,j) + & ! Additional terms if Newtonian heating
            w_on_u(theta,i,j) * trad_p(i,j) * C_s_square(i,j)/gz(j) + &
            p(i,j) * trad_p(i,j) * C_s_square(i,j) * Rho_s(i,j) * &
            (kappa) / P_s(i,j)! + &  ! last term for heating function
!             gascons/(1.d0-kappa)*Q_p(i,j)*&
!             (w_on_p(theta,i,j,k)/gz(j)-p(i,j)/C_p_square(i,j))
           CALL PZELSET(mat_evol,t,pos,desc_mat,-ic*alpha)
           t=t+1
          END DO
        END DO
      
      ! --------------------- w evolution ------------------------------------
        DO i=1,nlat
          DO j=1,nz-1
              alpha = u_on_w(u,i,j)*(-2.d0*omega*DCOS(phi(i))-2.d0*U_hs(i,j)/z_h(j)) + &
              v_on_w(v,i,j)*(-2.d0*V_hs(i,j)/z_h(j) + dW_dphi_hs(i,j)/z_h(j) ) + &
              w(i,j)*(-1.D0*V_hs(i,j)*dRho_dphi_hs(i,j)/(z_h(j)*Rho_hs(i,j)) + &
              dW_dr_hs(i,j) - W_hs(i,j)*dRho_dr_hs(i,j)/Rho_hs(i,j) ) + &
              m*w(i,j)*U_hs(i,j)/(z_h(j)*DCOS(phi(i))) + &
              dw_dphi_on_w(w,i,j,dphi)*V_hs(i,j)/z_h(j) + &
              dw_dr_on_w(w,i,j,dz)*W_hs(i,j) + &
              u_on_w(p,i,j)*(-dP_dr_hs(i,j)/(C_hs_square(i,j)*Rho_hs(i,j))) + &
              du_dr_on_w(p,i,j,dz) + &
              theta(i,j)*(dP_dr_hs(i,j)/(Rho_hs(i,j)*gz_h(j))) + &
              w(i,j)*tdrag_w(i,j)
            CALL PZELSET(mat_evol,t,pos,desc_mat,-ic*alpha)
            t=t+1
          END DO
        END DO
      
      
      ! --------------------- theta evolution ------------------------------------
        DO i=1,nlat
          DO j=1,nz-1
            alpha = v_on_w(v,i,j)*(gz_h(j)*dTheta_dphi_hs(i,j)/(z_h(j)*Theta_hs(i,j))) + &
            w(i,j)*N_hs_square(i,j) + &
            theta(i,j)*(V_hs(i,j)/(z(j)) *  &
            (dTheta_dphi_hs(i,j)/Theta_hs(i,j)-dRho_dphi_hs(i,j)/Rho_hs(i,j)) + &
            W_hs(i,j)* (N_hs_square(i,j)/gz_h(j) + 2.d0/z_h(j) - & ! derivative of g 
            dRho_dr_hs(i,j)/Rho_hs(i,j) ) ) + &
            ic*m*theta(i,j)*U_hs(i,j)/(z_h(j)*DCOS(phi(i))) + &
            dw_dphi_on_w(theta,i,j,dphi)*V_hs(i,j)/z_h(j) + &
            dw_dr_on_w(theta,i,j,dz)*W_hs(i,j) + & ! additional terms
            theta(i,j)*trad_theta(i,j) + &
            u_on_w(p,i,j) * gz_h(j) *Rho_hs(i,j)* trad_theta(i,j) * &
            kappa / P_hs(i,j) !+ &
            !p_on_w(p,i,j,k)*Q_theta(i,j)/Thetas_w(i,j) * kappa*gz_h(j)* &
            !Rho_w(i,j)*(p0/Ps_w(i,j))**(kappa)/Ps_w(i,j)

           CALL PZELSET(mat_evol,t,pos,desc_mat,-ic*alpha) 
           t=t+1
          END DO
        END DO

    END DO
    
 
  END SUBROUTINE



  
  
  SUBROUTINE set_trad(trad_theta, trad_p, Ps_w, Ps_p)
    
    DOUBLE PRECISION, INTENT(INOUT) :: trad_theta(nlat,0:nz), &
    trad_p(nlat,nz)

    DOUBLE PRECISION, INTENT(IN) :: Ps_w(nlat,0:nz), &
    Ps_p(nlat,nz) 

    DOUBLE PRECISION ::logP, P_high, P_low, trad_bottom
    
 
    INTEGER :: i,j,k

    IF (trad_type=='iroetal') THEN
      P_high=1.0d6
      P_low=1.d0
! First for theta, then for p
        DO j=1,nlat
          DO k=0,nz
            IF (Ps_w(j,k)<P_high) THEN
              IF (Ps_w(j,k)<P_low) THEN 
                logP=DLOG10(P_low/1.0D5)
              ELSE
                logP=DLOG10(Ps_w(j,k)/1.0D5)
              END IF  
              trad_theta(j,k)=1.d0/(10.d0**(5.4659686d0+1.4940124d0*logP+ &
              0.66079196d0*(logP**2.d0)+0.16475329d0*(logP**3.d0)+ &
              0.014241552d0*(logP**4.d0)) )
            ELSE
              trad_theta(j,k)=0.d0
            END IF
          END DO

          DO k=1,nz
            IF (Ps_p(j,k)<P_high) THEN
              IF (Ps_p(j,k)<P_low) THEN
                logP=DLOG10(P_low/1.0D5)
              ELSE
                logP=DLOG10(Ps_p(j,k)/1.0D5)
              END IF
              trad_p(j,k)=1.d0/(10.d0**(5.4659686d0+1.4940124d0*logP+ &
              0.66079196d0*(logP**2.d0)+0.16475329d0*(logP**3.d0)+ &
              0.014241552d0*(logP**4.d0)) )
            ELSE
              trad_p(j,k)=0.d0
            END IF
          END DO

        END DO

    ELSE IF (trad_type == 'komacek') THEN
    ! In that case, be careful that trad is actually in s-1 and 
    ! trad bottom in seconds.  
      P_high=1.0D6
      P_low=1000.d0
      trad_bottom = 1.0D7
        DO j=1,nlat
          DO k=0,nz
           IF (Ps_w(j,k)<P_low) THEN
             trad_theta(j,k) = trad_invert
           ELSE IF (Ps_w(j,k)<P_high) THEN
             trad_theta(j,k) =1.d0/ ( trad_bottom * (Ps_w(j,k) / P_high) ** &
             (DLOG((1.d0/trad_invert)/trad_bottom)/DLOG(P_low/P_high)) )
           ELSE
             trad_theta(j,k) =1.d0 / trad_bottom
           END IF
         END DO

         DO k=1,nz
           IF (Ps_p(j,k)<P_low) THEN
             trad_p(j,k) = trad_invert
           ELSE IF (Ps_p(j,k)<P_high) THEN
             trad_p(j,k) = 1.d0/ ( trad_bottom * (Ps_p(j,k) / P_high) ** &
             (DLOG((1.d0/trad_invert)/trad_bottom)/DLOG(P_low/P_high)) ) 
           ELSE
             trad_p(j,k) = 1.d0/ trad_bottom  
           END IF
         END DO

       END DO
   
    ELSE IF(trad_type == 'constnt') THEN
        DO j=1,nlat
          DO k=0,nz
            trad_theta(j,k) = trad_invert
          END DO
          
          DO k=1,nz
            trad_p(j,k) = trad_invert
          END DO
        END DO
    END IF 

  END SUBROUTINE






  SUBROUTINE set_tdrag(tdrag_u, tdrag_v,tdrag_w,Ps_u,Ps_v, Ps_w)

    DOUBLE PRECISION, INTENT(INOUT) :: tdrag_u(nlat,nz), &
    tdrag_v(0:nlat,nz),tdrag_w(nlat,0:nz)
    

    DOUBLE PRECISION, INTENT(IN) :: Ps_u(nlat,nz), &
    Ps_v(0:nlat,nz) , Ps_w(nlat,0:nz)
    

    DOUBLE PRECISION ::logP, P_high, P_low, tdrag_bottom


    INTEGER :: i,j,k


    IF (tdrag_type == 'komacek') THEN
      P_low = 1.0D6
      P_high = p0
      tdrag_bottom = 1.0D-6
        DO j=1,nlat
          DO k=1,nz     
            tdrag_u(j,k) = MAX(tdrag_invert,tdrag_bottom * (Ps_u(j,k)-P_low)/(P_high-P_low))
            tdrag_v(j,k) = MAX(tdrag_invert,tdrag_bottom * (Ps_v(j,k)-P_low)/(P_high-P_low))
            tdrag_w(j,k) = MAX(tdrag_invert,tdrag_bottom * (Ps_w(j,k)-P_low)/(P_high-P_low))
          END DO
          tdrag_w(j,0) = MAX(tdrag_invert,tdrag_bottom * (Ps_w(j,0)-P_low)/(P_high-P_low))   
        END DO
        DO k=1,nz
          tdrag_v(0,k) = MAX(tdrag_invert,tdrag_bottom * (Ps_v(0,k)-P_low)/(P_high-P_low))
        END DO
 
    ELSE IF (tdrag_type == 'constnt' .OR. tdrag_type == 'iro' ) THEN
        DO j=1,nlat
          DO k=1,nz
            tdrag_u(j,k) = tdrag_invert
            tdrag_v(j,k) = tdrag_invert
            tdrag_w(j,k) = tdrag_invert
          END DO
          tdrag_w(j,0) =   tdrag_invert
        END DO
        DO k=1,nz
          tdrag_v(0,k) = tdrag_invert
        END DO

    END IF

  END SUBROUTINE

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Interpolation from u 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

FUNCTION u_on_v(u,i,j)    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlat,nz)
    
    DOUBLE PRECISION :: u_on_v


        IF (i==0) THEN
          u_on_v=0.5d0*(2.5d0*u(1,j)-0.5d0*u(2,j))
        ELSE IF (i==nlat) THEN
          u_on_v=0.5d0*(2.5d0*u(nlat,j)-0.5d0*u(nlat-1,j))
        ELSE 
          u_on_v=0.5d0*(u(i,j)+u(i+1,j))
        END IF

    
    RETURN

  END FUNCTION    
 
    
  FUNCTION u_on_w(u,i,j)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlat,nz)
    
    DOUBLE PRECISION :: u_on_w

    

        u_on_w=0.5d0*(u(i,j)+u(i,j+1))

  
    RETURN
  END FUNCTION    
    
    
  FUNCTION du_dphi_on_u(u,i,j,dphi)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    
    DOUBLE PRECISION :: du_dphi_on_u



        IF (i==1) THEN
          du_dphi_on_u=0.5d0*(1.5d0*u(2,j)-1.5d0*u(1,j))/dphi
        ELSE IF (i==nlat) THEN
          du_dphi_on_u=0.5d0*(1.5d0*u(nlat,j)-1.5d0*u(nlat-1,j))/dphi
        ELSE 
          du_dphi_on_u=0.5d0*(u(i+1,j)-u(i-1,j))/dphi
        END IF

  
    RETURN
  END FUNCTION    
  
  
  FUNCTION du_dr_on_u(u,i,j,dr)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    
    DOUBLE PRECISION :: du_dr_on_u

    

        IF (j==1) THEN
          du_dr_on_u=0.5d0*(1.5d0*u(i,2)-1.5d0*u(i,1))/dr
        ELSE IF (j==nz) THEN
          du_dr_on_u=0.5d0*(1.5d0*u(i,nz)-1.5d0*u(i,nz-1))/dr
        ELSE 
          du_dr_on_u=0.5d0*(u(i,j+1)-u(i,j-1))/dr
        END IF

  
    RETURN
  END FUNCTION      
  
  
  FUNCTION du_dphi_on_v(u,i,j,dphi)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    
    DOUBLE PRECISION :: du_dphi_on_v

    

        IF (i==0) THEN
          du_dphi_on_v=0.5d0*(u(2,j)-u(1,j))/dphi
        ELSE IF (i==nlat) THEN
          du_dphi_on_v=0.5d0*(u(nlat,j)-u(nlat-1,j))/dphi
        ELSE 
          du_dphi_on_v=(u(i+1,j)-u(i,j))/dphi
        END IF

  
    RETURN
  END FUNCTION   
  
  FUNCTION du_dr_on_w(u,i,j,dr)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    
    DOUBLE PRECISION :: du_dr_on_w

        du_dr_on_w=(u(i,j+1)-u(i,j))/dr

  
    RETURN
  END FUNCTION      
   
   
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Interpolation from v 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  FUNCTION v_on_u(v,i,j)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(0:nlat,nz)
    
    DOUBLE PRECISION :: v_on_u

    

        v_on_u=0.5d0*(v(i-1,j)+v(i,j))

  
    RETURN
  END FUNCTION
  
  FUNCTION v_on_w(v,i,j)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(0:nlat,nz)
    
    DOUBLE PRECISION :: v_on_w


        v_on_w=0.25d0*(v(i-1,j)+v(i,j)+v(i-1,j+1)+v(i,j+1))

    RETURN
    
  END FUNCTION  

  FUNCTION dv_dphi_on_u(v,i,j,dphi)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    
    DOUBLE PRECISION :: dv_dphi_on_u

        dv_dphi_on_u=(v(i,j)-v(i-1,j))/dphi

  
    RETURN
  END FUNCTION
  
  
  
  FUNCTION dv_dphi_on_v(v,i,j,dphi)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    DOUBLE PRECISION :: dv_dphi_on_v

    

        IF (i==0) THEN
          dv_dphi_on_v=0.5d0*(1.5d0*v(1,j)-1.5d0*v(0,j))/dphi
        ELSE IF (i==nlat) THEN
          dv_dphi_on_v=0.5d0*(1.5d0*v(nlat,j)-1.5d0*v(nlat-1,j))/dphi
        ELSE
          dv_dphi_on_v=0.5d0*(v(i+1,j)-v(i-1,j))/dphi
        END IF

  
    RETURN
  END FUNCTION
  
  FUNCTION dv_dr_on_v(v,i,j,dr)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(0:nlat,nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    DOUBLE PRECISION :: dv_dr_on_v

    

        IF (j==1) THEN
          dv_dr_on_v=0.5d0*(1.5d0*v(i,2)-1.5d0*v(i,1))/dr
        ELSE IF (j==nz) THEN
          dv_dr_on_v=0.5d0*(1.5d0*v(i,nz)-1.5d0*v(i,nz-1))/dr
        ELSE
          dv_dr_on_v=0.5d0*(v(i,j+1)-v(i,j-1))/dr
        END IF

  
    RETURN
  END FUNCTION
  
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Interpolation from w 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  FUNCTION w_on_u(w,i,j)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: w(nlat,0:nz)
    
    DOUBLE PRECISION :: w_on_u

    

        w_on_u=0.5d0*(w(i,j)+w(i,j-1))

  
    RETURN
  END FUNCTION
  
  FUNCTION w_on_v(w,i,j)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: w(nlat,0:nz)
    
    DOUBLE PRECISION :: w_on_v

    

        IF (i==0) THEN
          w_on_v=0.25d0*(2.5d0*w(1,j-1)-0.5d0*w(2,j-1)+2.5d0*w(1,j)-0.5d0*w(2,j))
        ELSE IF (i==nlat) THEN 
          w_on_v=0.25d0*(2.5d0*w(nlat,j-1)-0.5d0*w(nlat-1,j-1)+2.5d0*w(nlat,j)-0.5d0*w(nlat-1,j))
        ELSE
          w_on_v=0.25d0*(w(i,j-1)+w(i+1,j-1)+w(i,j)+w(i+1,j))
        END IF

  
    RETURN
  END FUNCTION
  
  FUNCTION dw_dphi_on_w(w,i,j,dphi)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: w(nlat,0:nz)
    DOUBLE PRECISION, INTENT(IN) :: dphi
    DOUBLE PRECISION :: dw_dphi_on_w


        IF (i==1) THEN
          dw_dphi_on_w=0.5d0*(1.5d0*w(2,j)-1.5d0*w(1,j))/dphi
        ELSE IF (i==nlat) THEN
          dw_dphi_on_w=0.5d0*(1.5d0*w(nlat,j)-1.5d0*w(nlat-1,j))/dphi
        ELSE
          dw_dphi_on_w=0.5d0*(w(i+1,j)-w(i-1,j))/dphi
        END IF

  
    RETURN
  END FUNCTION
  
  FUNCTION dw_dr_on_u(w,i,j,dr)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: w(nlat,0:nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    DOUBLE PRECISION :: dw_dr_on_u
    

        dw_dr_on_u=(w(i,j)-w(i,j-1))/dr

  
    RETURN
  END FUNCTION
  
  FUNCTION dw_dr_on_w(w,i,j,dr)
    
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: w(nlat,0:nz)
    DOUBLE PRECISION, INTENT(IN) :: dr
    DOUBLE PRECISION :: dw_dr_on_w
    

        dw_dr_on_w=0.5d0*(w(i,j+1)-w(i,j-1))/dr

  
    RETURN
  END FUNCTION
  
END MODULE    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
