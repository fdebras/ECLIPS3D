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
  SUBROUTINE fill_matrix(nu_data, nv_data, nh_data)

    IMPLICIT NONE
    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Subroutine arguments
    !----------------------------------------------------------
    !----------------------------------------------------------

    !Size of input arrays
    INTEGER, INTENT(IN) :: nu_data
    INTEGER, INTENT(IN) :: nv_data
    INTEGER, INTENT(IN) :: nh_data

    !----------------------------------------------------------
    !----------------------------------------------------------
    !Local Variables
    !----------------------------------------------------------
    !----------------------------------------------------------


    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: u_array! Steady state
    ! functions on u grid

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: v_array! Steady state
    ! functions on v grid

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: h_array! Steady state
    ! functions on h grid

    !----------------------------------------------------------
    ! Steady state variables (see grid)
    !----------------------------------------------------------

    ! Steady on u
    DOUBLE PRECISION, DIMENSION(nlong,nlat) :: Us_u, dUs_dx_u, dUs_dy_u, &
    Vs_u

    ! Steady on v
    DOUBLE PRECISION, DIMENSION(nlong,0:nlat) :: Us_v, Vs_v, &
    dVs_dx_v, dVs_dy_v

    !Steady on h
    DOUBLE PRECISION, DIMENSION(nlong,nlat) :: Hs_h, dHs_dx_h,dHs_dy_h, &
    Us_h,dUs_dx_h, Vs_h, dVs_dy_h

    !----------------------------------------------------------
    ! Perturbed variables (see grid)
    !----------------------------------------------------------

    DOUBLE PRECISION :: u(nlong,nlat), v(nlong,0:nlat), h(nlong,nlat)

    !----------------------------------------------------------
    ! Useful quantities
    !----------------------------------------------------------

    DOUBLE PRECISION :: beta

    INTEGER :: i,j,k,pos! Loops integers

    DOUBLE PRECISION :: x_u(nlong),x_v(nlong),y_u(nlat), y_v(0:nlat)
    ! latitude full and half
    DOUBLE PRECISION :: dx,dy ! height and latitude steps
    DOUBLE PRECISION :: S(ntot) ! Just used for computing the matrix
    COMPLEX*16 :: alpha

    DOUBLE PRECISION :: L_nd !, T_nd, V_nd ! non dimensionalising quantities
    ! only need L in here
    COMPLEX*16 :: ic = (0.d0,1.d0)

    DOUBLE PRECISION :: mat_test(ntot,ntot)

    DOUBLE PRECISION :: ifreq(ntot), rfreq(ntot), test_Vr(ntot,ntot)

    INTEGER :: LWORK
    DOUBLE PRECISION,ALLOCATABLE :: WORK(:)
    INTEGER :: info_test

    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Subroutine
    !----------------------------------------------------------
    !----------------------------------------------------------


    LWORK = 10*ntot
    ALLOCATE(WORK(LWORK))

    open(unit=1,file=TRIM(DIRDATA) // "u_data.dat", access= &
    'sequential')

    open(unit=2,file=TRIM(DIRDATA) // "v_data.dat", access= &
    'sequential')

    open(unit=3,file=TRIM(DIRDATA) // "h_data.dat", access= &
    'sequential')

    ALLOCATE(u_array(nu_data))
    ALLOCATE(v_array(nv_data))
    ALLOCATE(h_array(nh_data))

    READ(1,*) u_array
    READ(2,*) v_array
    READ(3,*) h_array

    CLOSE(1)
    CLOSE(2)
    CLOSE(3)
    !----------------------------------------------------------
    ! Initialisation of the variables
    !----------------------------------------------------------

    !-----------------------------------------------------------
    ! Reading the input array according to grid
    k=1

    !------------
    !U
    DO i=1,nlong
      DO j=1,nlat
        Us_u(i,j)=u_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        dUs_dx_u(i,j)=u_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        dUs_dy_u(i,j)=u_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        Vs_u(i,j)=u_array(k)
        k=k+1
      END DO
    END DO

    !------------
    !V
    k=1
    DO i=1,nlong
      DO j=0,nlat
        Us_v(i,j)=v_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        Vs_v(i,j)=v_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        dVs_dx_v(i,j)=v_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=0,nlat
        dVs_dy_v(i,j)=v_array(k)
        k=k+1
      END DO
    END DO


    !------------
    !H
    k=1
    DO i=1,nlong
      DO j=1,nlat
        Us_h(i,j)=h_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        dUs_dx_h(i,j)=h_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        Vs_h(i,j)=h_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        dVs_dy_h(i,j)=h_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        Hs_h(i,j)=h_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        dHs_dx_h(i,j)=h_array(k)
        k=k+1
      END DO
    END DO

    DO i=1,nlong
      DO j=1,nlat
        dHs_dy_h(i,j)=h_array(k)
        k=k+1
      END DO
    END DO

    IF(myrow == 0 .AND. mycol == 0) print *, 'all read'

!    print *, MAXVAL(Hs_h)
!    print *, MAXVAL(dHs_dx_h)

!    print*, MINVAL(Hs_h)
!--------------------------------------------------------
! Now initialisation
!---------------------------------------------------------


    beta = 2.d0*omega*DCOS(phi)/rtot

    L_nd = DSQRT(DSQRT(g*height)/beta)

    !Defining height and latitude steps
    dx = (2.d0*pi*rtot/L_nd)/DBLE(nlong) ! x goes from east to west
    dy= (pi*rtot/(1.3d0*L_nd))/DBLE(nlat) ! y is a hemisphere wide

    ! Defining height and gravity
    DO i=1,nlong
      x_u(i) = (DBLE(i)-0.5d0)*dx
      x_v(i) = dble(i)*dx
    END DO

    DO j=1,nlat
      y_u(j) = (dble(j)-0.5d0)*dy
      y_v(j) = dble(j)*dy
    END DO

    y_v(0) = 0.d0

    !----------------------------------------------------------
    ! Seek for eigenvectors
    !----------------------------------------------------------

    ! To create the matrix, we will just set all the variables to
    ! 0 at each point of the grid and allow one to be non zero at one point
    ! of the grid. Then we will calculate the evolution coefficients from the
    ! perturbed equations due to this value and store it in the matrix from
    ! which we will get the eigenvectors and values.


    DO k=1,ntot
      S(k)=0.d0
    END DO

    DO pos=1,ntot
      S(pos)=1.d0

      u=0.d0
      v=0.d0
      h = 0.d0

      k=1
      DO i=1,nlong
        DO j=1,nlat
          u(i,j)=S(k)
          k=k+1
        END DO
      END DO

      DO i=1,nlong
        DO j=0,nlat
          v(i,j)=S(k)
          k=k+1
        END DO
      END DO

      DO i=1,nlong
        DO j=1,nlat
          h(i,j)=S(k)
          k=k+1
        END DO
      END DO

       DO i=1,nlong
         v(i,0) = 0.d0
    !     v(i,nlat) = 0.d0
       END DO
      k=1

      DO i=1,nlong
        DO j=1, nlat
          alpha = DCMPLX(0.d0,0.71710110656588177*h(i,j))+ &
          DCMPLX(Us_u(i,j)*du_dx_on_u(u,i,j,dx) + u(i,j)*dUs_dx_u(i,j) + &
          Vs_u(i,j)*du_dy_on_u(u,i,j,dy) + v_on_u(v,i,j)*dUs_dy_u(i,j) - &
          y_u(j)*v_on_u(v,i,j)+ &
          tdrag_invert*u(i,j),0.d0)

          ! In the matrix, we put -1/i * alpha so that the eigenvalues
          ! are exactly the frequencies
          CALL PZELSET(mat_evol,k,pos,desc_mat,alpha)
          mat_test(pos,k) = alpha
          k=k+1

        END DO
      END DO

      DO i=1,nlong
        DO j=0, nlat

          IF (j==0) THEN
!             alpha = Us_v(i,j)*dv_dx_on_v(v,i,j,dx)+u_on_v(u,i,j)*dVs_dx_v(i,j) + &
!             Vs_v(i,j)*dv_dy_on_v(v,i,j,dy)+v(i,j)*dVs_dy_v(i,j) + &
!             dh_dy_on_v(h,i,j,dy)+y_v(j)*u_on_v(u,i,j)+ &
!             tdrag_invert*v(i,j)

            alpha = DCMPLX(0.d0,0.d0)
          !ELSE IF (j==nlat) THEN
          !  alpha = 0.d0
          ELSE
            alpha = DCMPLX(Us_v(i,j)*dv_dx_on_v(v,i,j,dx)+u_on_v(u,i,j)*dVs_dx_v(i,j) + &
            Vs_v(i,j)*dv_dy_on_v(v,i,j,dy)+v(i,j)*dVs_dy_v(i,j) + &
            dh_dy_on_v(h,i,j,dy)+y_v(j)*u_on_v(u,i,j)+ &
            tdrag_invert*v(i,j),0.d0)
          END IF
          CALL PZELSET(mat_evol,k,pos,desc_mat,alpha)
          mat_test(pos,k) = alpha
          k=k+1
        END DO
      END DO

      DO i=1,nlong
        DO j=1, nlat
          alpha = DCMPLX(0.d0,Hs_h(i,j)*0.71710110656588177*u(i,j))+ &
          DCMPLX(u_on_h(u,i,j)*dHs_dx_h(i,j)+ &
          Hs_h(i,j)*dv_dy_on_h(v,i,j,dy)+v_on_h(v,i,j)*dHs_dy_h(i,j)+ &
          h(i,j)*dUs_dx_h(i,j) + Us_h(i,j)*dh_dx_on_h(h,i,j,dx)+ &
          h(i,j)*dVs_dy_h(i,j) + Vs_h(i,j)*dh_dy_on_h(h,i,j,dy) + &
          trad_invert*h(i,j),0.d0)

          ! In the matrix, we put -1/i * alpha so that the eigenvalues
          ! are exactly the frequencies
          CALL PZELSET(mat_evol,k,pos,desc_mat,alpha)
          mat_test(pos,k) = alpha
          k=k+1
        END DO
      END DO

    S(pos)=0.d0

    END DO
!        print *, MAXVAL(mat_test)
!    print *, MINVAL(mat_test)


!     CALL DGEEV('N','V',ntot,mat_test,ntot,rfreq,ifreq,test_Vr,ntot,&
!     test_Vr,ntot,WORK,LWORK,info_test)
!
!     print *, (info_test==0)
!     open(unit=11,file=TRIM(DIRDATA) // "test_freq.dat", access= &
!     'sequential')
!     open(unit=12,file=TRIM(DIRDATA) // "test_ifreq.dat", access= &
!     'sequential')
!     WRITE(11,*) rfreq
!     WRITE(12,*) ifreq
!
!     CLOSE(11)
!     CLOSE(12)

  888 FORMAT(10E14.6 )

  END SUBROUTINE


!   	DO i=1,nlat
!   	  DO j=1,nz-1
!   	    CALL PZELSET(mat_evol,k,pos,desc_mat,theta_t(i,j))
!   	    k=k+1
!   	  END DO
!   	END DO


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Interpolation from u
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  FUNCTION u_on_v(u,i,j)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat)

    DOUBLE PRECISION :: u_on_v
    IF (i==nlong) THEN
      IF (j==0) THEN
        u_on_v = 0.25d0*(2.5d0*u(nlong,1)-0.5d0*u(nlong,2)+ &
        2.5d0*u(1,1)-0.5d0*u(1,2))

      ELSE IF (j==nlat) THEN
        u_on_v= 0.25d0*(2.5d0*u(nlong,nlat)-0.5d0*u(nlong,nlat-1)+&
        2.5d0*(1.5d0*u(nlong,nlat)-0.5d0*u(nlong-1,nlat))- &
        0.5d0*(1.5d0*u(nlong,nlat-1)-0.5d0*u(nlong-1,nlat-1)))

      ELSE
        u_on_v = 0.25d0*(u(nlong,j)+1.5d0*u(nlong,j)-0.5d0*u(nlong-1,j)+ &
        u(nlong,j+1)+ 1.5d0*u(nlong,j+1)-0.5d0*u(nlong-1,j))
      END IF

    ELSE
      IF (j==0) THEN
        u_on_v = 0.25d0*(2.5d0*u(i,1)-0.5d0*u(i,2)+ &
        2.5d0*u(i+1,1)-0.5d0*u(i+1,2))

      ELSE IF (j==nlat) THEN
        u_on_v = 0.25d0*(2.5d0*u(i,nlat)-0.5d0*u(i,nlat-1)+&
        2.5d0*u(i+1,nlat)-0.5d0*u(i+1,nlat-1))

      ELSE
        u_on_v = 0.25d0*(u(i,j)+u(i+1,j)+u(i,j+1)+ u(i+1,j+1))
      END IF
    END IF

    RETURN

  END FUNCTION


  FUNCTION u_on_h(u,i,j)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat)

    DOUBLE PRECISION :: u_on_h


    IF (i==nlong) THEN
      u_on_h = 0.5d0*(u(nlong,j)+1.5d0*u(nlong,j)-0.5d0*u(nlong-1,j))
    ELSE
      u_on_h = 0.5d0*(u(i,j)+u(i+1,j))
    END IF

    RETURN

  END FUNCTION

  FUNCTION du_dx_on_u(u,i,j,dx)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat)
    DOUBLE PRECISION, INTENT(IN) :: dx

    DOUBLE PRECISION :: du_dx_on_u


    IF (i==1) THEN
      du_dx_on_u = 0.5d0*(u(2,j)-(1.5d0*u(1,j)-0.5d0*u(2,j)))/dx
    ELSE IF (i==nlong) THEN
      du_dx_on_u = 0.5d0*((1.5d0*u(nlong,j)-0.5d0*u(nlong-1,j))-u(nlong-1,j))/dx
    ELSE
      du_dx_on_u = 0.5d0*(u(i+1,j)-u(i-1,j))/dx
    END IF

    RETURN

  END FUNCTION

  FUNCTION du_dy_on_u(u,i,j,dy)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat)
    DOUBLE PRECISION, INTENT(IN) :: dy

    DOUBLE PRECISION :: du_dy_on_u


    IF (j==1) THEN
      du_dy_on_u = 0.5d0*(1.5d0*u(i,2)-1.5d0*u(i,1))/dy
    ELSE IF (j==nlat) THEN
      du_dy_on_u = 0.5d0*(1.5d0*u(i,nlat)-1.5d0*u(i,nlat-1))/dy
    ELSE
      du_dy_on_u = 0.5d0*(u(i,j+1)-u(i,j-1))/dy
    END IF

    RETURN

  END FUNCTION


  FUNCTION du_dx_on_h(u,i,j,dx)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: u(nlong,nlat)
    DOUBLE PRECISION, INTENT(IN) :: dx

    DOUBLE PRECISION :: du_dx_on_h

    IF (i==nlong) THEN
      du_dx_on_h = (1.5d0*u(nlong,j)-0.5d0*u(nlong-1,j)-u(nlong,j))/dx
    ELSE
      du_dx_on_h = (u(i+1,j)-u(i,j))/dx
    END IF

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
    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat)

    DOUBLE PRECISION :: v_on_u


    IF (i==1) THEN
      v_on_u = 0.25d0*(v(1,j)+1.5d0*v(1,j)-0.5d0*v(2,j)+ &
      v(1,j-1)+1.5d0*v(1,j-1)-0.5d0*v(2,j-1))
    ELSE
      v_on_u = 0.25d0*(v(i,j)+v(i-1,j)+v(i,j-1)+v(i-1,j-1))
    END IF

    RETURN

  END FUNCTION

  FUNCTION v_on_h(v,i,j)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat)

    DOUBLE PRECISION :: v_on_h

    v_on_h = 0.5d0*(v(i,j)+v(i,j-1))

    RETURN

  END FUNCTION

  FUNCTION dv_dx_on_v(v,i,j,dx)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat)
    DOUBLE PRECISION, INTENT(IN) :: dx

    DOUBLE PRECISION :: dv_dx_on_v


    IF (i==1) THEN
      dv_dx_on_v = 0.5d0*(v(2,j)-v(nlong,j))/dx
    ELSE IF (i==nlong) THEN
      dv_dx_on_v = 0.5d0*(v(1,j)-v(nlong-1,j))/dx
    ELSE
      dv_dx_on_v = 0.5d0*(v(i+1,j)-v(i-1,j))/dx
    END IF

    RETURN

  END FUNCTION

  FUNCTION dv_dy_on_v(v,i,j,dy)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat)
    DOUBLE PRECISION, INTENT(IN) :: dy

    DOUBLE PRECISION :: dv_dy_on_v


    IF (j==0) THEN
      dv_dy_on_v = 0.5d0*(1.5d0*v(i,1)-1.5d0*v(i,0))/dy
    ELSE IF (j==nlat) THEN
      dv_dy_on_v = 0.5d0*(1.5d0*v(i,nlat)-1.5d0*v(i,nlat-1))/dy
    ELSE
      dv_dy_on_v = 0.5d0*(v(i,j+1)-v(i,j-1))/dy
    END IF

    RETURN

  END FUNCTION


  FUNCTION dv_dy_on_h(v,i,j,dy)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: v(nlong,0:nlat)
    DOUBLE PRECISION, INTENT(IN) :: dy

    DOUBLE PRECISION :: dv_dy_on_h

    dv_dy_on_h = (v(i,j)-v(i,j-1))/dy

    RETURN
  END FUNCTION


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Interpolation from h
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  FUNCTION dh_dx_on_u(h,i,j,dx)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: h(nlong,nlat)
    DOUBLE PRECISION, INTENT(IN) :: dx

    DOUBLE PRECISION :: dh_dx_on_u

    IF (i==1) THEN
      dh_dx_on_u = (h(1,j)-(1.5d0*h(1,j)-0.5d0*h(2,j)))/dx
    ELSE
      dh_dx_on_u = (h(i,j)-h(i-1,j))/dx
    END IF

    RETURN
  END FUNCTION

  FUNCTION dh_dy_on_v(h,i,j,dy)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: h(nlong,nlat)
    DOUBLE PRECISION, INTENT(IN) :: dy

    DOUBLE PRECISION :: dh_dy_on_v

    IF (j==0) THEN
      dh_dy_on_v = (0.5d0*h(i,2)-0.5d0*h(i,1))/dy
    ELSE IF (j==nlat) THEN
      dh_dy_on_v = (0.5d0*h(i,nlat)-0.5d0*h(i,nlat-1))/dy
    ELSE
      dh_dy_on_v = (h(i,j+1)-h(i,j))/dy
    END IF

    RETURN
  END FUNCTION

  FUNCTION dh_dx_on_h(h,i,j,dx)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: h(nlong,nlat)
    DOUBLE PRECISION, INTENT(IN) :: dx

    DOUBLE PRECISION :: dh_dx_on_h

    IF (i==1) THEN
      dh_dx_on_h = 0.5d0*(h(2,j)-h(nlong,j))/dx
    ELSE IF (i==nlong) THEN
      dh_dx_on_h = 0.5d0*(h(1,j)-h(nlong-1,j))/dx
    ELSE
      dh_dx_on_h = 0.5d0*(h(i+1,j)-h(i-1,j))/dx
    END IF

    RETURN
  END FUNCTION

  FUNCTION dh_dy_on_h(h,i,j,dy)
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, INTENT(IN) :: h(nlong,nlat)
    DOUBLE PRECISION, INTENT(IN) :: dy

    DOUBLE PRECISION :: dh_dy_on_h


    IF (j==1) THEN
      dh_dy_on_h = 0.5d0*(1.5d0*h(i,2)-1.5d0*h(i,1))/dy
    ELSE IF (j==nlat) THEN
      dh_dy_on_h = 0.5d0*(1.5d0*h(i,nlat)-1.5d0*h(i,nlat-1))/dy
    ELSE
      dh_dy_on_h = 0.5d0*(h(i,j+1)-h(i,j-1))/dy
    END IF

    RETURN

  END FUNCTION


END MODULE
