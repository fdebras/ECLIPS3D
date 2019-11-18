PROGRAM eigenvectors_3D_sca

  IMPLICIT NONE

  INTEGER :: nlong
  INTEGER :: nlat
  INTEGER :: nz
  INTEGER :: ntot
  COMPLEX*16,ALLOCATABLE,  DIMENSION(:) :: p_temp
  COMPLEX*16,ALLOCATABLE, DIMENSION(:,:,:) :: p
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: pp
  COMPLEX*16 :: freq

  LOGICAL :: is_nm

  INTEGER :: zero_long,zero_z,zero_lat
  INTEGER, ALLOCATABLE :: zerolong(:),zeroz(:), zerolat(:)

  INTEGER :: i,j,k,a,eig

  CHARACTER(LEN=*), PARAMETER ::DIRDATA='../data/'
  OPEN(unit=1,file=DIRDATA // 'eig_sca.dat', &
    access='SEQUENTIAL')
  OPEN(unit=2, file=DIRDATA // 'frequency.dat', &
    access='SEQUENTIAL')
  OPEN(unit=10,file=DIRDATA // 'nummodes.dat', &
    access="SEQUENTIAL")

  OPEN(unit=4,file=DIRDATA // 'size.dat', &
    access="SEQUENTIAL")

  READ(4,*) nlong,nlat,nz,ntot
  CLOSE(4)
  ALLOCATE(p_temp(nlong*nlat*nz))
  ALLOCATE(p(nlong,nlat,nz))
  ALLOCATE(pp(nlong,nlat,nz))
  ALLOCATE(zerolong(0:nlong))
  ALLOCATE(zerolat(0:nlat))
  ALLOCATE(zeroz(0:nz))
  print *, ntot
    !----------------------------------------------------------
    ! Calculating the number of vertical and lateral zeros
    !----------------------------------------------------------
  READ(1,*)
    DO eig=1,ntot


      READ(1,9999) p_temp
      READ(2,9999) freq

      a=1
      DO i=1,nlong
        DO j=1,nlat
          DO k=1,nz
            p(i,j,k)=p_temp(a)
            a=a+1
          END DO
        END DO
      END DO


      pp=DBLE(p)+DIMAG(p)

      zerolong=0
      zerolat=0
      zeroz=0

      is_nm = verif(pp,nlong,nlat,nz)

      if (is_nm) THEN
      ! At each latitude-height, we calculate the number of zonal zeros.
      ! Then we fill an array which counts the number of latitude-height points
      ! having a certain number of zeros. The most common value
      ! will then be taken as the number of zonal zeros
      DO j=1,nlat-1
        DO k=1,nz
          zero_long=zeros(DBLE(p(:,j,k))+DIMAG(p(:,j,k)),nlong)
          zerolong(zero_long)=zerolong(zero_long)+1
        END DO
      END DO

      zero_long=0
      DO i=0,nlong
        IF (zerolong(i)>zerolong(zero_long)) THEN
          zero_long=i
        END IF
      END DO


      !Same for latitude
      DO i=1,nlong-1
        DO k=1,nz-1
         zero_lat=zeros(DBLE(p(i,:,k))+DIMAG(p(i,:,k)),nlat)
         zerolat(zero_lat)=zerolat(zero_lat)+1
        END DO
      END DO

      zero_lat=0
      DO j=0,nlat
        IF (zerolat(j)>zerolat(zero_lat)) THEN
          zero_lat=j
        END IF
      END DO


      ! Same for height
      DO i=1,nlong-1
        DO j=1,nlat
          zero_z=zeros(DBLE(p(i,j,:))+DIMAG(p(i,j,:)),nz-1)
          zeroz(zero_z)=zeroz(zero_z)+1
        END DO
      END DO

      zero_z=0
      DO k=0,nz
        IF (zeroz(k)>zeroz(zero_z)) THEN
          zero_z=k
        END IF
      END DO

      ! Same for longitude




      IF  (zero_long<=3 .and. zero_lat<=3 .and. zero_z<=3) THEN

        print *, "n=",eig , "freq= ", DBLE(freq), 'im_freq=', DIMAG(freq), &
            'kz=',  zero_z, 'klat=',zero_lat, 'klong=', zero_long
       write(10,*) eig, freq
!         IF ((eig==3361) .or. (eig==3362) .or. (eig==7560) .or. (eig==7561) .or. (eig==9194)) THEN
!           WRITE(44,*) zero_lat,zero_z
!           WRITE(44,*) freq_r(k),freq_i(k)
!           WRITE(44,888) REALPART(u)
!           WRITE(44,888) IMAGPART(v)
!           WRITE(44,888) REALPART(p)
!           WRITE(44,888) REALPART(theta)
!           WRITE(44,888) IMAGPART(w)
         END IF

      END IF
      !print *, zero_lat


    END DO
    DEALLOCATE(zeroz)
    DEALLOCATE(zerolat)
    DEALLOCATE(zerolong)
    DEALLOCATE(p)
    DEALLOCATE(pp)
    DEALLOCATE(p_temp)
    WRITE(10,*), 0, (0.d0, 0.d0)
    CLOSE(10)
9999 FORMAT( E15.8,E15.8 )
    CONTAINS


     FUNCTION zeros(func,n)
    !----------------------------------------------------------
    !----------------------------------------------------------
    ! Function to calculate the number of zeros in a 1D array
    !----------------------------------------------------------
    !----------------------------------------------------------
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: func(n)

    INTEGER :: zeros

    DOUBLE PRECISION :: eps ! We just need a constraint to avoid constant zeros functions
    ! and negligible oscillations around a zero
    INTEGER :: i ! for loop
    DOUBLE PRECISION :: plusminus ! REAL equal to +-1 to see if there is a change in sign
    DOUBLE PRECISION :: ONE=1.D0
    zeros=0

    ! Initialisation : choice of epsilon
    eps=ABS(func(1))
    DO i=2,n
      eps=MAX(eps,ABS(func(i)))
    END DO
    eps=0.0001*eps

    ! Then we take plusminus equal to the opposite of the sign of the first
    ! terms of the function
    i=1
    DO WHILE (DABS(func(i))<eps .AND. i<=n)
      i=i+1
    END DO
    plusminus=-SIGN(ONE,func(i))

    ! If the function changes sign, there is a zero ! And plusminus changes
    ! sign again
    DO WHILE (i<=n)
      IF (func(i)*plusminus>eps) THEN ! positive only if func(i) changed sign
        zeros=zeros+1
        plusminus=-plusminus
      END IF
      i=i+1
    END DO


    RETURN
  END FUNCTION

FUNCTION verif(func,nx,ny,nz)

    INTEGER, INTENT(IN) :: nx,ny,nz
    DOUBLE PRECISION, INTENT(IN) :: func(nx,ny,nz)

    INTEGER :: i,j,k,t,z

    LOGICAL :: verif


    verif=.FALSE.
    t=0
    DO i=1,nx
      DO j=1,ny
        DO k=1,nz
          IF (DABS(func(i,j,k))>1d-7) THEN
            t=t+1
          END IF
        END DO
      END DO
    END DO

    IF (t>(nx*ny*nz/3.d0)) THEN
      verif=.TRUE.
    END IF

    z=0

    DO k=1,nz
    DO i=1,nx
      t=0
      DO j=2,ny-1
        IF (func(i,j,k)*func(i,j-1,k)<0.d0) THEN
          t=t+1
        END IF
      END DO
      IF (t>3) THEN
        z=z+1
      END IF
    END DO
    END DO

    IF (z>5) THEN
      verif=.FALSE.
    END IF
    RETURN
  END FUNCTION

END PROGRAM
