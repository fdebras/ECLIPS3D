PROGRAM find_eigenvectors
  
  
  IMPLICIT NONE
  
  
  INTEGER, PARAMETER :: nlong=40
  INTEGER, PARAMETER :: nlat=20
  INTEGER, PARAMETER :: ntot=2*nlong*nlat+nlong*(nlat+1)
  COMPLEX*16, DIMENSION(nlong*nlat) :: p_temp 
  COMPLEX*16, DIMENSION(nlong,nlat) :: p
  DOUBLE PRECISION, DIMENSION(nlong,nlat) :: pp
  COMPLEX*16 :: freq
  
  COMPLEX*16 :: test
  
  LOGICAL :: is_nm
  
  INTEGER :: zero_long,zero_lat
  INTEGER :: zerolong(0:nlong), zerolat(0:nlat)
  
  CHARACTER(LEN=300) :: DIR

  INTEGER :: i,j,a,eig

  DIR = '/Users/florian/Desktop/test/ECLIPS3D/trunk/2D_shallow/data/'
   
  OPEN(unit=1,file=TRIM(DIR) // 'eig_sca.dat', &
    access='SEQUENTIAL')
  OPEN(unit=2, file=TRIM(DIR) // 'fort.43', & 
    access='SEQUENTIAL')
  OPEN(unit=5,file=TRIM(DIR) // 'nummodes.dat', &
    access="SEQUENTIAL")
  READ(1,*)
  

  DO eig=1,ntot

    READ(1,9999) p_temp
    READ(2,9999) freq
  
    a=1
    DO i=1,nlong
      DO j=1,nlat
        p(i,j)=p_temp(a)
        a=a+1
      END DO
    END DO

    pp=DBLE(p)!+DIMAG(p)
  
    is_nm=verif(pp,nlong,nlat)
    IF(is_nm)  THEN
      zerolong=0
      zerolat=0
	! At each latitude, we calculate the number of vertical zeros.
	! Then we fill an array which counts the number of latitudinal points 
	! having a certain number of zeros. The most common value 
	! will then be taken as the number of vertical zeros
      DO i=1,nlat
        zero_long=zeros(pp(:,i),nlong)
        zerolong(zero_long)=zerolong(zero_long)+1
      END DO
  
      zero_long=0
      DO j=0,nlong
        IF (zerolong(j)>zerolong(zero_long)) THEN
          zero_long=j
        END IF
      END DO 
  
	! Same for longitude
      DO j=1,nlong
        zero_lat=zeros(pp(j,:),nlat)
        zerolat(zero_lat)=zerolat(zero_lat)+1
      END DO
  
      zero_lat=0
      DO i=0,nlat
        IF (zerolat(i)>zerolat(zero_lat)) THEN
          zero_lat=i
        END IF
      END DO 
!       print *, MAXVAL(pp)-MINVAL(pp)
!      print *,'k= ', eig, 'freq =', freq
      IF ((zero_long<=4 .AND. zero_lat<=3)) THEN !.OR. (DIMAG(freq)<-1.d-5)) THEN  
        print *,'k= ', eig, 'freq =', freq
        WRITE(5,*) eig, freq
      END IF
    END IF
    
  END DO
  WRITE(5,*) 0, (0.d0 , 0.d0)
  
  
  CLOSE(5)
  CLOSE(1)
  CLOSE(2)
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
    DOUBLE PRECISION :: plusminus ! real equal to +-1 to see if there is a change in sign
    zeros=0
  
    ! Initialisation : choice of epsilon
    eps=DABS(func(1))
    DO i=2,n
      eps=MAX(eps,DABS(func(i)))
    END DO
    eps=0.0001d0*eps
    
    ! Then we take plusminus equal to the opposite of the sign of the first 
    ! terms of the function
    i=1
    DO WHILE (DABS(func(i))<eps .AND. i<=n)
      i=i+1
    END DO
    plusminus=-SIGN(1.d0,func(i))
    
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

FUNCTION verif(func,nx,ny)
  
    INTEGER, INTENT(IN) :: nx,ny
    DOUBLE PRECISION, INTENT(IN) :: func(nx,ny)
    
    INTEGER :: i,j,k,t

    LOGICAL :: verif
    
    
    verif=.FALSE.
    k=0
    DO i=1,nx
      DO j=1,ny
      
        IF (DABS(func(i,j))>1.d-12) THEN
          k=k+1
        END IF
        
      END DO
    END DO
    
    IF (k>(nx*ny/2.d0)) THEN
      verif=.TRUE.
    END IF
    
    k=0
    
    DO i=1,nx
      t=0
      DO j=1,ny-1
        IF (func(i,j)*func(i,j-1)<0.d0) THEN
          t=t+1
        END IF
      END DO
      IF (t>4) THEN
        k=k+1
      END IF
    END DO 
    
    IF (k>10) THEN
      verif=.FALSE.
    END IF
    RETURN
  END FUNCTION

  
END PROGRAM
