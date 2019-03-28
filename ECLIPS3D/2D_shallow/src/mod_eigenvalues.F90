MODULE mod_eigenvalues

  USE mod_data
  USE mod_init_para
  USE mod_init_matrix
  USE mod_fill_matrix
  IMPLICIT NONE
  INTEGER, DIMENSION(50) :: descVR, descVL
  COMPLEX*16, DIMENSION(:), ALLOCATABLE:: freq
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE:: vl, vr
  
  
CONTAINS 

  SUBROUTINE evalues
  
  INTEGER, DIMENSION(9) :: descTemp, descW
  INTEGER*8 :: LWORK
  
  INTEGER :: irow, icol, proc_row,proc_col
  
  INTEGER  :: liwork=1000
  INTEGER :: i,j,k,ki, mput
  LOGICAL, DIMENSION(ntot) :: sel_dum  
  INTEGER :: IA,JA
  
  INTEGER :: info
  
  INTEGER :: ILO, IHI
  INTEGER, DIMENSION(:), ALLOCATABLE :: scale
  
  COMPLEX*16 :: ONE=(1.d0,0.d0)
  COMPLEX*16 :: ZERO=(0.d0,0.d0)
  
  COMPLEX*16, DIMENSION(ntot-1):: tau
  COMPLEX*16 :: alpha, shift
  
  COMPLEX*16, DIMENSION(:), ALLOCATABLE  :: WORK, IWORK

  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: vr_temp
  
  DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: RWORK  
  EXTERNAL BLACS_GRIDINFO, DESCINIT, PZELSET, PDHSEQR, PDGEBAL, &
  PDGEHRD, BLACS_GRIDEXIT, BLACS_EXIT

  !CALL INIT_M
  


  ALLOCATE(IWORK(LIWORK))
  print *, nlld
  
  CALL DESCINIT(descVR, ntot,ntot,nb,nb,0,0,info_txt,nlld,info)
  CALL DESCINIT(descVL, ntot,ntot,nb,nb,0,0,info_txt,nlld,info)

  ALLOCATE(vr(nlld,nlld))
  ALLOCATE(vl(nlld,nlld))
  ALLOCATE(freq(ntot))
  
  IF (myrow==0 .and. mycol==0) THEN
    print *, 'Begin eigenvalues finder !'
  END IF
  
  CALL PZLASET(' ', ntot,ntot,ZERO, ONE, vr, 1,1, descVR)
  
  ILO=1
  IHI=ntot

  IA=1
  JA=1

!  CALL PZLAPRNT( N, N, LM, 1, 1, DESCLM, 0, 0, 'A', 6, WORK )
  CALL PZGEHRD(ntot,ILO,IHI,mat_evol,IA,JA,desc_mat,tau,IWORK,-1,info)
  LWORK=IWORK(1)
  
  ALLOCATE(WORK(LWORK))
  CALL PZGEHRD(ntot,ILO,IHI,mat_evol,IA,JA,desc_mat,tau,WORK,LWORK,info)
  DEALLOCATE(WORK)
  IF (info==0) THEN

  ELSE 
    print *, "PZGEHRD FAILURE"
  END IF
  
  !CALL PDORMHR('R','N',n,n,ILO,IHI,M,1,1,descM,tau,vr,1,1,descVR,work,lwork,info)
  CALL PZUNMHR('R','N',ntot,ntot,ILO,IHI,mat_evol,IA,JA,desc_mat,tau,VR,IA,JA,descVR,iwork,-1,info)
  
  LWORK=IWORK(1)
  ALLOCATE(WORK(LWORK))
  CALL PZUNMHR('R','N',ntot,ntot,ILO,IHI,mat_evol,IA,JA,desc_mat,tau,VR,IA,JA,descVR,work,LWORK,info)
  DEALLOCATE(WORK)
  !PRINT *, REALPART(VR)
  IF (info==0) THEN

  ELSE 
    print *, "PZUNMHR FAILURE"
  END IF
  !print *, REALPART(LM)
  
   DO i=1,ntot
     DO j=1,ntot
       IF (j<i-1) THEN
         CALL PZELSET(mat_evol,i,j,desc_mat,ZERO)
       END IF
     END DO
   END DO
   

  IF (myrow==0 .and. mycol==0) THEN
    print *, 'Diagonal killed'
  END IF
    
  call PZLAHQR(.TRUE.,.TRUE.,ntot,ILO,IHI,mat_evol,desc_mat,freq,1,ntot,VR,descVR, &
  iwork,-1,iwork,LIWORK,info)
  LWORK=IWORK(1)
  
  print *, LWORK  
  
  ALLOCATE(WORK(LWORK))
  call PZLAHQR(.TRUE.,.TRUE.,ntot,ILO,IHI,mat_evol,desc_mat,freq,1,ntot,VR,descVR, &
  work,LWORK,iwork,LIWORK,info)
  
  IF (info==0) THEN
  ELSE 
    print *, "PZLAHQR FAILURE"
  END IF

  IF (myrow==0 .and. mycol==0) THEN
    print *, 'Last routine'
  END IF
  
  
    ALLOCATE(RWORK(LWORK))
  CALL PZTREVC('R','B',sel_dum,ntot,mat_evol,desc_mat,VL,descVL,VR,descVR,ntot,mput,work,rwork,info)
  IF (info==0) THEN  
  ELSE 
    print *, "PZTREVC FAILURE"
  END IF


!   DEALLOCATE(WORK)
!   
!   CALL DESCINIT(descTemp,ntot,ntot,nb,nb,0,0,info_txt,nlld,info)
!   CALL DESCINIT(descW, ntot,1,nb,1,0,0,info_txt,nlld,info)
!   
!   ALLOCATE(WORK(ntot))
!   ALLOCATE(vr_temp(nlld,nlld))
! 
!   CALL PZLASET(' ', ntot,ntot, ZERO, ZERO, vr_temp, 1,1,descTemp)
!   
!   
!   
!     
!   ! Store the diagonal elements of mat_evol
!   
!   DO i=1,ntot
!     CALL INFOG2L( i, i, desc_mat, nprow, npcol, myrow, mycol, irow, &
!                       icol, proc_row, proc_col)
!     IF ((myrow.eq.proc_row) .AND. (mycol.eq.proc_col)) THEN
!       WORK(i)=mat_evol(irow,icol)
!     END IF
!   END DO  
!   
! 
! 
!   ! Perform the operation (A-lambda*I + ( 0 0 0 0 // 000// ...//0 0 1 0000 ))x = (00 ..1.. 00 )
!   ! at each matrix size.
!   
!   DO k=ntot,1,-1
!     SHIFT = ZERO
!     CALL INFOG2L( k, k, desc_mat, nprow, npcol, myrow, mycol, irow, &
!                       icol, proc_row, proc_col)
!     IF ((myrow.eq.proc_row) .AND. (mycol.eq.proc_col)) THEN
!       shift=mat_evol(irow,icol)
!     END IF
!     
!     ! Store shift in all processors
!     CALL ZGSUM2D( info_txt, 'ALL', ' ', 1, 1, shift, 1, -1, -1 )
!     DO ki=1,k
!       ! Shift the diagonal elements
!       CALL INFOG2L( ki, ki, desc_mat, nprow, npcol, myrow, mycol, irow, &
!                       icol, proc_row, proc_col)
!       IF ((myrow.eq.proc_row) .AND. (mycol.eq.proc_col)) THEN
!         mat_evol(irow,icol)=mat_evol(irow,icol)-shift
!       END IF 
!     END DO
!     
!     ! Keep one on the last diagonal element
!     CALL INFOG2L( k, k, desc_mat, nprow, npcol, myrow, mycol, irow, &
!                       icol, proc_row, proc_col)
!     IF ((myrow.eq.proc_row) .AND. (mycol.eq.proc_col)) THEN
!       mat_evol(irow,icol)=ONE                 
!     END IF 
!     
!     ! Put One on the diagonal element of vr_temp
!     CALL INFOG2L( k, k, descTemp, nprow, npcol, myrow, mycol, irow, &
!                       icol, proc_row, proc_col)
!     IF ((myrow.eq.proc_row) .AND. (mycol.eq.proc_col)) THEN
!       vr_temp(irow,icol)=ONE                 
!     END IF 
!     
!     ! Solve the equation on the kth column of vr_temp
!     CALL PZTRSV ('U', 'N','N',k,mat_evol,1,1,desc_mat,vr_temp,1,k,descTemp,1)
!     
!     
!     ! Put back the good diagonal elements of T 
!     
!     DO ki=1,k
!       CALL INFOG2L( ki, ki, desc_mat, nprow, npcol, myrow, mycol, irow, &
!                       icol, proc_row, proc_col)
!       IF ((myrow.eq.proc_row) .AND. (mycol.eq.proc_col)) THEN
!         mat_evol(irow,icol)=work(ki)
!       END IF 
!       
!       
!     END DO
!     CALL PZGEMV('N', ntot,k,ONE,vr,1,1,descVR,vr_temp,1,k,descTemp,1,ZERO, &
!     vr,1,k, descVR,1)
!     
!   END DO
    
  IF (myrow==0 .and. mycol==0) THEN
    print *, 'Job finished'
  END IF
  !PRINT *, "eigenvalues" 
  !print* , REALPART(freq_r)
  !print *, IMAGPART(freq_r)
  !print *, "Eigenvectors"
  !PRINT *, REALPART(VR)	
  !print *, 'lol'
 ! print *, IMAGPART(VR)


  ! CALL PZLAPRNT( 10, 10, Vr, 1, 1, DESCVR, 0, 0, 'VR', 6, WORK )
  
  IF (myrow==0 .and. mycol==0) THEN
   !  print *, 'eigenvalues'
    ! print *, DBLE(freq)   
    open(unit=43,file=TRIM(DIRDATA) // 'fort.43', &
    access='SEQUENTIAL')
    DO i=1,ntot
      WRITE(43,9999) freq(i)
    END DO
    CLOSE(43)
    
    open(unit=44,file=TRIM(DIRDATA) // 'fort.44', &
    access='SEQUENTIAL')
    WRITE(44,*) DBLE(freq), AIMAG(freq)
    close(44)
  END IF

  CALL PZLAWRITE(TRIM(DIRDATA) // 'scalapack.dat',ntot, &
       ntot,VR, 1,1,descVR,0,0,WORK)

  CALL PZLAWRITE(TRIM(DIRDATA) // 'eig_sca.dat',nlong*nlat, &
       ntot,VR, 1,1,descVR,0,0,WORK)
  
  
!   
!   
! 
!   CALL PZLAPRNT(1,1,VR,15,15,DESCVR,0,0,'Vr',6,WORK)
!   CALL PZELGET('A', ' ',alpha,VR,15,15,descVR)  
!   CALL PZLAPRNT(1,1,VR,15,15,DESCVR,0,0,'Vrbis',6,WORK)
!   print *, alpha
  
  
  !DEALLOCATE(freq,vl,vr)  
  !
  DEALLOCATE(IWORK)
  DEALLOCATE(WORK)
  !DEALLOCATE(LM)
  

9999 FORMAT( E15.8,E15.8 )

  END SUBROUTINE
  
END MODULE
