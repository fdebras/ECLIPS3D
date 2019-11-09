MODULE mod_eigenvalues

  USE mod_init_para
  USE mod_data
  USE mod_init_matrix
  USE mod_fill_matrix
  USE mod_writing
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


  CALL DESCINIT(descVR, ntot,ntot,nb,nb,0,0,info_txt,nlld,info)
  CALL DESCINIT(descVL, ntot,ntot,nb,nb,0,0,info_txt,nlld,info)

  ALLOCATE(vr(nlld,nlld))
  ALLOCATE(vl(nlld,nlld))
  ALLOCATE(freq(ntot))

  IF (myrow==0 .and. mycol==0) THEN
    print *, 'Begin eigenvalues finder !'
    print *, nlld
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


  IF (myrow==0 .and. mycol==0) THEN
    print *, 'Job finished'
  END IF

  IF (myrow==0 .and. mycol==0) THEN
   open(unit=43,file=TRIM(DIRDATA) // 'frequency.dat', &
    access='SEQUENTIAL')
    DO i=1,ntot
      WRITE(43,9999) freq(i)
    END DO
    CLOSE(43)
  END IF

  CALL PZLAWRITE(TRIM(DIRDATA) // 'scalapack.dat',ntot, &
       ntot,VR, 1,1,descVR,0,0,WORK)

  CALL PZLAWRITE(TRIM(DIRDATA) // 'eig_sca.dat',nlat*nz, &
       ntot,VR, nz*(2*nlat+1)+1,1,descVR,0,0,WORK)

 !
  DEALLOCATE(IWORK)
  DEALLOCATE(WORK)


9999 FORMAT( E15.8,E15.8 )

  END SUBROUTINE

END MODULE
