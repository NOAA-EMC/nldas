!
!========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale
!  D                                      L  Land Modeling and Data
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! trnsfm.f:
!
! DESCRIPTION:
!  Transforms from grid space to catchment space (opt=1) and catchment 
!  space to grid space (opt=2) using the pointer array.
!  Transforms from grid space to catchment space (opt=3) and catchment 
!  space to grid space (opt=4) using the transformation code.
!
! REVISION HISTORY:
! 4  Apr 2000: Jeffrey Walker; Initial Code
! 20 Jun 2000: Jeffrey Walker; Added pointer option for transformations
! 27 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!                not a hard-wired undefined value of -999.9
! 15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new 
!                GRIB libraries
! 14 Jan 2003: Urszula Jambor; Added deallocation statements near end of 
!                routine and changed pointer variables to allocatable.
!=========================================================================

      SUBROUTINE TRNSFM(OPT,LDAS,GRIDSP,CATCHSP,NFF)

! Declare modules and data structures
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      TYPE (ldasdec) LDAS
      INTEGER :: OPT,NFF(LDAS%NCATM)
      REAL :: GRIDSP(LDAS%NC,LDAS%NR),CATCHSP(LDAS%NCATM)

!=== Local variables =====================================================

      INTEGER :: I,J,N,NX,NY
      INTEGER, allocatable :: MASK1(:,:),MASK2(:,:),NUM(:)
      INTEGER, parameter :: NOCAT=-9999
      REAL, pointer :: CATCHSC(:)
      LOGICAL*1 :: FIRST1=.true.,FIRST2=.true.
      SAVE 

!=== End Variable List ===================================================
      ALLOCATE (CATCHSC(LDAS%NCAT))
      ALLOCATE (NUM(LDAS%NCAT))

!=== Transform using pointer array =======================================

      IF (OPT.EQ.1 .OR. OPT.EQ.2) THEN
       IF (FIRST1) THEN
        ALLOCATE (MASK1(LDAS%NC,LDAS%NR))
        OPEN(10,FILE='ldas.point.dat',FORM='FORMATTED',STATUS='UNKNOWN')
        DO I=1,LDAS%NR
         READ (10,*) (MASK1(J,I), J=1,LDAS%NC) 
        END DO
        CLOSE (10,STATUS='KEEP')
        FIRST1=.false.
       END IF

!=== Transform from grid to catchment ====================================

       IF (OPT .EQ. 1) THEN
        DO N=1,LDAS%NCAT
         CATCHSC(N)=0.0
         NUM(N)=0
        END DO
        DO I=1,LDAS%NR
         DO J=1,LDAS%NC
          IF (MASK1(J,I).NE.NOCAT .AND. GRIDSP(J,I).NE.LDAS%UDEF) THEN
           CATCHSC(MASK1(J,I))=CATCHSC(MASK1(J,I))+GRIDSP(J,I)
           NUM(MASK1(J,I))=NUM(MASK1(J,I))+1
          END IF
         END DO
        END DO

!=== Re-map to LDAS domain
        DO N=1,LDAS%NCATM
         IF (NUM(NFF(N)) .GT. 0) THEN
          CATCHSP(N)=CATCHSC(NFF(N))/FLOAT(NUM(NFF(N)))
         ELSE
          CATCHSP(N)=LDAS%UDEF
          PRINT *, 'Error: Undefined forcing for catchment ',N
         END IF
        END DO

!=== Transform from catchment to grid ====================================

       ELSE IF (OPT .EQ. 2) THEN

!=== Map from LDAS domain
        DO N=1,LDAS%NCAT
         CATCHSC(N)=LDAS%UDEF
        END DO
        DO N=1,LDAS%NCATM
         CATCHSC(NFF(N))=CATCHSP(N)
        END DO

!=== Perform transformation
        DO I=1,LDAS%NR
         DO J=1,LDAS%NC
          IF (MASK1(J,I) .NE. NOCAT) THEN
           GRIDSP(J,I)=CATCHSC(MASK1(J,I))
          ELSE
           GRIDSP(J,I)=LDAS%UDEF
          END IF
         END DO
        END DO
       END IF

!=== Transform using transformation code =================================

      ELSE IF (OPT .EQ. 3 .OR. OPT .EQ. 4) THEN
       OPEN(10,FILE='ldas.grd',FORM='FORMATTED',STATUS='UNKNOWN')
       READ (10,*) NX
       READ (10,*) NY
       CLOSE(10,STATUS='KEEP')
       IF (NX.NE.LDAS%NC .OR. NY.NE.LDAS%NR) THEN
        WRITE (*,*) 'Transformation specification data mismatch'
        STOP
       END IF

!=== Transform from grid to catchment ====================================

       IF (OPT .EQ. 3) THEN

!=== Write grid data
        OPEN(10,FILE='griddata.dat',FORM='FORMATTED',STATUS='UNKNOWN')
        DO I=1,LDAS%NR
         WRITE(10,*) (DBLE(GRIDSP(J,I)), J=1,LDAS%NC) 
        END DO
        CLOSE(10,STATUS='KEEP')

!=== Perform transformation
        CALL SYSTEM('rm grd2cat.log')
        CALL SYSTEM('mpirun -np 1 grd2cat > grd2cat.log') 

!=== Read catchment data
        OPEN(10,FILE='catchdata.dat',FORM='FORMATTED',STATUS='OLD')
        READ(10,*) (CATCHSC(N), N=1,LDAS%NCAT)
        CLOSE(10,STATUS='KEEP')

!=== Re-map to LDAS domain
        DO N=1,LDAS%NCATM
         CATCHSP(N)=CATCHSC(NFF(N))
        END DO

!=== Transform from catchment to grid ====================================

       ELSE IF (OPT .EQ. 4) THEN

!=== Map from LDAS domain
        DO N=1,LDAS%NCAT
         CATCHSC(N)=LDAS%UDEF
        END DO
        DO N=1,LDAS%NCATM
         CATCHSC(NFF(N))=CATCHSP(N)
        END DO

!=== Write catchment data
        OPEN(10,FILE='catchdata.dat',FORM='FORMATTED',STATUS='UNKNOWN')
        WRITE(10,*) (CATCHSC(N), N=1,LDAS%NCAT) 
        CLOSE(10,STATUS='KEEP')

!=== Perform transformation
        CALL SYSTEM('rm cat2grd.log')
        CALL SYSTEM('mpirun -np 1 cat2grd > cat2grd.log') 

!=== Read grid data mask
        IF (FIRST2) THEN
         OPEN(10,FILE='ldas.g.mask',FORM='FORMATTED',STATUS='OLD')
c because transformation is currently for 1/4 degree
         ALLOCATE (MASK2(LDAS%NC/2,LDAS%NR/2))
         DO I=1,LDAS%NR/2
          READ(10,*) (MASK2(J,I), J=1,LDAS%NC/2)
         END DO
         CLOSE (10,STATUS='KEEP')
         FIRST2=.FALSE.
        END IF        

!=== Read grid data and apply mask
        OPEN(10,FILE='griddata.dat',FORM='FORMATTED',STATUS='OLD')

c because transformation is currently for 1/4 degree
        DO I=1,LDAS%NR/2
         READ(10,*) (GRIDSP(J,I), J=1,LDAS%NC/2) 
        END DO
        DO I=1,LDAS%NR/2
         DO J=1,LDAS%NC/2
          IF (MASK2(J,I).NE.1) GRIDSP(J,I)=LDAS%UDEF
         END DO
        END DO

        do i=ldas%nr,2,-2
         do j=ldas%nc,2,-2
          gridsp(j,i)=gridsp(j/2,i/2)
          gridsp(j,i-1)=gridsp(j/2,i/2)
          gridsp(j-1,i)=gridsp(j/2,i/2)
          gridsp(j-1,i-1)=gridsp(j/2,i/2)
         end do
        end do

c for transformation at 1/8 degree
!      DO I=1,LDAS%NR
!       READ(10,*) (GRIDSP(J,I), J=1,LDAS%NC) 
!      END DO

        CLOSE(10,STATUS='KEEP')
       ELSE
        WRITE(*,*) 'Invalid Transformation Option'
        STOP
       END IF
      END IF

      DEALLOCATE (CATCHSC)
      DEALLOCATE (NUM)
      IF (ALLOCATED(MASK1)) DEALLOCATE(MASK1)
      IF (ALLOCATED(MASK2)) DEALLOCATE(MASK2)

      END SUBROUTINE TRNSFM

