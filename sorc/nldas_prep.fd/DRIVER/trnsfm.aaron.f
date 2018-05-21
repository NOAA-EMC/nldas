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
! trnsfm.aaron.f:
!
! DESCRIPTION:
!  Transforms from catchment space to grid space using the pointer array.
!
! REVISION HISTORY:
! 27 Feb 2001: Brian Cosgrove; Initial code based heavily on trnsfm.f
! 15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new 
!                GRIB libraries
! 14 Jan 2003: Urszula Jambor; Added deallocation statements near end of 
!                routine and changed pointer variables to allocatable.
!=========================================================================

      SUBROUTINE TRNSFMAARON(LDAS,GRIDSP,CATCHSP)

! Declare modules and data structures
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      TYPE (ldasdec) LDAS
      REAL :: GRIDSP(LDAS%NC,LDAS%NR),CATCHSP(5018)

!=== Local variables =====================================================

      INTEGER :: I,J,N,NX,NY
      INTEGER, allocatable :: MASK1(:,:),MASK2(:,:),NUM(:)
      INTEGER, parameter :: NOCAT=-9999
      REAL, pointer :: CATCHSC(:)
      LOGICAL*1 :: FIRST1=.true.,FIRST2=.true.
      SAVE 

!=== End Variable List ===================================================
      ALLOCATE (CATCHSC(5018))
      ALLOCATE (NUM(5018))

!=== Transform using pointer array =======================================
       IF (FIRST1) THEN
        ALLOCATE (MASK1(LDAS%NC,LDAS%NR))
        OPEN(10,FILE='ldas.point.mod.dat',FORM='FORMATTED',STATUS=
     &  'UNKNOWN')
        DO I=1,LDAS%NR
         READ (10,*) (MASK1(J,I), J=1,LDAS%NC) 
        END DO
        CLOSE (10,STATUS='KEEP')
        FIRST1=.false.
       END IF

!=== Transform from catchment to grid ====================================
!=== Perform transformation
        DO I=1,LDAS%NR
         DO J=1,LDAS%NC
          IF (MASK1(J,I) .NE. NOCAT) THEN
           GRIDSP(J,I)=CATCHSP(MASK1(J,I))
          ELSE
           GRIDSP(J,I)=LDAS%UDEF
          END IF
         END DO
        END DO

!=== Transform using transformation code =================================

      DEALLOCATE (CATCHSC)
      DEALLOCATE (NUM)
      IF (ALLOCATED(MASK1)) DEALLOCATE(MASK1)
      IF (ALLOCATED(MASK2)) DEALLOCATE(MASK2)

      END SUBROUTINE TRNSFMAARON

