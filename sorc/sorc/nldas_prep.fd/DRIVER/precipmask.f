!=========================================================================
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
! precipmask.f:
!
! DESCRIPTION:
!  This subroutine reads in precipmask information
!
! REVISION HISTORY:
!  06 Sep 2001: Brian Cosgrove; Initial Code
!=========================================================================

      SUBROUTINE precipmask(LDAS,GRID)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! GRID variables
      IMPLICIT NONE
      type (ldasdec) LDAS
      type (griddec) grid(ldas%nc,ldas%nr)

	INTEGER R,C


C=== Open weighting file used in making merged precip
C=== Forces use of ETA/S4 precip in Canada and Mexico
      IF (LDAS%PRECIPMASK.EQ.1) THEN
        OPEN(32,FILE=LDAS%PRECIPWEIGHT,
     &  STATUS='OLD')
        DO R=1,LDAS%NR
         DO C=1,LDAS%NC
          READ (32,*) GRID(C,R)%PRECIPWEIGHT
         ENDDO
        ENDDO
        CLOSE(32)
      ENDIF
      END
