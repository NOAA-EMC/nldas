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
! statsg.f:
!
! DESCRIPTION:
!  Calculates statistics for a given grid variable
!
! REVISION HISTORY:
! April 6 2000:  Jon Radakovich; Initial code
!=========================================================================

      SUBROUTINE STATSG(UDEF,VAR,NC,NR,MEAN,STDEV,MIN,MAX)

!=== Local variables =====================================================
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      TYPE (ldasdec) LDAS
      INTEGER :: NC,NR,C,R,N
      REAL,DIMENSION(NC,NR) :: VAR
      REAL :: SUM,MEAN,DEV,STDEV,MIN,MAX,UDEF

!=== End Variable List ===================================================
      SUM=0.
      MEAN=0.
      DEV=0.
      STDEV=0.
      MIN=100000.
      MAX=-100000.
      N=0

      DO R=1,NR
       DO C=1,NC
        IF(VAR(C,R).NE.UDEF)THEN
         SUM=SUM+VAR(C,R)
         IF(VAR(C,R).GT.MAX)MAX=VAR(C,R)
         IF(VAR(C,R).LT.MIN)MIN=VAR(C,R)
         N=N+1
        ENDIF
       ENDDO
      ENDDO

      MEAN=SUM/FLOAT(N)

      DO R=1,NR
       DO C=1,NC
        IF(VAR(C,R).NE.UDEF)THEN
         DEV=DEV+(VAR(C,R)-MEAN)**2
        ENDIF
       ENDDO
      ENDDO

      STDEV=(DEV*(FLOAT(N)-1)**(-1))**(0.5)


      RETURN
      END
      
      
