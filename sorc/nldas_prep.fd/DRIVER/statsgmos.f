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
! January 10 2001: Brian Cosgrove; Modified code so that it 
!                  only includes land points in calculations
!=========================================================================

      SUBROUTINE STATSGMOS(LDAS,GRID,UDEF,VAR,NC,NR,MEAN,STDEV,MIN,MAX)
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE GRID_MODULE      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      TYPE (ldasdec) LDAS
      TYPE (GRIDDEC) GRID(LDAS%NC,LDAS%NR)

!=== Local variables =====================================================
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
        IF (GRID(C,R)%MASK.NE.0.0) THEN
        IF(VAR(C,R).NE.UDEF)THEN
         SUM=SUM+VAR(C,R)
         IF(VAR(C,R).GT.MAX)MAX=VAR(C,R)
         IF(VAR(C,R).LT.MIN)MIN=VAR(C,R)
         N=N+1
        ENDIF
        ENDIF
       ENDDO
      ENDDO

	IF (N.GT.0) THEN
      MEAN=SUM/FLOAT(N)
      DO R=1,NR
       DO C=1,NC
        IF (GRID(C,R)%MASK.NE.0.0) THEN
        IF(VAR(C,R).NE.UDEF)THEN
         DEV=DEV+(VAR(C,R)-MEAN)**2
         if (dev.eq.0.0) then
c          print *,'PROBLEM, DEV=0 in stats, setting to .01'
         DEV=0.01
	 endif
        ENDIF
        ENDIF
       ENDDO
      ENDDO
        if (float(n).eq.1) then
c        print *,'PROBLEM, N=1 in stats, setting to 2'
        N=2
        endif

      STDEV=(DEV*(FLOAT(N)-1)**(-1))**(0.5)
	ELSE
	MEAN=-777
	DEV=-777
	STDEV=-777
	ENDIF
      RETURN
      END
