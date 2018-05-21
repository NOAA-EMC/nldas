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
! stats.f:
!
! DESCRIPTION:
!  Calculates statistics for a given variable
!
! REVISION HISTORY:
! Nov 11 1999:  Jon Radakovich; Initial code
!=========================================================================

      SUBROUTINE STATS(VAR,UDEF,NCH,MEAN,STDEV,MIN,MAX)

!=== Local variables =====================================================
      INTEGER :: NCH,T
      REAL :: VAR(NCH)
      REAL :: MEAN,DEV,STDEV,MIN,MAX,UDEF,VSUM

!=== End Variable List ===================================================
      VSUM=0.
      MEAN=0.
      DEV=0.
      STDEV=0.
      MIN=100000.
      MAX=-100000.
      DO T=1,NCH
       IF(VAR(T).NE.UDEF)THEN
        VSUM=VSUM+VAR(T)
        IF(VAR(T).GT.MAX)MAX=VAR(T)
        IF(VAR(T).LT.MIN)MIN=VAR(T)
       ENDIF
      ENDDO
      IF(VSUM.EQ.0.)THEN
       MAX=0.
       MIN=0.
      ENDIF
      
      MEAN=VSUM/FLOAT(NCH)
      
      DO T=1,NCH
       IF(VAR(T).NE.UDEF)THEN
        DEV=DEV+(VAR(T)-MEAN)**2
       ENDIF
      ENDDO
      STDEV=(DEV*(FLOAT(NCH)-1)**(-1))**(0.5)
      RETURN
      END
      
      
