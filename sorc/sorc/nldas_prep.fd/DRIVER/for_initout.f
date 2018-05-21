!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://LDAS.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! for_initout.f: 
!
! DESCRIPTION:
!  Initialize LDAS forcing output arrays
!
! REVISION HISTORY:
! 4 Nov. 1999: Jon Radakovich; Initial Code
! 4 Sep. 2001: Brian Cosgrove; Added more variables to be initialized
!=========================================================================

      SUBROUTINE FOR_INITOUT(LDAS, GRID)
 
! Declare modules and data structures
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      TYPE (LDASdec) LDAS
      TYPE (griddec) GRID(LDAS%NC,LDAS%NR)

      INTEGER I,C,R
      DO c=1,LDAS%NC
       DO r=1,LDAS%NR
        DO i=1,LDAS%NF
         IF(ldas%feta.eq.1)then
          GRID(c,r)%ETADATA1(i)=0.0
          GRID(c,r)%ETADATA2(i)=0.0
	 ENDIF
         IF((ldas%fncep.eq.1).or.(ldas%fnasa.eq.1))then
          GRID(c,r)%NCEPDATA1(i)=0.0
          GRID(c,r)%NCEPDATA2(i)=0.0
	 ENDIF 
        ENDDO
       ENDDO
      ENDDO
       
      GRID%NESDATA1=0.0
      GRID%NESDATA2=0.0
      GRID%PINKDATA1=0.0
      GRID%PINKDATA2=0.0

      GRID%TOTALTMP=0.
      GRID%TOTALSPFH=0.
      GRID%TOTALDSWRF=0.
      GRID%TOTALDLWRF=0.
      GRID%TOTALUGRD=0.
      GRID%TOTALVGRD=0.
      GRID%TOTALPRES=0.
      GRID%TOTALAPCP=0.
      GRID%TOTALACPCP=0.

      GRID%TOTALPAR=0.
      GRID%TOTALODSWRF=0.
      GRID%TOTALOAPCP=0.


      GRID%COUNTFOR=0

      
      END SUBROUTINE FOR_INITOUT
