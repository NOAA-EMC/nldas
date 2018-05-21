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
! mos_initout.f: 
!
! DESCRIPTION:
!  Initialize Mosaic output arrays
!
! REVISION HISTORY:
! 4 Nov. 1999: Paul Houser; Initial Code
!=========================================================================

      SUBROUTINE MOS_INITOUT(LDAS,MOS)
 
! Declare modules and data structures
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE mos_module        
      IMPLICIT NONE
      TYPE (ldasdec) LDAS
      TYPE (mosdec) MOS(LDAS%NCH)

!=== Local variables =====================================================
      INTEGER :: T,N                ! Tile loop counter
      INTEGER :: m

!=== End Variable List ===================================================

      DO T=1,LDAS%NCH
       ALLOCATE (MOS(T)%RETURN(LDAS%MOS_NRET))
       ALLOCATE (MOS(T)%TOTRET(LDAS%MOS_NRET))
      ENDDO

      MOS%COUNT=0

      DO N=1,LDAS%MOS_NRET
       DO T=1,LDAS%NCH
        MOS(T)%RETURN(N)=0.0
        MOS(T)%TOTRET(N)=0.0
       ENDDO
      ENDDO
      
      END SUBROUTINE MOS_INITOUT













