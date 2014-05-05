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
! noah_initout.f: 
!
! DESCRIPTION:
!  Initialize NOAH output arrays
!
! REVISION HISTORY:
! 4 Nov. 1999: Paul Houser; Initial Code
! 28 Apr. 2002: K. Arsenault; Modified to NOAH LSM 2.5 code to LDAS
!=========================================================================

      SUBROUTINE NOAH_INITOUT(LDAS,NOAH)
 
! Declare modules and data structures
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE noah_module      ! NOAH LSM module  
      IMPLICIT NONE
      TYPE (ldasdec) LDAS
      TYPE (noahdec) NOAH(LDAS%NCH)

!=== Local variables =====================================================
      INTEGER :: T,N                ! Tile loop counter
      INTEGER :: m

!=== End Variable List ===================================================

      DO T=1,LDAS%NCH
       ALLOCATE (NOAH(T)%RETURN(LDAS%NOAH_NRET))
       ALLOCATE (NOAH(T)%TOTRET(LDAS%NOAH_NRET))
      ENDDO

      NOAH%COUNT=0

      DO N=1,LDAS%NOAH_NRET
       DO T=1,LDAS%NCH
        NOAH(T)%RETURN(N)=0.0
        NOAH(T)%TOTRET(N)=0.0
       ENDDO
      ENDDO
      
      END SUBROUTINE NOAH_INITOUT













