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
! FORCE2TILE.f: 
!
! DESCRIPTION:
!  Transfer forcing from grid to tile space.
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Code
!  28 Jan 2002: Jon Gottschalck; Added option for different number of forcing variables  
!=========================================================================

      SUBROUTINE FORCE2TILE(LDAS,TILE,GRID)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (tiledec) TILE(LDAS%NCH)
      type (griddec) GRID(LDAS%NC,LDAS%NR)   

!=== Local Variables =====================================================
      INTEGER :: F,C,R,T,I,J     ! Loop counters
      INTEGER :: NFORCE          ! Number of forcing variables
!=== End Variable Definition =============================================

!=== Looping through the right number of forcing variables
      IF (LDAS%TSCOUNT .EQ. 0 .AND. LDAS%STARTCODE .EQ. 4) THEN
        NFORCE = LDAS%NMIF
      ELSE
        NFORCE = LDAS%NF
      ENDIF

      DO F=1,NFORCE
        DO T=1,LDAS%NCH
          TILE(T)%FORCING(F)=GRID(TILE(T)%COL,TILE(T)%ROW)%FORCING(F)
        ENDDO
      ENDDO

      RETURN
      END



































