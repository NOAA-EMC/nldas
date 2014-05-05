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
! tile_module.f: 
!
! DESCRIPTION:
!  LDAS non-model-specific tile variables only.
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial code
!  22 Aug 2000: Brian Cosgrove; Modified code for output of
!               standard LDAS output variables--added CC and AC,
!               the canopy and aerodynamic conductance
!=========================================================================

      MODULE tile_module 

      IMPLICIT NONE
      public tiledec

      type tiledec

!=== LDAS Non-Model-Specific Tile Variables ==============================

!=== LDAS Tile Definition Variables ======================================
      INTEGER :: COL        !Grid Column of Tile
      INTEGER :: ROW        !Grid Row of Tile
      REAL    :: FGRD       !Fraction of Grid covered by tile 
      REAL    :: LON        !Longitude of tile
      REAL    :: LAT        !Latitude of tile
      INTEGER :: VEGT       !Vegetation Type of Tile
      INTEGER :: SOILT      !Soil Type of Tile (May need DIM in future)
      INTEGER :: PVEG       !Predominance of vegetation clas in grid
      REAL    :: CC         !Canopy Conductance of Tile
      REAL    :: AC         !Aerodynamic Conductance of Tile

      REAL, pointer :: FORCING(:) !LDAS forcing array in tile space

      end type
      end module tile_module










