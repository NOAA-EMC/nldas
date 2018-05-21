#include <misc.h>

module drv_tilemodule 
!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
!  Module for tile space variable specification.
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!  21 Mar 2002: Brian Cosgrove; added initial condition vars for init from
!               Dag Lohmanns LDAS conditions
!=========================================================================
! $Id: drv_tilemodule.F90,v 1.1.1.1 2003/02/06 16:10:44 jgottsch Exp $
!=========================================================================

  use precision
  use clm1_varpar, only : nlevsoi
  implicit none

  public clm_tiledec
  type clm_tiledec

!=== TILE SPACE User-Defined Parameters ====================================

     integer  :: col           ! Grid Column of Tile
     integer  :: row           ! Grid Row of Tile
     integer  :: vegt          ! Vegetation Type of Tile
     integer  :: pveg          ! Predominance of vegetation clas in grid
     real(r8) :: fgrd          ! Fraction of grid covered by a given veg type (%/100)

     real(r8) :: sand(nlevsoi) ! Percent sand in soil (vertically average)
     real(r8) :: clay(nlevsoi) ! Percent clay in soil (vertically average)

     real(r8) :: scalez        ! Soil layer thickness discretization (m)
     real(r8) :: hkdepth       ! length scale for Ksat decrease(m)
     real(r8) :: roota         ! Temporary CLM vegetation parameter
     real(r8) :: rootb         ! Temporary CLM vegetation parameter

     real(r8) :: IT1           !Initial Temp 1 from Dag/NOAA
     real(r8) :: IT2           !Initial Temp 2 from Dag/NOAA
     real(r8) :: IWET(nlevsoi) !Initial Soil Volumetric Mst based on Dag/NOAA wetness data


!=== End Variable List ===================================================

  end type clm_tiledec

end module drv_tilemodule
