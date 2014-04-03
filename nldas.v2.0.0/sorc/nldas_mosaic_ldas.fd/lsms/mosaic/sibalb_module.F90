!-------------------------------------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
! Released October 2005
!
! See SOFTWARE DISTRIBUTION POLICY for software distribution policies
!
! The LIS source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See COPYRIGHT.TXT for copyright details.
!
!-------------------------------------------------------------------------
!
! DESCRIPTION:
!  In order to use the MOSAIC-SiB albedo calculation, SiB vegetation types
!  are mapped to UMD vegetation types.
!  For both SiB and UMD, 4 coefficient arrays are defined for calculating 
!  the albedos of soil in the subroutine sibalb, representing the 
!  following bands:
!  visible, direct solar radiation,
!  infra-red, direct solar radiation,
!  visible, diffuse solar radiation, and
!  infra-red, diffuse solar radiation.
!
! MOSAIC/SiB ITYP: Vegetation type as follows:
!                  1:  BROADLEAF EVERGREEN TREES
!                  2:  BROADLEAF DECIDUOUS TREES
!                  3:  NEEDLELEAF TREES
!                  4:  GROUND COVER
!                  5:  BROADLEAF SHRUBS
!                  6:  DWARF TREES (TUNDRA)
!                  7:  BARE SOIL
!
! UMD ITYP: Vegetation type as follows:
!                  1.  Evergreen Needleleaf Forest
!                  2.  Evergreen Broadleaf Forest
!                  3.  Deciduous Needleleaf Forest
!                  4.  Deciduous Broadleaf Forest
!                  5.  Mixed Cover
!                  6.  Woodland
!                  7.  Wooded Grassland
!                  8.  Closed Shrubland
!                  9.  Open Shrubland
!                 10.  Grassland
!                 11.  Cropland
!                 12.  Bare Ground
!                 13.  Urban and Built-Up
!
!===========================================================================
! REVISION HISTORY:
!  04 Apr 2001: Urszula Jambor; Initial code, using old calc_albedo.f scheme
!===========================================================================

      MODULE sibalb_module 

      IMPLICIT NONE
      public sibalbdec

      type sibalbdec 

!=== Constants used in albedo calculations: =========

!=== Original SiB albedo coefficients and SiB mapped to UMD coefficients
!=== dimension (NLAI=14, 2, NTYPS=9)
!=== dimension (NLAI=14, 2, UMDNTYPS=13)

	real, dimension(14, 2, 9) :: ALVDRold, BTVDRold, GMVDRold
	real, dimension(14, 2, 9) :: ALIDRold, BTIDRold, GMIDRold

	real, dimension(14, 2, 13) :: ALVDR, BTVDR, GMVDR
	real, dimension(14, 2, 13) :: ALIDR, BTIDR, GMIDR

      end type !sibalbdec
      
      type (sibalbdec) sib
      
      end module sibalb_module


