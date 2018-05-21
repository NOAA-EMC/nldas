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

!BOP
!
! !MODULE: noahdrv_module.F90 
!
! !DESCRIPTION:
!  Module for runtime specific Noah variables
!
! !REVISION HISTORY:
!
! 14 Oct 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
module noahdrv_module
! !ARGUMENTS:
  type noahdrvdec
     integer :: noahopen          !Keeps track of opening files
     integer :: numout            !Counts number of output times for Noah
     integer :: noah_nvegp        !Number of static vegetation parameter
     integer :: noah_nsoilp       !Number of static soil parameters
     integer :: noah_zst          !Number of Zobler soil classes
     integer :: noah_gflag        !Time flag to update gfrac files
     integer :: noah_albtime      !Time flag to update albedo files
     integer :: noah_aflag        !Time flag to update albedo files 
     integer :: noah_albdchk      !Day check to interpolate alb values
     integer :: noah_gfracdchk    !Day check to interpolate gfrac value
     integer :: albedo_type       !albedo source: 1=LIS, 2=GSWP-2
     integer :: gfrac_type        !gfrac source: 1=LIS, 2=GSWP-2
     integer :: mxsnalb_type      !mxsnalb source: 1=LIS, 2=GSWP-2
     integer :: tbot_type         !tbot source: 1=LIS, 2=GSWP-2
     integer :: varid(24)         !For netcdf output 
     character*40 :: noah_rfile   !NOAH Active Restart File
     character*40 :: noah_vfile   !NOAH Static Vegetation Parameter File
     character*40 :: noah_sfile   !NOAH Soil Parameter File
     character*40 :: noah_mgfile  !NOAH Monthly Veg. Green Frac.
     character*40 :: noah_albfile !NOAH Quart. Snow-free albedo
     character*50 :: noah_mxsnal  !NOAH GLDAS max snow albedo
     character*50 :: noah_tbot    !NOAH GLDAS Bottom Temp
     real*8 :: noah_gfractime     !Time flag to update gfrac files
     real :: noah_ism             !NOAH Initial Soil Moisture (m3/m3)
     real :: noah_it              !NOAH Initial Soil Temperature (K)
     real :: writeintn            !NOAH Output Interval (hours)
  end type noahdrvdec
!EOC
end module noahdrv_module
