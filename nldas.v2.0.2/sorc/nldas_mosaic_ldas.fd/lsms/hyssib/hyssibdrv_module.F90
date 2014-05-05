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
! !MODULE: hyssibdrv_module.f:
!
! !DESCRIPTION:
!  Module for runtime specific HY-SSiB variables
!
! !REVISION HISTORY:
! 14 Oct 2003: Sujay Kumar, Initial Version
! 23 Apr 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      MODULE hyssibdrv_module
!EOP
      type hyssibdrvdec
!BOC
      integer :: hyssibopen          !Keeps track of opening files
      integer :: numout            !Counts number of output times for HY-SSiB
      integer :: hyssib_nvegp        !Number of static vegetation parameter
      integer :: hyssib_nvegip       !Number of monthly vegetation parameter
      integer :: hyssib_flgres       !Restart flag
      integer :: hyssib_nsoilp       !Number of static soil parameters
      integer :: hyssib_zst          !Number of Zobler soil classes
      integer :: hyssib_gflag        !Time flag to update gfrac files
      integer :: hyssib_albtime      !Time flag to update albedo files
      integer :: hyssib_aflag        !Time flag to update albedo files 
      integer :: hyssib_albdchk      !Day check to interpolate alb values
      integer :: hyssib_gfracdchk    !Day check to interpolate gfrac value
      integer :: STATEVAR_AVG        !HY-SSiB Instantaneous (=0) or Time-Averaged (=1) output for state variables
      integer :: varid(46)
      CHARACTER*40 :: HYSSIB_RFILE   !HY-SSiB Active Restart File
      CHARACTER*40 :: HYSSIB_AFILE   !HY-SSiB Albedo and Radiation Parameter File
      CHARACTER*40 :: HYSSIB_VFILE   !HY-SSiB Static Vegetation Parameter File
      CHARACTER*40 :: HYSSIB_SFILE   !HY-SSiB Soil Parameter File
      CHARACTER*40 :: HYSSIB_MGFILE  !HY-SSiB Monthly Veg. Green Frac.
      CHARACTER*40 :: HYSSIB_ALBFILE !HY-SSiB Quart. Snow-free albedo
      CHARACTER*40 :: HYSSIB_GALBRES !HY-SSiB Quart. Snow-free albedo resolution
      CHARACTER*50 :: HYSSIB_MXSNAL  !HY-SSiB LIS max snow albedo
      CHARACTER*50 :: HYSSIB_TBOT    !HY-SSiB LIS Bottom Temp
      CHARACTER*50 :: HYSSIB_TOPOSTD !HY-SSiB LIS Standard Dev. of Topography
      character*80 :: ofile 
      REAL*8 :: HYSSIB_GFRACTIME     !Time flag to update gfrac files
      REAL :: HYSSIB_ISM             !HY-SSiB Initial Soil Moisture (m3/m3)
      REAL :: HYSSIB_IT              !HY-SSiB Initial Soil Temperature (K)
      REAL :: WRITEINT              !HY-SSiB Output Interval (hours)
      end type hyssibdrvdec
!EOC
      end MODULE hyssibdrv_module

