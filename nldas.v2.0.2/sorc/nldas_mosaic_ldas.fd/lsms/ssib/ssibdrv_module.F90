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
! !MODULE: ssibdrv_module.f:
!
! !DESCRIPTION:
!  Module for runtime specific SSiB variables
!
! !REVISION HISTORY:
! 14 Oct 2003: Sujay Kumar, Initial Version
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
! 19 May 2004: David Mocko, made compatible with SiB-lings
!
! !INTERFACE:
      MODULE ssibdrv_module
!EOP
      type ssibdrvdec
!BOC
      integer :: ssibopen          !Keeps track of opening files
      integer :: numoutnh          !Counts number of output times for SSiB
      integer :: ssib_nvegp        !Number of static vegetation parameter
      integer :: ssib_nvegip       !Number of monthly vegetation parameter
      integer :: ssib_flgres       !Restart flag
      integer :: ssib_nsoilp       !Number of static soil parameters
      integer :: ssib_zst          !Number of Zobler soil classes
      integer :: ssib_gflag        !Time flag to update gfrac files
      integer :: ssib_albtime      !Time flag to update albedo files
      integer :: ssib_aflag        !Time flag to update albedo files 
      integer :: ssib_albdchk      !Day check to interpolate alb values
      integer :: ssib_gfracdchk    !Day check to interpolate gfrac value
      integer :: STATEVAR_AVG      !Instantaneous (=0) or Time-Averaged (=1) output for state variables
      CHARACTER*40 :: SSIB_RFILE   !SSIB Active Restart File
      CHARACTER*40 :: SSIB_MFILE   !SSIB model init. restart file
      CHARACTER*40 :: SSIB_VFILE   !SSIB Static Vegetation Parameter File
      CHARACTER*40 :: SSIB_SFILE   !SSIB Soil Parameter File
      CHARACTER*40 :: SSIB_MGFILE  !SSIB Monthly Veg. Green Frac.
      CHARACTER*40 :: SSIB_ALBFILE !SSIB Quart. Snow-free albedo
      CHARACTER*50 :: SSIB_MXSNAL  !SSIB GLDAS max snow albedo
      CHARACTER*50 :: SSIB_TBOT    !SSIB GLDAS Bottom Temp
      REAL*8 :: SSIB_GFRACTIME     !Time flag to update gfrac files
      REAL :: SSIB_ISM             !SSIB Initial Soil Moisture (m3/m3)
      REAL :: SSIB_IT              !SSIB Initial Soil Temperature (K)
      REAL :: WRITEINTN            !SSIB Output Interval (hours)
      end type ssibdrvdec
!EOC
      end MODULE ssibdrv_module

