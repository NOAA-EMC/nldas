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
! !MODULE: geosdrv_module.F90 
!
! !DESCRIPTION:
!  Module containing runtime specific GEOS variables
!
! !REVISION HISTORY:
!
! 11 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
module geosdrv_module
! !ARGUMENTS:
  type geosdrvdec
     integer :: ncold, nrold   !AWIPS 212 dimensions
     integer :: nmif
     character*40 :: geosdir   !GEOS Forcing Directory
     real*8 :: geostime1,geostime2
     real*8 :: griduptime
     logical :: gridchange
  end type geosdrvdec
!EOC
end module geosdrv_module
