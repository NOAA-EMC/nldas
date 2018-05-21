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
! !MODULE: gswpdrv_module.F90 
!
! !DESCRIPTION:
!  Module containing runtime specific GSWP variables
!
! !REVISION HISTORY:
!
! 20Feb2004; Sujay Kumar, Initial Version
! 
! !INTERFACE:
module gswpdrv_module
! !ARGUMENTS:
  type gswpdrvdec
     integer :: ncold, nrold   !AWIPS 212 dimensions
     integer :: nmif
     real*8 :: gswptime1,gswptime2
     character*100 :: tair
     character*100 :: qair
     character*100 :: psurf
     character*100 :: wind
     character*100 :: rainf
     character*100 :: snowf
     character*100 :: swdown
     character*100 :: lwdown
     character(len=40) :: albedo
     character(len=40) :: gfrac
     !character*40 :: gswpdir   !GEOS Forcing Directory
  end type gswpdrvdec
!EOC
end module gswpdrv_module
