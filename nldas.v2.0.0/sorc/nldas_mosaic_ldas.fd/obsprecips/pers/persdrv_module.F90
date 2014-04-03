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
! !MODULE: persdrv_module.F90 
!
! !DESCRIPTION:
!  Module for runtime specific PERS variables
!
! !REVISION HISTORY:
!
! 11 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
module persdrv_module
! !ARGUMENTS:
  type persdrvdec
     integer :: ncold, nrold   !AWIPS 212 dimensions
     character*40 :: persdir   !PERS Forcing Directory
     real*8 :: perstime
  end type persdrvdec
!EOC
end module persdrv_module
