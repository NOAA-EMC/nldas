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
! !MODULE: bergdrv_module.F90 
!
! !DESCRIPTION:
!  Module containing runtime specific BERG variables
!
! !REVISION HISTORY:
!
! 26 Jan 2004; Sujay Kumar, Initial Version
! 
! !INTERFACE:
module bergdrv_module
! !ARGUMENTS:
  type bergdrvdec
     integer      :: ncold, nrold   !AWIPS 212 dimensions
     integer      :: nmif
     real*8       :: fmodeltime1,fmodeltime2
     integer      :: remask1d(720*360)
     character*40 :: emaskfile !1/2deg Reanal-BERG Land-Sea Mask File
  end type bergdrvdec
!EOC
end module bergdrv_module
