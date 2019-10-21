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
! !MODULE: nldas2drv_module.F90 
!
! !DESCRIPTION:
!  Module containing runtime specific NLDAS2 variables
!
! !REVISION HISTORY:
!
! 02Feb2004; Sujay Kumar, Initial Version
! 
! !INTERFACE:
module nldas2drv_module
! !ARGUMENTS:
  type nldas2drvdec
     integer :: ncold, nrold   !AWIPS 212 dimensions
     integer :: nmif
     character*40 :: nldas2dir   !NLDAS2 Forcing Directory
     character*40 :: blend   !NLDAS2 Forcing Directory
     real*8 :: nldas2time1,nldas2time2
  end type nldas2drvdec
!EOC
end module nldas2drv_module
