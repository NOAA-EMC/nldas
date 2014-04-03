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
! !MODULE: nldas2domain_module.F90
! 
! !DESCRIPTION: 
!  Contains routines and variables that define the native domain
!  for NLDAS2 model forcing. 
! 
! !INTERFACE:
module nldas2domain_module
! !USES:
  use nldas2drv_module
! !ARGUMENTS:
  type(nldas2drvdec) :: nldas2drv
  integer :: mi
!EOP
contains
  
!BOP
!
! !ROUTINE: defnatnldas2.F90
! 
! !DESCRIPTION: 
!  Defines the gridDesc array describing the native forcing resolution 
!  for NLDAS2 data. 
!
! !REVISION HISTORY: 
! 02Feb2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine defnatnldas2(gridDesci)
! !USES: 
    use lisdrv_module, only: lis
    use time_manager, only : date2time
    implicit none
!EOP
    real, intent(inout) :: gridDesci(50)
!BOC
    call readnldas2crd(nldas2drv,gridDesci)
    mi = nldas2drv%ncold*nldas2drv%nrold
!EOC
  end subroutine defnatnldas2
end module nldas2domain_module
