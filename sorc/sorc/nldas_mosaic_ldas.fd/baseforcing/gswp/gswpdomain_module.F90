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
! !MODULE: gswpdomain_module.F90
! 
! !DESCRIPTION: 
!  Contains routines and variables that define the native domain
!  for GSWP model forcing. 
! 
! !INTERFACE:
module gswpdomain_module
! !USES:
  use gswpdrv_module
! !ARGUMENTS:
  type(gswpdrvdec) :: gswpdrv
  integer :: mi
!EOP
contains
  
!BOP
!
! !ROUTINE: defnatgswp.F90
! 
! !DESCRIPTION: 
!  Defines the gridDesc array describing the native forcing resolution 
!  for GSWP data. 
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine defnatgswp(gridDesci)
! !USES: 
    use lisdrv_module, only: lis
    use time_manager, only : date2time
    implicit none
!EOP
    real, intent(inout) :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    
!BOC
    call readgswpcrd(gswpdrv,gridDesci)
    mi = gswpdrv%ncold*gswpdrv%nrold
    gswpdrv%gswptime1 = 3000.0
    gswpdrv%gswptime2 = 0.0
!EOC
  end subroutine defnatgswp
end module gswpdomain_module
