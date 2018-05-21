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
! !MODULE: ecmwfdomain_module.F90
! 
! !DESCRIPTION: 
!  Contains routines and variables that define the native domain
!  for ECMWF model forcing. 
! 
! !INTERFACE:
module ecmwfdomain_module
! !USES:
  use ecmwfdrv_module
! !ARGUMENTS:
  type(ecmwfdrvdec) :: ecmwfdrv
  integer :: mi
!EOP
contains
  
!BOP
!
! !ROUTINE: defnatecmwf.F90
! 
! !DESCRIPTION: 
!  Defines the gridDesc array describing the native forcing resolution 
!  for ECMWF data. 
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine defnatecmwf(gridDesci)
! !USES: 
    use lisdrv_module, only: lis
    implicit none
!EOP
    real, intent(inout) :: gridDesci(50)
!BOC
    call readecmwfcrd(ecmwfdrv,gridDesci)
    mi = ecmwfdrv%ncold*ecmwfdrv%nrold

!EOC
  end subroutine defnatecmwf
end module ecmwfdomain_module
