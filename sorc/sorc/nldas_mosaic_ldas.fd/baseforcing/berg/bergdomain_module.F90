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
! !MODULE: bergdomain_module.F90
! 
! !DESCRIPTION: 
!  Contains routines and variables that define the native domain for 
!  BERG model forcing
! 
! !INTERFACE:
module bergdomain_module
! !USES:
  use bergdrv_module
!EOP  
  type(bergdrvdec) :: bergdrv
  integer :: mi
contains
  
!BOP
!
! !ROUTINE: defnatberg.F90
! 
! !DESCRIPTION: 
!  Defines the gridDesc array describing the native forcing resolution 
!  for BERG data. 
!
! !REVISION HISTORY: 
! 26Jan2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine defnatberg(gridDesci)
! !USES: 
    implicit none
! !ARGUMENTS:
    real, intent(inout) :: gridDesci(50)
!EOP
!BOC
    call readbergcrd(bergdrv)
    gridDesci(1) = 0
    gridDesci(2) = bergdrv%ncold
    gridDesci(3) = bergdrv%nrold
    gridDesci(4) = -89.750
    gridDesci(5) = -179.750
    gridDesci(6) = 128
    gridDesci(7) = 89.750
    gridDesci(8) = 179.750
    gridDesci(9) = 0.500
    gridDesci(10) = 0.500
    gridDesci(20) = 255
    call readbergmask()
    mi = bergdrv%ncold*bergdrv%nrold
    print*,'defnatberg mi ',mi
!EOC
  end subroutine defnatberg
end module bergdomain_module
