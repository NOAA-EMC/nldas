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
! !MODULE: gdasdomain_module.F90
! 
! !DESCRIPTION: 
!  Contains routines and variables that define the native domain for 
!  GDAS model forcing
! 
! !INTERFACE:
module gdasdomain_module
! !USES:
  use gdasdrv_module
!EOP  
  type(gdasdrvdec) :: gdasdrv
  integer :: mi
contains
  
!BOP
!
! !ROUTINE: defnatgdas.F90
! 
! !DESCRIPTION: 
!  Defines the gridDesc array describing the native forcing resolution 
!  for GDAS data. 
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine defnatgdas(gridDesci)
! !USES: 
    use lisdrv_module, only: lis
    use time_manager, only : date2time
    implicit none
! !ARGUMENTS:
    real, intent(inout) :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
!EOP
!BOC
    call readgdascrd(gdasdrv)
    gridDesci(1) = 4
    gridDesci(2) = 384
    gridDesci(3) = 190
    gridDesci(4) = 89.277
    gridDesci(5) = 0
    gridDesci(6) = 128
    gridDesci(7) = -89.277
    gridDesci(8) = -0.938
    gridDesci(9) = 0.938
    gridDesci(10) = 95
    gridDesci(20) = 255
    mi = gdasdrv%ncold*gdasdrv%nrold

    yr1 = 2000
    mo1 = 01
    da1 = 24
    hr1 = 12
    mn1 = 0; ss1 = 0
    call date2time( gdasdrv%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
    
    yr1 = 2002     !grid update time
    mo1 = 10
    da1 = 29
    hr1 = 12
    mn1 = 0; ss1 = 0
    call date2time(gdasdrv%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

    yr1 = 2005     !grid update time
    mo1 = 05
    da1 = 31
    hr1 = 12
    mn1 = 0; ss1 = 0
    call date2time(gdasdrv%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

    gdasdrv%gridchange1 = .true.
    gdasdrv%gridchange2 = .true.
    gdasdrv%gridchange3 = .true.
!EOC
  end subroutine defnatgdas
end module gdasdomain_module
