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
! !MODULE: geosdomain_module.F90
! 
! !DESCRIPTION: 
!  Contains routines and variables that define the native domain
!  for GEOS model forcing. 
! 
! !INTERFACE:
module geosdomain_module
! !USES:
  use geosdrv_module
! !ARGUMENTS:
  type(geosdrvdec) :: geosdrv
  integer :: mi
!EOP
contains
  
!BOP
!
! !ROUTINE: defnatgeos.F90
! 
! !DESCRIPTION: 
!  Defines the gridDesc array describing the native forcing resolution 
!  for GEOS data. 
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine defnatgeos(gridDesci)
! !USES: 
    use lisdrv_module, only: lis
    use time_manager, only : date2time
    implicit none
!EOP
!    integer, intent(inout) :: gridDesci(200)
    real, intent(inout) :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    
!BOC
    call readgeoscrd(geosdrv,gridDesci)
    mi = geosdrv%ncold*geosdrv%nrold
    yr1 = 2002  !grid update time
    mo1 = 10
    da1 = 01
    hr1 = 0; mn1 = 0; ss1 = 0
    call date2time(geosdrv%griduptime,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
!EOC
  end subroutine defnatgeos
end module geosdomain_module
