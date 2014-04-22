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
#include "misc.h"
!BOP
!
! !MODULE: gswp_module.F90
! 
! !DESCRIPTION:
!   This module contains useful routines that generates indices 
!   for reading the GSWP data
! 
! !REVISION HISTORY: 
!  24Feb04    Sujay Kumar  Initial Specification
! 
! !INTERFACE:
module gswp_module
    
contains
  subroutine getgswp_monindex(yr,mo,index)
    integer, intent(out) :: index
    integer, intent(in) :: yr, mo
    integer :: k 
    logical :: leap 
    index = 0
    index = index + (yr-1982)*12 + mo
  end subroutine getgswp_monindex
  
#if 0 
  subroutine getgswp_timeindex(yr,mo,da,hr,index)
    
    implicit none
    integer, intent(in)  :: yr, mo, da, hr
    integer, intent(out) :: index
    integer :: ryr, rmo, rda, rhr, days1(12),days2(12),yrs(13)
    integer :: k
    data yrs /184,365,366,365,365,365,366,365,365,365,366,365,365/ !from 1983-94
    data days1 /31,28,31,30,31,30,31,31,30,31,30,31/
    data days2 /31,29,31,30,31,30,31,31,30,31,30,31/
    logical :: leap 
    ryr = 1982
    rmo = 7
    rda = 1
    rhr = 0
    
    index = 0
    if(yr.gt.1982) then 
       if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
            .or.(mod(yr,400).eq.0))then 
          leap = .true.                   
       else
          leap = .false.
       endif
       k = mo
       index = index + hr+ (da-1)*24 
       do while(k.gt.1)
          if(leap) then 
             index = index + days2(k-1) *24  !current year's
          else 
             index = index + days1(k-1) *24  !current year's
          endif
          k = k-1
       enddo
       ! now add the missing years..
       k = yr-1982
       do while(k.ge.1) 
          index = index + yrs(k)*24
          k = k-1
       enddo
       index = index/3
    else
       if(mo.gt.7) then 
          k=mo
          index = index+hr+(da-1)*24
          do while(k .gt.7) 
             index = index + days1(k-1)*24
             k = k-1
          end do
       else
          index = hr+(da-1)*24
       endif
       index = index/3
    endif
    index = index+1
    !convert times to 3 hour index
    
  end subroutine getgswp_timeindex
#endif  
  subroutine getgswp_timeindex(yr,mo,da,hr,mn,ss,index)

    implicit none
    integer, intent(in)  :: yr, mo, da, hr, mn, ss
    integer, intent(out) :: index
    integer :: ryr, rmo, rda, rhr, days1(12)
    integer :: tmp_da, tmp_hr
    data days1 /31,28,31,30,31,30,31,31,30,31,30,31/
    logical :: leap 
    ryr = 1982
    rmo = 7
    rda = 1
    rhr = 3
  
    index = 0
    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then
       leap = .true.                   
    else
       leap = .false.
    endif

    if ( leap ) then
       days1(2) = 29
    endif

    tmp_da = da
    tmp_hr = hr

    ! Note:  GSWP forcing data are written into monthly files with
    ! a frequency of 3 hours.  The first entry corresponds to 03:00:00
    ! of the first day of the given month.  The last entry corresponds
    ! to 00:00:00 of the first day of the next month.
    !
    ! E.g.; for Tair_cru198207.nc, the data run from 1982-07-01T03:00:00
    ! through 1982-08-01T00:00:00, inclusive.
    !
    ! So, when you are at hour 0 on the first day of the month,
    ! reset the day to the last day of the previous month, and reset the
    ! hour from 0 to 24.  This will compute the correct index for the last
    ! entry in the forcing file.
    !
    ! E.g.; 1982-08-01T00:00:00 should be re-written as 1982-07-31T24:00:00.
    if ( tmp_da == 1 .and. tmp_hr == 0 .and. mn == 0 .and. ss == 0) then
       if ( mo == 1 ) then
          tmp_da = days1(12)
       else
          tmp_da = days1(mo-1)
       endif
       tmp_hr = 24
    endif

    index = (tmp_da-1)*8 + (tmp_hr+mn/60.0+ss/3600)/3
    !index = index + 1

  end subroutine getgswp_timeindex
end module gswp_module
