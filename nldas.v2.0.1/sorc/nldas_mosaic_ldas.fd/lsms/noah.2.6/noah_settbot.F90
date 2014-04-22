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
! !ROUTINE: setnoahp.F90
!
! !DESCRIPTION:
!  This subroutine retrieves NOAH bottom temperature.
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added NOAH LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!  24 Jun 2005: James Geiger; Merged in Sujay Kumar's GSWP-2 support
!
! !INTERFACE:
subroutine noah_settbot
! !USES:
   use lisdrv_module, only : tile,lis
   use noah_varder      ! NOAH tile variables
   use lis_openfileMod
   use lis_indices_module
!EOP      
   implicit none

   integer           :: i                  !loop counters
   real, allocatable :: placetbot(:,:)
   real, allocatable :: placetbot1d(:)
!BOC  

!-----------------------------------------------------------------------
! Read in bottom temperature fields and adjust for elevation differnce
! with either Eta (NLDAS) or GDAS (GLDAS) correction datasets. 
!-----------------------------------------------------------------------

   if ( noahdrv%tbot_type == 1 ) then

      allocate(placetbot(lis_nc_data,lis_nr_data))

      call lis_open_file(12, file=noahdrv%noah_tbot,   & 
                         access='direct',status='old', &
                         form='unformatted', recl=4, script='gettbot.pl')

      call lis_read_file(12,placetbot)

      close(12)
  
      call lis_log_msg('MSG: noah_settbot -- Read TBOT file')

      do i = 1, lis%d%nch
         if ( placetbot(tile(i)%col,tile(i)%row-lis_tnroffset) /= -9999.0 ) then
            noah(i)%tempbot = placetbot(tile(i)%col, tile(i)%row-lis_tnroffset)
         endif
      enddo

      deallocate(placetbot)

   elseif ( noahdrv%tbot_type == 2 ) then

      call lis_log_msg('MSG: noah_settbot -- Opening GSWP TBOT File: '// &
                       trim(noahdrv%noah_tbot))

      allocate(placetbot1d(lis%d%nch))

      open(unit=12, file=noahdrv%noah_tbot, form='unformatted')
      read(12) placetbot1d
      close(12)

      call lis_log_msg('MSG: noah_settbot -- Read TBOT file')

      do i = 1, lis%d%nch
         noah(i)%tempbot = placetbot1d(i)
      enddo

      close(12)

      deallocate(placetbot1d)

   else

      call lis_log_msg("ERR: noah_settbot -- Don't know how to read tbot.  "//&
                       "Please check the definition of noahdrv%tbot_type.")
      call endrun

    endif

end subroutine noah_settbot

