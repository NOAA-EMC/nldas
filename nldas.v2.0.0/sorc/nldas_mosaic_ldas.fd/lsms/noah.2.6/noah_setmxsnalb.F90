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
! !ROUTINE: set_mxsnalb
!
! !DESCRIPTION:
!  This subroutine retrieves NOAH max snow albebo.
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added NOAH LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!  24 Jun 2005: James Geiger; Merged in Sujay Kumar's GSWP-2 support
!
! !INTERFACE:
subroutine noah_setmxalb
! !USES:
   use lisdrv_module, only : tile,lis
   use noah_varder      ! NOAH tile variables
   use lis_openfileMod
   use lis_indices_module
!EOP      
   implicit none

!=== Local Variables =====================================================

   integer           :: i           !loop counters
   real, allocatable :: tmpalb(:,:)
   real, allocatable :: tmpalb1d(:)
!BOC  

!-----------------------------------------------------------------------
! The MAX SNOW ALBEDO field is opened and read in here:
!-----------------------------------------------------------------------

   if ( noahdrv%mxsnalb_type == 1 ) then

      allocate(tmpalb(lis_nc_data,lis_nr_data))

      call lis_open_file(12,file=noahdrv%noah_mxsnal,        & 
                         access='direct',form="unformatted", &
                         recl=4, script='getmaxsnalb.pl')
  
      call lis_read_file(12,tmpalb)

      close(12)

      do i = 1, lis%d%nch
         if ( tmpalb(tile(i)%col, tile(i)%row-lis_tnroffset) /= -9999.00 ) then
            noah(i)%mxsnalb = tmpalb(tile(i)%col, tile(i)%row-lis_tnroffset)
         endif
      enddo
  
      deallocate(tmpalb)

   elseif ( noahdrv%mxsnalb_type == 2 ) then

      allocate(tmpalb1d(lis%d%nch))

      call lis_log_msg('MSG: noah_setmxsnalb -- Opening GSWP MAXSNALB File: '//&
                       trim(noahdrv%noah_mxsnal))

      open(unit=12,file=noahdrv%noah_mxsnal, form='unformatted')     

      read(12) tmpalb1d

      close(12)

      do i = 1, lis%d%nch
         noah(i)%mxsnalb = tmpalb1d(i) / 100.0
      enddo
  
      deallocate(tmpalb1d)

   else

      call lis_log_msg("ERR: noah_setmxsnalb -- " //          &
                       "Don't know how to read mxsnalb.  " // &
                       "Please check the definition of noahdrv%mxsnalb_type.")
      call endrun

   endif

!EOC
end subroutine noah_setmxalb

