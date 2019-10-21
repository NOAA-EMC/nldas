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
! !ROUTINE: vic_writerst.F90
!
! !DESCRIPTION:
!  This program writes restart files for VIC.  This
!   includes all relevant water/energy storages, tile information,
!   and time information. 
!
!
! !REVISION HISTORY:
! 
! 17 Nov 2003; Justin Sheffield : Initial Version
! 21 Nov 2003; Sujay Kumar : Added the time manager restart. 
! 
! !INTERFACE:
#include "misc.h"
subroutine vic_writerestart
! !USES:
  use lisdrv_module, only: lis,tile
  use vic_varder, only : vicdrv
  use time_manager
  use spmdMod
!EOP
  character*80 :: arch_file
  character*80 :: fdir
  character(len=40) :: cdate       !date char string
  character(len=40) :: cyear       !year char string
  character(len=40) :: cnhdate     !no hour date char string
  character(len=40) :: ccode
  integer :: day                    !day (1 -> 31)
  integer :: mon                    !month (1 -> 12)
  integer :: yr                     !year (0 -> ...)
  integer :: sec, testsec           !seconds into current day
!BOC
  if ( ( lis%t%gmt == ( 24 - vicdrv%writeintvic ) ) &
       .or. lis%t%endtime == 1 ) then

     call gather_model_state(vicdrv%vic_snowband, &
                             vicdrv%vic_nlayer,   &
                             vicdrv%vic_nnode)

     if ( masterproc ) then
!       call get_curr_date (yr, mon, day, sec) 
       yr  = lis%t%yr
       mon = lis%t%mo
       day = lis%t%da
       sec = lis%t%ss

       testsec = sec / 3600

       write(cdate,'(i4.4,i2.2,i2.2,i2.2)') yr,mon,day,testsec
       write(cnhdate,'(i4.4,i2.2,i2.2)') yr,mon,day
       write(cyear,'(i4)') yr
       write(ccode,'(i3)') lis%o%expcode

       fdir = "mkdir -p "//trim(lis%o%odir)//"/EXP"//trim(adjustl(ccode))//&
              "/VIC/"//trim(cyear)//"/"//trim(cnhdate)

       call system(trim(fdir))

       arch_file = trim(lis%o%odir)//"/EXP"//trim(adjustl(ccode))//"/VIC/"//&
                   trim(cyear)//"/"//trim(cnhdate)//"/LIS.E"//              &
                   trim(adjustl(ccode))//"."//trim(cdate)//".VICrst"

!<kluge>
! We no longer write time info into the restart file
!        OPEN(40,FILE="time.rst",FORM='unformatted') 
!        call timemgr_write_restart(40)
!</kluge>

!<kluge>
! We currently need all process of a GDS-based (OPENDAP) run to call
! the write_model_state routine.  This is to work around a problem
! with single variable output mode.  In general only the master process
! writes data.  See write_model_state routine.
#if ( defined OPENDAP ) && ( defined SPMD )
endif
#endif
!</kluge>
        call write_model_state(trim(vicdrv%VIC_RFILE), & 
             len(trim(vicdrv%VIC_RFILE)), &
             lis%d%glbnch,         &
             vicdrv%vic_snowband, vicdrv%vic_nlayer, vicdrv%vic_nnode)
        call write_model_state(trim(arch_file), & 
             len(trim(arch_file)), &
             lis%d%glbnch,         &
             vicdrv%vic_snowband, vicdrv%vic_nlayer, vicdrv%vic_nnode)
!<kluge>
#if ( defined OPENDAP ) && ( defined SPMD )
if ( masterproc ) then
#endif
!</kluge>

!<kluge>
! We no longer write time info into the restart file
!        close(40)
!</kluge>
     endif
  endif
!EOC
end subroutine vic_writerestart
