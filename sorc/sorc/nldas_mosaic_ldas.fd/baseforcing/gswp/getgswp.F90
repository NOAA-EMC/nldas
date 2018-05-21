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
! !ROUTINE: getgswp.F90
!  
! !DESCRIPTION: 
!
!  Opens, reads, and interpolates GSWP forcing.  
!
!    TIME1 = most recent past data\\
!    TIME2 = nearest future data \\
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! \subsection{Core Functions of getgswp}
!  \begin{description}
!  \item[tick]  
!      Determines GSWP data times
!  \item[gswpfile]
!      Puts together appropriate file name for 3 hour intervals
!  \item[readgswp]
!      Interpolates GSWP data to LDAS grid
!  \end{description}
!
! !REVISION HISTORY:
!
! 20Feb2004; Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine getgswp()
! !USES:
  use lisdrv_module, only : lis
  use time_manager
  use spmdMod
  use tile_spmdMod
  use baseforcing_module, only: glbdata1,glbdata2,glbdata3
  use gswpdomain_module, only : gswpdrv
!EOP
  implicit none
!==== Local Variables=======================
  integer :: try, ferror
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1,ts1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2,ts2
  integer :: yr3,mo3,da3,hr3,mn3,ss3,doy3,ts3
  real*8 :: time1,time2,time3,dumbtime1,dumbtime2
  real*8 :: timenow
  character*80 :: name
  character*40 :: elevfile, fpart1, fpart2
  real :: gmt1,gmt2,gmt3
  integer :: movetime      ! 1=move time 2 data into time 1
  integer :: nforce     ! GSWP forcing file time, # forcing variables
  integer :: nstep,ierr
  integer :: gridDesci(50)
!BOC  
  if ( masterproc ) then
     nstep = get_nstep(lis%t)
  endif

!-------------------------------------------------------------------
! Determine the correct number of forcing variables
!-------------------------------------------------------------------
  if ( nstep == 0 ) then
     nforce = gswpdrv%nmif
  else
     nforce = lis%f%nf
  endif
  lis%f%findtime1=0
  lis%f%findtime2=0
  lis%f%shortflag = 2
  lis%f%longflag=2             !Time averaged LW 
  movetime=0
!-------------------------------------------------------------------
! Determine Required GSWP Data Times 
! (The previous hour & the future hour)
!-------------------------------------------------------------------
  yr1=lis%t%yr    !Time now
  mo1=lis%t%mo
  da1=lis%t%da
  hr1=lis%t%hr
  mn1=lis%t%mn
  ss1=0
  ts1=0        
  
  call tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
      
  yr1=lis%t%yr    !Previous Hour
  mo1=lis%t%mo
  da1=lis%t%da
  hr1=3*((lis%t%hr)/3)
  mn1=0
  ss1=0
  ts1=0
  call tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr2=lis%t%yr    !Next Hour
  mo2=lis%t%mo
  da2=lis%t%da
  hr2=3*((lis%t%hr)/3)
  mn2=0
  ss2=0
  ts2=3*60*60
  call tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

  yr3=lis%t%yr    !Uber-Next Hour
  mo3=lis%t%mo
  da3=lis%t%da
  hr3=3*((lis%t%hr)/3)
  mn3=0
  ss3=0
  ts3=2*3*60*60
  call tick(time3,doy3,gmt3,yr3,mo3,da3,hr3,mn3,ss3,ts3)

  ! For GSWP, we must read past (previous hour), current (next hour), and
  ! future (uber-next hour) forcing time-steps.
  !
  ! So, for GSWP, lis%f%findtime1 means update past and current forcing data.
  !
  ! lis%f%findtime2 means update the future forcing data.
  ! 
  ! gswpdrv%gswptime1 is the update time for the past forcing time-step,
  ! which is only needed at initialization.
  !
  ! gswpdrv%gswptime2 is the update time for the current forcing time-step.
  !
  ! When it is time to update the current forcing time-step, we actually shift
  ! the old current to past, the old future to current, and we read in
  ! the new future.
  if ( timenow > gswpdrv%gswptime2 ) then
     movetime        = 1
     lis%f%findtime2 = 1
  endif
  
  if ( nstep == 0 .or. nstep == 1 .or. lis%f%rstflag == 1 ) then 
     lis%f%findtime1 = 1
     lis%f%findtime2 = 1
     glbdata1        = 0
     glbdata2        = 0
     glbdata3        = 0
     movetime        = 0
     lis%f%rstflag   = 0
  endif

  if ( lis%f%findtime1 == 1 ) then 
     call lis_log_msg('MSG: getgswp -- reading time1 data')
     print*,'getgswp',yr1,mo1,da1,hr1,mn1,ss1
     order = 1
     call readgswp(order,yr1,mo1,da1,hr1,mn1,ss1)
     gswpdrv%gswptime1 = time1

     call lis_log_msg('MSG: getgswp -- reading time2 data')
     print*,'getgswp',yr2,mo2,da2,hr2,mn2,ss2
     order = 2   
     call readgswp(order,yr2,mo2,da2,hr2,mn2,ss2)
     gswpdrv%gswptime2 = time2
  endif

  if ( movetime == 1 ) then
     gswpdrv%gswptime1=gswpdrv%gswptime2
     lis%f%findtime2=1 
     do f=1,nforce
        do c=1,lis%d%ngrid
           glbdata1(f,c)=glbdata2(f,c)
           glbdata2(f,c)=glbdata3(f,c)
        enddo
     enddo
  endif 

  if ( lis%f%findtime2 == 1 ) then 
     call lis_log_msg('MSG: getgswp -- reading time3 data')
     print*,'getgswp',yr3,mo3,da3,hr3,mn3,ss3
     order = 3   
     call readgswp(order,yr3,mo3,da3,hr3,mn3,ss3)
     gswpdrv%gswptime2 = time2
  endif
  
end subroutine getgswp
     
