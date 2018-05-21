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
! !ROUTINE:  getecmwf.F90: 
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 10 LDAS forcing fields from operational 
!  ECMWF model output.  
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data 
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!  Forcing fields are read in depending on data type.  Time integrated 
!  data are treated differently than instantaneous data. 
!
! !REVISION HISTORY:
!  18 Jun 2003: Urszula Jambor; original code based on getreanlecmwf.F90
! 
! !INTERFACE:
subroutine getecmwf()
! !USES:
  use lisdrv_module, only : lis, gindex     
  use baseforcing_module, only : glbdata1, glbdata2   
  use time_manager
  use ecmwfdomain_module, only : ecmwfdrv
!EOP
  implicit none
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer :: try, ferror
  integer :: c,r,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1,ts1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2,ts2
  real*8  :: time1,time2,dumbtime1,dumbtime2
  real    :: gmt1,gmt2
  integer :: movetime      ! 1=move time 2 data into time 1
  integer :: nforce  ! # forcing variables
!=== End Variable Definition =======================

  !=== Assumption will be not to find or move any data
  lis%f%findtime1=0
  lis%f%findtime2=0
  movetime=0

  if(get_nstep(lis%t).eq.0 ) then
    nforce = ecmwfdrv%nmif
  else
    nforce = lis%f%nf
  endif
  if (get_nstep(lis%t) .eq. 1.or.lis%f%rstflag.eq.1) then
    lis%f%findtime1=1
    lis%f%findtime2=1
    glbdata1 = 0
    glbdata2 = 0
    movetime=0
    lis%f%rstflag = 0
  endif
					    
  !=== Determine Required Data Times (The previous hour & the future hour)
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

  !=== Check if time interval boundary was crossed
  if(lis%t%time.gt.ecmwfdrv%ecmwftime2) then
     movetime=1
     lis%f%findtime2=1
  endif

  !=== Establish fmodeltime1
  if (lis%f%findtime1==1) then  !need to get new time1 from the past
     order=1   !Get data for glbdata1
     ferror = 0
     try = 0
     ts1 = -24*60*60
     do
        if ( ferror /= 0 ) then
           exit
        end if
        try = try+1
        call retecmwf(order,yr1,mo1,da1,hr1,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           ecmwfdrv%ecmwftime1=time1
        else  !ferror still=0, so roll back one day
           call tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        end if
        if ( try > ndays ) then 
           print *, 'ERROR: ECMWF data gap exceeds 10 days'
           STOP
        end if
     end do
  endif

  !=== Find new time2 value and tranfer time2 data to time1 
  if(movetime.eq.1) then
     ecmwfdrv%ecmwftime1=ecmwfdrv%ecmwftime2
     lis%f%findtime2=1 !to ensure getting new time2 data
     do f=1,lis%f%nforce
        do c=1,lis%d%ngrid
          glbdata1(f,c)=glbdata2(f,c)
        enddo
     enddo  
  endif  ! if movetime=1
  
  if(lis%f%findtime2.eq.1) then ! need new time2 data
     order=2   !Get data for glbdata2
     ferror = 0
     try = 0
     ts2 = -24*60*60
     do
        if ( ferror /= 0 ) exit
        try = try+1
        call retecmwf(order,yr2,mo2,da2,hr2,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           print*, 'reset ecmwftime2 to time2'
           ecmwfdrv%ecmwftime2=time2
        else  !ferror still=0, so roll back one day
           call tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        end if
        if ( try > ndays ) then
           print *, 'ERROR: ECMWF data gap exceeds 10 days'
           STOP
        end if
     end do
  endif
  
end subroutine getecmwf
   






