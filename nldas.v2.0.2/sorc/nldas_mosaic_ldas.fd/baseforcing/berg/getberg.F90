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
! !ROUTINE: getberg.F90
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 6-hrly, 1/2 degree Reanalysis 
!  ECMWF forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data 
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! !REVISION HISTORY:
!  11 Apr 2002: Urszula Jambor; original code based on getgeos.f
!  22 Oct 2002: Urszula Jambor; Limited SW forcing processing to 
!               land-only grid points
!  24 Nov 2003: Sujay Kumar; Included ECMWF code in LIS
! !INTERFACE:
subroutine getberg()
! !USES:
  use lisdrv_module, only : lis    ! LIS non-model-specific 1-D variables
  use baseforcing_module, only : glbdata1,glbdata2
  use bergdomain_module, only : bergdrv
  use time_manager
!EOP
  implicit none
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer, parameter :: nforce=9  ! # forcing variables
  integer :: try, ferror
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1,ts1,bdoy,byr,bmo
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2,ts2,bda,bhr,bmn
  real*8  :: time1,time2,dumbtime1,dumbtime2
  real*8  :: timenow
  real    :: gmt1,gmt2
  integer :: movetime      ! 1=move time 2 data into time 1
	 
  !=== Assumption will be not to find or move any data
  lis%f%findtime1=0
  lis%f%findtime2=0
  movetime=0
  !=== Determine Required Data Times (The previous hour & the future hour)

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
  hr1=6*((lis%t%hr)/6)
  mn1=0
  ss1=0
  ts1=0
  call tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2=lis%t%yr    !Next Hour
  mo2=lis%t%mo
  da2=lis%t%da
  hr2=6*((lis%t%hr)/6)
  mn2=0
  ss2=0
  ts2=6*60*60

  call tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)	 

  if(timenow.gt.bergdrv%fmodeltime2) then
     movetime=1
     lis%f%findtime2=1
  endif

  if(get_nstep(lis%t).eq.0 .or. get_nstep(lis%t).eq.1 & 
       .or. lis%f%rstflag.eq.1) then  !beginning of the run
     lis%f%findtime1=1
     lis%f%findtime2=1
     glbdata1 = 0
     glbdata2 = 0
     movetime=0
     lis%f%rstflag = 0
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
        call retberg(order,yr1,mo1,da1,hr1,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           bergdrv%fmodeltime1=time1
        else  !ferror still=0, so roll back one day
           call tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        end if
        if ( try > ndays ) then 
           print *, 'ERROR: berg data gap exceeds 10 days on file 1'
           STOP
        end if
     end do
  endif
  
  !  Repeat for time 2
  
  if(movetime.eq.1) then !transfer time2 data to time1
     bergdrv%fmodeltime1=bergdrv%fmodeltime2	
     lis%f%findtime2=1 !include to ensure getting new time2 data
     
     do f=1,nforce
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
        call retberg(order,yr2,mo2,da2,hr2,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           bergdrv%fmodeltime2=time2
        else  !ferror still=0, so roll back one day
           call tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        end if
        if ( try > ndays ) then 
           print *, 'ERROR: berg data gap exceeds 10 days on file 2'
           STOP
        end if
     end do
  endif

end subroutine getberg
	   






