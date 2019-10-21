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
! !ROUTINE: getnldas.F90
!
! !DESCRIPTION:
!  Opens, reads, and interpolates NCEP-LDAS forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!  The strategy for missing data is to backwards up to 10 days to get
!  forcing at the same time of day.
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  20 Dec 1999: Paul Houser; Allow for Eta Data Overwrite by NCEP data
!  27 Apr 2000: Brian Cosgrove; Turned zenith angle weighting back on.
!               Changed times supplied to ZTERP from GMT1 and GMT2 to 
!               LDAS%NLDASTIME1 and LDAS%NLDASTIME2
!   4 May 2000: Added 15 minutes to GMT1 and GMT2 to accurately
!               reflect the valid times of the NCEP radiation data
!  18 Aug 2000: Brian Cosgrove; Fixed error in date calculations so that
!               3600 second (1hr) timestep may be used.
!               Added code to make any calculated radiation forcing
!               values undefined if both ncep time 1 and ncep time 2
!               radiation values are undefined
!  27 Feb 2001: Brian Cosgrove; Added CZM into call for ZTERP subroutine
!  07 Mar 2001: Brian Cosgrove; Added code to allow for use of NASA-LDAS data 
!  04 Sep 2001: Brian Cosgrove; Changed tempgmt1,tempgmt2 to real to match
!               tick.f call, changed file name construction.
!  21 Aug 2002: Brian Cosgrove; Removed code that adjusted for 15 minute
!               offset in radiation fields supplied by NOAA or NASA
!               NLDAS realtime or retrospective standard forcing files.
!               This offset no longer exists as it is now dealt with 
!               during forcing file creation through zenith angle correction.
!               The need for this fix was just discovered...unfortunately
!               simulations before the date of this fix using 
!               this fortran subroutine have incorrectly shifted
!               radiation data.
! 02Feb 2004 : Sujay Kumar ; Initial Version in LIS
! 
! !INTERFACE:
subroutine getnldas()
! !USES:
  use lisdrv_module, only : lis,grid
  use baseforcing_module, only : glbdata1,glbdata2
  use nldasdomain_module, only : nldasdrv
  use time_manager, only : tick
!EOP  
  implicit none
!=== Local Variables =====================================================
  integer :: c,r,f,ferror,try,zdoy
  real*8  :: time1,time2,timenow
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1,ts1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2,ts2
  character*80 :: name
  real :: gmt1,gmt2
  integer:: movetime     ! 1=move time 2 data into time 1  

!=== End Variable Definition =============================================
  try=-999

!====Assumption will be not to find or move any data
  lis%f%findtime1=0
  lis%f%findtime2=0
  movetime=0
  
!=== Determine Required NCEP Data Times (The previous hour and the future hour)
  yr1 = lis%t%yr
  mo1=lis%t%mo
  da1=lis%t%da
  hr1=lis%t%hr
  mn1=lis%t%mn
  ss1=0
  ts1=0
  call tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr1 = lis%t%yr
  mo1=lis%t%mo
  da1=lis%t%da
  hr1=lis%t%hr
  mn1=0
  ss1=0
  ts1=0

  call tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2=lis%t%yr    !next hour
  mo2=lis%t%mo
  da2=lis%t%da
  hr2=lis%t%hr
  mn2=0
  ss2=0
  ts2=60*60
        !print *,'-1 yr2,mo2,da2,hr2',yr2,mo2,da2,hr2
  call tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        !print *,'0 yr2,mo2,da2,hr2',yr2,mo2,da2,hr2


  if(timenow.gt.nldasdrv%nldastime2) then 
     movetime = 1
     lis%f%findtime2 = 1
  endif
!  print*, timenow, nldasdrv%nldastime2, lis%f%findtime1, lis%f%findtime2

  if(lis%t%tscount.eq.1 .or.lis%f%rstflag.eq.1  ) then    !beginning of the run	
	print *,'at beginning of run'
      lis%f%findtime1=1
      lis%f%findtime2=1
      movetime=0
      lis%f%rstflag = 0
   endif
   
   if(movetime.eq.1) then
	print *,'moving time2 into time1'
      nldasdrv%nldastime1=nldasdrv%nldastime2
      do f=1,lis%f%nf
         do c=1,lis%d%ngrid
            glbdata1(f,c)=glbdata2(f,c)
         enddo
      enddo
   endif    !end of movetime=1
   
   if(lis%f%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
      ferror=0
      try=0  
      ts1=-60*60*24
      do 
         if ( ferror /= 0 ) exit
         try=try+1
         call ncepfile(name,nldasdrv%nldasdir,yr1,mo1,da1,hr1)
        print *,'yr1,mo1,da1,hr1',yr1,mo1,da1,hr1
         print *,'getting file1.. ',name
         call retnldas(1,name,ferror,1,0)
         if(ferror.eq.1) nldasdrv%nldastime1=time1
         call tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
         if(try.gt.11)then
            write(*,*)'error: ncep data gap exceeds 10 days on file 1'
            stop
         endif
      enddo
!=== end of data search
   endif   !end of lis%f%findtime=1	   	


   if(lis%f%findtime2.eq.1) then  
!=== the following looks back 10 days, at the same hour to fill data gaps.
      ferror=0
      try=0  
      ts2=-60*60*24
      do 
         if ( ferror /= 0 ) exit
         try=try+1
        !print *,'1 yr2,mo2,da2,hr2',yr2,mo2,da2,hr2
         call ncepfile(name,nldasdrv%nldasdir,yr2,mo2,da2,hr2)
	!print *,'2 yr2,mo2,da2,hr2',yr2,mo2,da2,hr2
         print *,'getting file2.. ',name
         call retnldas(2,name,ferror,1,0)
        !print *,'3 time2,yr2,mo2,da2,hr2',time2,yr2,mo2,da2,hr2,'ferror=',ferror
         if(ferror.eq.1)nldasdrv%nldastime2=time2
         if(ferror.eq.0)call tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        !print *,'4 time2,yr2,mo2,da2,hr2',time2,yr2,mo2,da2,hr2
         if(try.gt.11)then
            write(*,*)'error: ncep data gap exceeds 10 days on file 2'
            stop
         endif
      enddo
!=== end of data search
   endif   ! end of findtime2=1a
   return
 end subroutine getnldas

!BOP
! !ROUTINE: ncepfile.f: 
!
! !DESCRIPTION:
!  This subroutine puts together the ncep data filename
!
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
! !INTERFACE:
 subroutine ncepfile(name,ncepdir,yr,mo,da,hr)
!EOP
   implicit none

   character*80, intent(out)::  name
   character*40, intent(in) ::  ncepdir
   integer, intent(in)      :: yr,mo,da,hr

   integer                  :: i, c
   character*1              :: fname(80),fbase(80),fsubs(80)
   character*1              ::  ftime(10),fdir(15)

   !=== end variable definition =============================================

   !=== put together filename
91 format(a4,i3,a11,i3)
92 format(80a1)
93 format(a80)
94 format(i4,i2,i2,i2)
95 format(10a1)
96 format(a40)
97 format(a14) 
89 format(a14)
98 format(a1,i4,a1,i4,i2,i2,a1)
99 format(15a1)
   open(90,file='temp',form='formatted',access='direct',recl=80)

   write(90,98,rec=1)'/',yr,'/',yr,mo,da,'/'
   read(90,99,rec=1)fdir
   do i=1,15
      if(fdir(i).eq.(' '))fdir(i)='0'
   enddo

   write(90,94,rec=1)yr,mo,da,hr
   read(90,95,rec=1)ftime
   do i=1,10
      if(ftime(i).eq.(' '))ftime(i)='0'
   enddo

!   write(90,97,rec=1)'.FORCING.GRB'
   write(90,97,rec=1)'.lsmforce_noaa'
   read(90,92,rec=1) (fsubs(i),i=1,14)

   write(90,96,rec=1) ncepdir                       
   read(90,92,rec=1) (fbase(i),i=1,80)

   c=0
   do i=1,80
      if(fbase(i).eq.(' ').and.c.eq.0)c=i-1
   enddo

   write(90,92,rec=1)(fbase(i),i=1,c), (fdir(i),i=1,15), & 
        (ftime(i),i=1,10),(fsubs(i),i=1,14 ) 
   read(90,93,rec=1)name
   close(90)

   return
 end subroutine ncepfile






































