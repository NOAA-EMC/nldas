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
! !ROUTINE: getcmap.F90
!
! !DESCRIPTION:
!  Opens and reads global precipitation forcing
!
!    CTIME  = Current time\\
!    FTIMENRL  = Nearest future data for NRL data\\
!    FTIMEHUFF = Nearest future data for HUFFMAN data\\
!    FTIMEPERS = Nearest future data for PERSIANN data\\
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 10 Oct 2001: Jon Gottschalck; Modified to adjust convective precip
!               using a ratio of the model convective / total ratio
! 30 Jul 2002: Jon Gottschalck; Added PERSIANN and HUFFMAN global observed precip data sources
!
! !INTERFACE:
subroutine getcmap()
! !USES:
  use lisdrv_module, only : lis, gindex
  use time_manager
  use cmapdomain_module, only : cmapdrv, conserv_cmap_interp_input
  implicit none
!EOP
   
!==== Local Variables=======================
  integer :: c, r, flag1, flag2                                     ! Program flow flags used in HUFFMAN precip section
  integer :: ferror_cmap,ferror_huff, ferror_pers      ! Error flags for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1                ! Time parameters for current LDAS time
  integer :: doy5, yr5, mo5, da5, hr5, mn5, ss5, ts5                ! Time parameters for CMAP boundary end time
  integer :: endtime_cmap   ! 1=get a new file 
  real*8  :: ctime,ftime_cmap  ! Current LDAS time and end boundary times for precip data sources 
  real*8  :: datatime, gap, breaktime, fnametime                    ! Times used in HUFFMAN to determine data and filename boundaries (see below)
  real    :: gmt1, gmt2, gmt3, gmt4, kgmt3, mgmt3, gmt5             ! GMT times for current LDAS time and end boundary times for precip data sources
  character(len=80) :: name,name_cmap ! Filename variables for precip data sources
  integer :: index
  real :: gridDesci(50)
!=== End Variable Definition =======================
!BOC
  endtime_cmap = 0
!------------------------------------------------------------------------
! Set parameter to measure 1.5 hour time offset when using HUFFMAN
!------------------------------------------------------------------------
  gap = 0.0001712328767098370
!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end time)
! Model current time
!------------------------------------------------------------------------
  yr1 = lis%t%yr  !current time
  mo1 = lis%t%mo
  da1 = lis%t%da
  hr1 = lis%t%hr
  mn1 = lis%t%mn
  ss1 = 0
  ts1 = 0
  call tick( ctime, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   
!------------------------------------------------------------------------ 
! CMAP product end time
!------------------------------------------------------------------------
  yr5 = lis%t%yr  !end accumulation time data
  mo5 = lis%t%mo
  da5 = lis%t%da
  hr5 = 6*(lis%t%hr/6)
  mn5 = 0
  ss5 = 0
  ts5 = 6*60*60
  call tick( ftime_cmap, doy5, gmt5, yr5, mo5, da5, hr5, mn5, ss5, ts5 )

  if ( ctime > cmapdrv%griduptime1 .and. & 
       ctime < cmapdrv%griduptime2 .and. & 
       cmapdrv%gridchange1 ) then 

     call lis_log_msg('MSG: getcmap -- changing cmap grid to 2002-2005')

     cmapdrv%ncold = 768
     cmapdrv%nrold = 384
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 768
     gridDesci(3) = 384
     gridDesci(4) = 89.462
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.462
     gridDesci(8) = -0.469
     gridDesci(9) = 0.469
     gridDesci(10) = 192
     gridDesci(20) = 255

     call conserv_cmap_interp_input(gridDesci,lis%d%gridDesc,&
                                    lis%d%lnc*lis%d%lnr)
     cmapdrv%gridchange1 = .false.

  elseif ( ctime > cmapdrv%griduptime2 .and. & 
           cmapdrv%gridchange2 ) then 

     call lis_log_msg('MSG: getcmap -- changing cmap grid to 2005-')

     cmapdrv%ncold = 1152
     cmapdrv%nrold = 576
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1152
     gridDesci(3) = 576
     gridDesci(4) = 89.761
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.761
     gridDesci(8) = -0.313
     gridDesci(9) = 0.313
     gridDesci(10) = 288
     gridDesci(20) = 255

     call conserv_cmap_interp_input(gridDesci,lis%d%gridDesc,&
                                    lis%d%lnc*lis%d%lnr)
     cmapdrv%gridchange2 = .false.

  endif
!------------------------------------------------------------------------
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  if ( lis%f%gpcpsrc.eq.4.and. get_nstep(lis%t).eq. 1 & 
  .or.lis%f%rstflag .eq. 1) then 
     endtime_cmap = 1
     lis%f%rstflag = 0
  endif
!------------------------------------------------------------------------
! Check for and get CMAP CPC Precipitation data
!------------------------------------------------------------------------
  if (lis%f%gpcpsrc==4) then
     if ( ctime > cmapdrv%cmaptime ) endtime_cmap = 1
     if ( endtime_cmap == 1 ) then  !get new time2 data
        ferror_cmap = 0
        call cmapfile( name, cmapdrv%cmapdir, yr5, mo5, da5, hr5 )
        print*, 'Getting new CMAP CPC precip data',name
        call glbprecip_cmap( name, ferror_cmap, hr5 )
        cmapdrv%cmaptime = ftime_cmap
     endif  !need new time2
  endif
  return
!EOC
end subroutine getcmap

!BOP
! !ROUTINE: cmapfile
!
! !DESCRIPTION: This subroutine puts together CMAP file name for
!               6 hour file intervals
!
! !INTERFACE:
subroutine cmapfile( name, cmapdir, yr, mo, da, hr)
!EOP
  implicit none

!==== Local Variables=======================

  character(len=80) :: name, cmapdir
  integer :: yr, mo, da, hr
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: doy
  real    :: gmt
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(10), fsubs2(4)

!=== End Variable Definition ===============

!=== formats for filename segments
!BOC
91 format (a4)
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, i2)
95 format (10a1)
96 format (a40)
97 format (a10)
98 format (a1, i4, i2, a1)
99 format (8a1)
!------------------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!------------------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 6*(hr/6)  !hour needs to be a multiple of 6 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.

  open(unit=90, file='temp', form='formatted', access='direct', recl=80)
  write(90, 96, rec=1) cmapdir
  read(90, 92, rec=1) (fbase(i), i=1,80)

  write(90, 98, rec=1) '/', uyr, umo, '/'
  read(90, 99, rec=1) fdir
  do i = 1, 8
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(90, 97, rec=1) 'cmap_gdas_'
  read (90, 92, rec=1) (fsubs(i), i=1,10)

  write(90, 94, rec=1) uyr, umo, uda, uhr
  read(90, 95, rec=1) ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  write(90, 94, rec=1) uyr, umo, uda, uhr
  read(90, 95, rec=1) ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  write(90, 91, rec=1) '.grb'
  read (90, 92, rec=1) (fsubs2(i), i=1,4)
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(90, 92, rec=1) (fbase(i), i=1,c), (fdir(i), i=1,8),  &
                       (fsubs(i), i=1,10),(ftime(i), i=1,10), &
                       (fsubs2(i), i=1,4)

  read(90, 93, rec=1) name

  close(90)
  return
!EOC
end subroutine cmapfile

