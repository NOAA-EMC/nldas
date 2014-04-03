!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! glbobspcp.F90: 
!
! DESCRIPTION:
!  Opens and reads global precipitation forcing
!
!    CTIME  = Current time
!    FTIMENRL  = Nearest future data for NRL data
!    FTIMEHUFF = Nearest future data for HUFFMAN data
!    FTIMEPERS = Nearest future data for PERSIANN data
!
! REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 10 Oct 2001: Jon Gottschalck; Modified to adjust convective precip
!               using a ratio of the model convective / total ratio
! 30 Jul 2002: Jon Gottschalck; Added PERSIANN and HUFFMAN global observed precip data sources
!=========================================================================

subroutine glbobspcp(ldas, grid)

  use ldas_module      ! LDAS non-model-specific 1-D variables
  use grid_module      ! LDAS non-model-specific grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc, ldas%nr)	
	
!==== Local Variables=======================
  integer :: c, r, flag1, flag2                                     ! Program flow flags used in HUFFMAN precip section
  integer :: ferror_nrl, ferror_huff, ferror_pers, ferror_cmap      ! Error flags for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1                ! Time parameters for current LDAS time
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2                ! Time parameters for NRL boundary end time
  integer :: doy3, yr3, mo3, da3, hr3, mn3, ss3, ts3                ! Time parameters for HUFFMAN boundary end time
  integer :: kdoy3, kyr3, kmo3, kda3, khr3, kmn3, kss3, kts3        ! Time parameters for datatime (used for HUFFMAN data, see below)
  integer :: mdoy3, myr3, mmo3, mda3, mhr3, mmn3, mss3, mts3        ! Time parameters for fnametime (used for HUFFMAN data, see below)
  integer :: doy4, yr4, mo4, da4, hr4, mn4, ss4, ts4                ! Time parameters for PERSIANN boundary end time
  integer :: doy5, yr5, mo5, da5, hr5, mn5, ss5, ts5                ! Time parameters for CMAP boundary end time
  integer :: endtime_nrl,endtime_pers,endtime_huff,endtime_cmap     ! 1=get a new file 
  real*8  :: ctime,ftime_nrl,ftime_huff,ftime_pers,ftime_cmap       ! Current LDAS time and end boundary times for precip data sources 
  real*8  :: datatime, gap, breaktime, fnametime                    ! Times used in HUFFMAN to determine data and filename boundaries (see below)
  real    :: gmt1, gmt2, gmt3, gmt4, kgmt3, mgmt3, gmt5             ! GMT times for current LDAS time and end boundary times for precip data sources
  real    :: ratio(ldas%nc, ldas%nr)                                ! Ratio of model convective precip to total model precip
  character(len=80) :: name,name_trmm,name_huff,name_pers,name_cmap ! Filename variables for precip data sources

!=== End Variable Definition =======================

!=== Assumption will be made to not find any data
  endtime_nrl  = 0
  endtime_pers = 0
  endtime_huff = 0
  endtime_cmap = 0

!=== Set parameter to measure 1.5 hour time offset when using HUFFMAN
  gap = 0.0001712328767098370

!=== Determine required observed precip data times (current, accumulation end time)
!=== Model current time
  yr1 = ldas%yr  !current time
  mo1 = ldas%mo
  da1 = ldas%da
  hr1 = ldas%hr
  mn1 = ldas%mn
  ss1 = 0
  ts1 = 0
  call tick( ctime, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

!=== NRL products end time
  yr2 = ldas%yr  !end accumulation time data
  mo2 = ldas%mo
  da2 = ldas%da
  hr2 = 6*(ldas%hr/6)
  mn2 = 0
  ss2 = 0
  ts2 = 6*60*60
  call tick( ftime_nrl, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

88 format(f14.8,1x,f14.8,1x,f14.8,1x,i3,1x,f14.8,1x,i2)
89 format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2,1x,i6)	

!=== HUFFMAN product end time
  yr3 = ldas%yr  !end accumulation time data
  mo3 = ldas%mo
  da3 = ldas%da
  hr3 = 3*(ldas%hr/3)
  mn3 = 0
  ss3 = 0
  ts3 = 3*60*60
  call tick( ftime_huff, doy3, gmt3, yr3, mo3, da3, hr3, mn3, ss3, ts3 )

!=== Determine times for HUFFMAN program flow
!=== breaktime: used to determine where in 1.5 hour time offset current time falls
!=== datatime:  used as the end boundary time for using the data
!=== fnametime: used to get proper filename
  breaktime = ftime_huff - ctime
  datatime  = ftime_huff
  fnametime = ftime_huff

!=== Need to adjust above end boundary time for HUFFMAN to account for a 1.5 hour time offset forward in time
!=== Need to addjust for time when data is valid and for the filename itself (these are different)
!=== The first if/endif determines which 1.5 half of the 3 hour interval it is.
!=== If first half of 3 hour interval must set times back 1.5 and 3 hours. If the time is 0 hours GMT
!=== must go back these values plus 24 hours due to a peculiarity in time.f 
!=== For the second half, simply move forward three 3 hours.
!=== When in second half of the 3 hour interval must have a special condition for first time step for proper program flow.
!=== The flag variables are used for program flow from one half of the 3 hour interval to the other.

  if (ldas%gpcpsrc(3) == 1) then

!  write(*,*) "BREAK/GAP TIME = ",breaktime, gap
  if (breaktime .ge. gap) then
!    write(*,*) "MGMT1.0 = ", mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3
!     write(*,*) "KGMT1.0 = ", kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3
    call time2date( datatime, kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3 )
    call time2date( fnametime, mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3 )
    flag1 = 1
!    print*, "section 1 ",kgmt3, mgmt3, flag1, flag2
!     write(*,*) "MGMT1.1 = ", mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3
!     write(*,*) "KGMT1.1 = ", kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3
     if (khr3 == 24) khr3 = 0
     if (mhr3 == 24) mhr3 = 0
     if (kgmt3 .eq. 0.0 .and. flag2 .eq. 2) then
      kts3 = -25.5*60*60
      call tick( datatime, kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3, kss3, kts3 )
      mts3 = -27*60*60
      call tick( fnametime, mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3, mss3, mts3 )
     else
      kts3 = -1.5*60*60
      call tick( datatime, kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3, kss3, kts3 )
      mts3 = -3*60*60
      call tick( fnametime, mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3, mss3, mts3 )
     endif
     flag2 = 1
!     write(*,*) "MGMT1.2 = ", mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3
!     write(*,*) "KGMT1.2 = ", kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3
  else
    if (ldas%tscount .eq. 1) then
     call time2date( datatime, kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3 )
     call time2date( fnametime, mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3 )
!     write(*,*) "3KGMT3.3 = ", kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3
!     write(*,*) "4MGMT3.3 = ", mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3
     if (kgmt3 .eq. 0) then
      mts3 = -24*60*60
      call tick( fnametime, mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3, mss3, mts3 )
      kts3 = -22.5*60*60
      call tick( datatime, kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3, kss3, kts3 )
     else
      mts3 = 0
      call tick( fnametime, mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3, mss3, mts3 )
      kts3 = 1.5*60*60
      call tick( datatime, kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3, kss3, kts3 )
     endif
    else
     flag1 = 2
!     print*, "section 2 ", kgmt3, mgmt3, flag1, flag2
!     write(*,*) "MGMT2.1 = ", mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3
!     write(*,*) "KGMT2.1 = ", kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3
      if (flag2 .eq. 1) then
       mts3 = 3*60*60
       call tick( fnametime, mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3, mss3, mts3 )
       kts3 = 3*60*60
       call tick( datatime, kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3, kss3, kts3 )
      endif
      flag2 = 2
!     print*, "section 2 ", kgmt3, mgmt3, flag1, flag2
!     write(*,*) "MGMT2.2 = ", mdoy3, mgmt3, myr3, mmo3, mda3, mhr3, mmn3
!     write(*,*) "KGMT2.2 = ", kdoy3, kgmt3, kyr3, kmo3, kda3, khr3, kmn3
    endif
  endif

!  write(*,*) "CTIME = ",ctime," flag1/2 = ",flag1, flag2
   
   endif
  
!=== PERSIANN product end time  
  yr4 = ldas%yr  !end accumulation time data
  mo4 = ldas%mo
  da4 = ldas%da
  hr4 = 1*(ldas%hr/1)
  mn4 = 0
  ss4 = 0
  ts4 = 1*60*60
  call tick( ftime_pers, doy4, gmt4, yr4, mo4, da4, hr4, mn4, ss4, ts4 )
  
!=== CMAP product end time
  yr5 = ldas%yr  !end accumulation time data
  mo5 = ldas%mo
  da5 = ldas%da
  hr5 = 6*(ldas%hr/6)
  mn5 = 0
  ss5 = 0
  ts5 = 6*60*60
  call tick( ftime_cmap, doy5, gmt5, yr5, mo5, da5, hr5, mn5, ss5, ts5 )
  
!=== Compute ratio between convective model precip and total model precip
!=== so that it can be applied to the observed global precip
     do c = 1,ldas%nc
      do r = 1,ldas%nr
       if (grid(c,r)%forcing(8) .ne. 0.0 .and.       &
&          grid(c,r)%forcing(8) .ne. ldas%udef .and. &
&          grid(c,r)%forcing(9) .ne. ldas%udef) then
         ratio(c,r) = grid(c,r)%forcing(9) / grid(c,r)%forcing(8) 
         if (ratio(c,r) .gt. 1.0) ratio(c,r) = 1.0
         if (ratio(c,r) .lt. 0.0) ratio(c,r) = 0.0
       else
         ratio(c,r) = 0.0
       endif
      enddo
     enddo

!=== Ensure that data is found during first time step

  if ( ldas%gpcpsrc(1) .gt. 0 .and. ldas%tscount == 1 ) endtime_nrl = 1
  if ( ldas%gpcpsrc(2) .gt. 0 .and. ldas%tscount == 1 ) endtime_pers = 1
  if ( ldas%gpcpsrc(3) .gt. 0 .and. ldas%tscount == 1 ) endtime_huff = 1
  if ( ldas%gpcpsrc(4) .gt. 0 .and. ldas%tscount == 1 ) endtime_cmap = 1

!=== Check to see if current time has crossed boundary end time

!=== Check for and get NRL observed Precipitation data
  if (ldas%gpcpsrc(1) .gt. 0) then
   if ( ctime > ldas%nrltime ) endtime_nrl = 1
    if ( endtime_nrl == 1 ) then  !get new time2 data
       print*, 'Getting new NRL satellite precip data'
       ferror_nrl = 0
       call nrlfile( name, ldas%nrldir, yr2, mo2, da2, hr2 )
       call nrlfile_trmm( name_trmm, ldas%nrldir, yr2, mo2, da2, hr2 )
       call glbprecip_nrl( ldas, grid, name, name_trmm, ferror_nrl )
       ldas%nrltime = ftime_nrl
    endif  !need new time2
  endif

!=== Overriding current forcing precipiattion data with latest NRL data if relevant
  if (ferror_nrl /= 0 .and. ldas%gpcpsrc(1) .gt. 0) then
   do c = 1, ldas%nc
     do r = 1, ldas%nr
       if (grid(c,r)%obsprecip .ne. -1.0) then
	 grid(c,r)%forcing(8) = grid(c,r)%obsprecip / 3600.0
	 grid(c,r)%forcing(9) = ratio(c,r) * grid(c,r)%forcing(8)
       endif
      enddo
    enddo	    

 endif

!=== Check for and get PERSIANN observed Precipitation data
  if (ldas%gpcpsrc(2) .gt. 0) then
!===Added if/endif since PERSIANN is hourly data and if used > for timestep 1 hr it would skip every other hour of precip
   if (ldas%ts < 3600) then
     if ( ctime > ldas%perstime ) endtime_pers = 1
   else
     if ( ctime >= ldas%perstime ) endtime_pers = 1
   endif
     if ( endtime_pers == 1 ) then  !get new time2 data
      print*, 'Getting new PERSIANN satellite precip data'
      ferror_pers = 0
      call persfile( name_pers, ldas%persdir, yr4, mo4, da4, hr4 )
      call glbprecip_pers( ldas, grid, name_pers, ferror_pers )
      ldas%perstime = ftime_pers
     endif  !need new time2
  endif

!=== Overriding current forcing precipiattion data with latest PERSIANN data if relevant
  if (ferror_pers /= 0 .and. ldas%gpcpsrc(2) .gt. 0) then
   do c = 1, ldas%nc
     do r = 1, ldas%nr
       if (grid(c,r)%obsprecip .ne. -1.0) then
	  grid(c,r)%forcing(8) = grid(c,r)%obsprecip / 3600.0
	  grid(c,r)%forcing(9) = ratio(c,r) * grid(c,r)%forcing(8)
	endif
      enddo
    enddo
  endif

!=== Check for and get HUFFMAN observed Precipitation data
  if (ldas%gpcpsrc(3) .gt. 0) then
   if ( ctime > ldas%hufftime ) then      
     endtime_huff = 1
     if ( endtime_huff == 1 ) then  !get new time2 data
      print*, 'Getting new HUFFMAN satellite precip data', endtime_huff
      ferror_huff = 0
      call hufffile( name_huff, ldas%huffdir, myr3, mmo3, mda3, mhr3 )
      call glbprecip_huff( ldas, grid, name_huff, ferror_huff )
       ldas%hufftime = datatime
       endif
     endif  !need new time2
  endif		      
  
!=== Overriding current forcing precipiattion data with latest HUFFMAN data if relevant
  if (ferror_huff /= 0 .and. ldas%gpcpsrc(3) .gt. 0) then
     do c = 1, ldas%nc
      do r = 1, ldas%nr
       if (grid(c,r)%obsprecip .ne. -1.0) then
	   grid(c,r)%forcing(8) = grid(c,r)%obsprecip / 3600.0
	   grid(c,r)%forcing(9) = ratio(c,r) * grid(c,r)%forcing(8)
       endif					       
      enddo
     enddo
     
  endif
  
!=== Check for and get CMAP CPC Precipitation data
  if (ldas%gpcpsrc(4) .gt. 0) then
    if ( ctime > ldas%cmaptime ) endtime_cmap = 1
      if ( endtime_cmap == 1 ) then  !get new time2 data
        print*, 'Getting new CMAP CPC precip data'
	ferror_cmap = 0
	call cmapfile( name, ldas%cmapdir, yr5, mo5, da5, hr5 )
        call glbprecip_cmap( ldas, grid, name, ferror_cmap, hr5 )
	ldas%cmaptime = ftime_cmap
      endif  !need new time2
    endif			   

!=== Overriding current forcing precipitation data with latest CMAP data if relevant
  if (ferror_cmap /= 0 .and. ldas%gpcpsrc(4) .gt. 0) then
     do c = 1, ldas%nc
       do r = 1, ldas%nr
	if (grid(c,r)%obsprecip .ne. -1.0) then
	  grid(c,r)%forcing(8) = grid(c,r)%obsprecip / 3600.0
	  grid(c,r)%forcing(9) = ratio(c,r) * grid(c,r)%forcing(8)
	endif
       enddo
     enddo
   endif

  return

end subroutine glbobspcp

!==============================================================
!
!  DESCRIPTION: This subroutine puts together NRL file name for
!               6 hour file intervals
!==============================================================

subroutine nrlfile( name, nrldir, yr, mo, da, hr)

  implicit none

!==== Local Variables=======================

  character(len=80) :: name, nrldir
  integer :: yr, mo, da, hr
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: doy
  real    :: gmt
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(8)

!=== End Variable Definition ===============

!=== formats for filename segments
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, i2)
95 format (10a1)
96 format (a40)
97 format (a8)
98 format (a1, i4, i2, a1)
99 format (8a1)

!=== Make variables for the time used to create the file
!=== We don't want these variables being passed out
  uyr = yr
  umo = mo
  uda = da
  uhr = 6*(hr/6)  !hour needs to be a multiple of 6 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.

  open(unit=90, file='temp', form='formatted', access='direct', recl=80)
  write(90, 96, rec=1) nrldir 
  read(90, 92, rec=1) (fbase(i), i=1,80)

  write(90, 98, rec=1) '/', uyr, umo, '/'
  read(90, 99, rec=1) fdir
  do i = 1, 8
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(90, 94, rec=1) uyr, umo, uda, uhr
  read(90, 95, rec=1) ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  write(90, 97, rec=1) '00.6.geo'
  read (90, 92, rec=1) (fsubs(i), i=1,8)

!sets c as the last character position of fbase
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(90, 92, rec=1) (fbase(i), i=1,c), (fdir(i), i=1,8),  &
                       (ftime(i), i=1,10), (fsubs(i), i=1,8)

  read(90, 93, rec=1) name
  
  close(90)
  return

end subroutine nrlfile

!==============================================================
!
!  DESCRIPTION: This subroutine puts together NRL file name for
!               6 hour file intervals (SSMITRMM)
!==============================================================

  subroutine nrlfile_trmm( name, nrldir, yr, mo, da, hr)

  implicit none
  
!==== Local Variables=======================
  
  character(len=80) :: name, nrldir
  integer :: yr, mo, da, hr
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: doy
  real    :: gmt
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(13)
			      
!=== End Variable Definition ===============
			      
!=== formats for filename segments
  92 format (80a1)
  93 format (a80)
  94 format (i4, i2, i2, i2)
  95 format (10a1)
  96 format (a40)
  97 format (a13)
  98 format (a1, i4, i2, a1)
  99 format (8a1)
				      
!=== Make variables for the time used to create the file
!=== We don't want these variables being passed out
      uyr = yr
      umo = mo
      uda = da
      uhr = 6*(hr/6)  !hour needs to be a multiple of 6 hours
      umn = 0
      uss = 0
      ts1 = -24*60*60 !one day interval to roll back date.

      open(unit=90, file='temp', form='formatted',access='direct', recl=80)
      write(90, 96, rec=1) nrldir
      read(90, 92, rec=1) (fbase(i), i=1,80)
		     
      write(90, 98, rec=1) '/', uyr, umo, '/'
      read(90, 99, rec=1) fdir
      do i = 1, 8
       if ( fdir(i) == ' ' ) fdir(i) = '0'
      end do
      
     write(90, 94, rec=1) uyr, umo, uda, uhr
     read(90, 95, rec=1) ftime
     do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
     end do
				    
     write(90, 97, rec=1) '00.6.ssmitrmm'
     read (90, 92, rec=1) (fsubs(i), i=1,13)

!sets c as the last character position of fbase
     c = 0
     do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
     end do
			 
     write(90, 92, rec=1) (fbase(i), i=1,c),(fdir(i),i=1,8),  &
     (ftime(i), i=1,10), (fsubs(i), i=1,13)
     
     read(90, 93, rec=1) name
	  
     close(90)
     return
		    
   end subroutine nrlfile_trmm
   
!==============================================================
!
!  DESCRIPTION: This subroutine puts together HUFFMAN file name for
!               3 hour file intervals
!==============================================================

subroutine hufffile( name, huffdir, yr, mo, da, hr)

  implicit none
  
!==== Local Variables=======================

  character(len=80) :: name, huffdir 
  integer :: yr, mo, da, hr
  integer :: i, c 
  integer :: uyr, umo, uda, uhr, umn, uss, ts1 
  integer :: doy 
  real :: gmt 
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(4)
  character*1 :: fprefix(7)

!=== End Variable Definition ===============
!=== formats for filename segments
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, i2)
95 format (10a1)
96 format (a40)
97 format (a4)
98 format (a1, i4, i2, a1)
99 format (8a1)
89 format (a7)

!=== Make variables for the time used to create the file
!=== We don't want these variables being passed out

  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.
  
  open(unit=90, file='temp', form='formatted',access='direct', recl=80)
  write(90, 96, rec=1) huffdir
  read(90, 92, rec=1) (fbase(i), i=1,80)
  
  write(90, 98, rec=1) '/', uyr, umo, '/'
  read(90, 99, rec=1) fdir
  do i = 1, 8
   if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do
			   
  write(90, 94, rec=1) uyr, umo, uda, uhr
  read(90, 95, rec=1) ftime
  do i = 1, 10
    if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do
		
  write(90, 97, rec=1) '.bin'
  read (90, 92, rec=1) (fsubs(i), i=1,4)
  
  
  write(90, 89, rec=1) '3B42RT.'
  read (90, 92, rec=1) (fprefix(i), i=1,7)
     
!sets c as the last character position of fbase
  c = 0
  do i = 1, 80
   if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do
			 			 
  write(90, 92, rec=1) (fbase(i),i=1,c),(fdir(i),i=1,8),(fprefix(i), i=1,7),  &
    (ftime(i), i=1,10), (fsubs(i), i=1,4)
     
  read(90, 93, rec=1) name
	  
  close(90)
  return
		    
end subroutine hufffile

!==============================================================
!
!  DESCRIPTION: This subroutine puts together PERSIANN file name for
!               1 hour file intervals
!==============================================================

subroutine persfile( name, persdir, yr, mo, da, hr)

  implicit none
  
!==== Local Variables=======================
  
  character(len=80) :: name, persdir
  integer :: yr, mo, da, hr
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: doy
  real :: gmt
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(4)

!=== End Variable Definition ===============
!=== formats for filename segments
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, i2)
95 format (10a1)
96 format (a40)
97 format (a4)
98 format (a1, i4, i2, a1)
99 format (8a1)
89 format (a7)

!=== Make variables for the time used to create the file
!=== We don't want these variables being passed out

  uyr = yr
  umo = mo
  uda = da
  uhr = 1*(hr/1)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.

  open(unit=90, file='temp', form='formatted',access='direct', recl=80)
  write(90, 96, rec=1) persdir
  read(90, 92, rec=1) (fbase(i), i=1,80)
  
  write(90, 98, rec=1) '/', uyr, umo, '/'
  read(90, 99, rec=1) fdir
  do i = 1, 8
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do
  
  write(90, 94, rec=1) uyr, umo, uda, uhr
  read(90, 95, rec=1) ftime
  do i = 1, 10
   if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do
  
  write(90, 97, rec=1) '.bin'
  read (90, 92, rec=1) (fsubs(i), i=1,4)
  
!===sets c as the last character position of fbase
  c = 0
  do i = 1, 80
    if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do
	   
  write(90, 92, rec=1) (fbase(i),i=1,c),(fdir(i),i=1,8),  &
                       (ftime(i), i=1,10), (fsubs(i), i=1,4)
    
  read(90, 93, rec=1) name
      
  close(90)
    
  return
  
  end subroutine persfile

!==============================================================
!
!  DESCRIPTION: This subroutine puts together CMAP file name for
!               6 hour file intervals
!==============================================================

subroutine cmapfile( name, cmapdir, yr, mo, da, hr)

  implicit none

!==== Local Variables=======================

  character(len=80) :: name, cmapdir
  integer :: yr, mo, da, hr
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: doy
  real    :: gmt
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(8), fsubs2(4)

!=== End Variable Definition ===============

!=== formats for filename segments
91 format (a4)
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, i2)
95 format (10a1)
96 format (a40)
97 format (a10)
98 format (a1, i4, i2, a1)
99 format (8a1)

!=== Make variables for the time used to create the file
!=== We don't want these variables being passed out
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

!sets c as the last character position of fbase
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

end subroutine cmapfile

