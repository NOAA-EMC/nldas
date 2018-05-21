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
! !ROUTINE: getgdas.F90
!  
! !DESCRIPTION: 
! 
!  Opens, reads, and interpolates NCEP-GDAS forcing.  
!
!    TIME1 = most recent past data\\
!    TIME2 = nearest future data \\
!
!  The idea is to open either the 00 or 03 forecast file associated with
!  the most recent GDAS assimilation (available every 6 hours).
!  Precipitation rates and radiation fluxes will be taken from the F03 and 
!  F06 files, since averages are provided.\\
!  - if that fails, the strategy for missing data is to go backwards up to
!    10 days to get forcing at the same time of day.
!
! \subsection{Core Functions of getgdas}
!  \begin{description}
!  \item[tick]  
!      Determines GDAS data times
!  \item[gdasfile]
!      Puts together appropriate file name for 3 hour intervals
!  \item[gdasfilef06]
!      Puts together appropriate file name for 6 hour intervals
!  \item[retgdas]
!      Interpolates GDAS data to LDAS grid
!  \end{description}
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  25 Oct 1999: Jared Entin; Significant F90 Revision
!  11 Apr 2000: Brian Cosgrove; Fixed name construction error 
!               in Subroutine ETA6HRFILE 
!  27 Apr 2000: Brian Cosgrove; Added correction for use of old shortwave
!               data with opposite sign convention from recent shortwave data.
!               Added capability to use time averaged shortwave and longwave data.
!               Altered times which are passed into ZTERP--used to be GMT1 and GMT2,
!               now they are LDAS%ETATIME1 and LDAS%ETATIME2
!  11 May 2000: Brian Cosgrove; Added checks for SW values that are too high
!               due to zenith angle weighting...if too high, use linear weighting.
!               Also, if cos(zen) less than .1, then use linear weighting to
!               avoid computed values of SW that are too high.
!  18 May 2000: Brian Cosgrove; Corrected line of code in ETAEDASNAME which
!               assigned wrong year directory variable when constructing
!               EDAS filename
!  26 May 2000: Jared Entin; Changed numerical bound of the TRY variable
!               to fix a rollback problem.
!  5 June 2000: Brian Cosgrove; Fixed a problem with the correction of the negative
!               radiation sign convention. Prior to fix, was not correcting negative
!               values of -999.9...now it changes all negative values to positive ones.
!  18 Aug 2000: Brian Cosgrove; Fixed undefined value problem in check for negative
!               radiation values over land points.
!  08 Dec 2000: Urszula Jambor; Rewrote geteta.f in fortran90 to use GDAS in GLDAS
!  15 Mar 2001: Jon Gottschalck; Slight change to handle more forcing parameters at
!               time step 0.
!  09 Apr 2001: Urszula Jambor; Added capability of using DAAC forcing data every
!               6 hours, rather than every 3 hours.
!  30 May 2001: Urszula Jambor; Changed forcing used: T,q,u fields taken 
!               from F00 & F03 files, radiation and precip. fields taken 
!               from F06 & F03 (F03 fields are subtracted out from F06)
! !INTERFACE:
#include "misc.h"
subroutine getgdas()
! !USES:
  use lisdrv_module, only : lis, gindex      
  use baseforcing_module, only: glbdata1, glbdata2
  use time_manager
  use gdasdomain_module, only : gdasdrv
  use bilinear_interpMod, only : bilinear_interp_input
  use conserv_interpMod, only : conserv_interp_input
  use spmdMod
  use lis_indices_module, only: lis_nc_working, lis_nr_working
!EOP
  implicit none
  integer :: c, r, f,i
  integer :: ferror
  integer :: try
  integer, parameter :: ndays = 10  ! # of days to look back for forcing data
  integer :: order     ! 1 indicates lesser interpolation boundary time
                       ! 2 indicates greater interpolation boundary time
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2
  integer :: zdoy
  integer :: movetime  ! 1=move time2 into time1
  integer :: nforce    ! Number of forcing variables for model init. option
  integer :: tindex
  real*8  :: timenow, time1, time2
  real*8  :: dumbtime1, dumbtime2
  real    :: gmt1, gmt2
  character(len=80) :: name, nameF06="null"
  real :: gridDesci(50)
  integer :: nstep

!BOC
  lis%f%F00_flag = 0 
  lis%f%F06_flag = 0 

  lis%f%findtime1=0
  lis%f%findtime2=0
  movetime=0

!-----------------------------------------------------------------
! Determine the correct number of forcing variables
!-----------------------------------------------------------------
  if ( masterproc ) then
     nstep = get_nstep(lis%t)
  endif
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
  call MPI_BCAST(nstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ferror)
#endif
  if ( nstep == 0 ) then 
     nforce = gdasdrv%nmif
  else
     nforce = lis%f%nf
  endif
  if ( nstep == 1 .or. lis%f%rstflag == 1) then
    lis%f%findtime1=1
    lis%f%findtime2=1
    glbdata1 = 0
    glbdata2 = 0
    movetime=0        ! movetime is not properly set at time-step = 1
    lis%f%rstflag = 0
  endif 


!-----------------------------------------------------------------
! Determine required GDAS data times 
! (previous assimilation, current & future assimilation hours)
! The adjustment of the hour and the direction will be done
! in the subroutines that generate the names
!-----------------------------------------------------------------
  yr1 = lis%t%yr  !current time
  mo1 = lis%t%mo
  da1 = lis%t%da
  hr1 = lis%t%hr
  mn1 = lis%t%mn
  ss1 = 0
  ts1 = 0
  call tick( timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

  yr1 = lis%t%yr  !previous assimilation/forecast hour
  mo1 = lis%t%mo
  da1 = lis%t%da
  if ( lis%f%F06_flag == 0 ) then
     hr1 = 3*(int(real(lis%t%hr)/3.0))
  else
     hr1 = 6*(int(real(lis%t%hr)/6.))
  end if
  mn1 = 0
  ss1 = 0
  ts1 = 0
  call tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

  yr2 = lis%t%yr  !next assimilation/forecast hour
  mo2 = lis%t%mo
  da2 = lis%t%da
  if ( lis%f%F06_flag == 0 ) then
     hr2 = 3*(int(real(lis%t%hr)/3.0))
  else
     hr2 = 6*(int(real(lis%t%hr)/6.0))
  end if
  mn2 = 0
  ss2 = 0
  if ( lis%f%F06_flag == 0 ) then
     ts2 = 3*60*60
  else
     ts2 = 6*60*60
  end if
  call tick( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )
!-----------------------------------------------------------------
! Use these if need to roll back time.
!-----------------------------------------------------------------
  dumbtime1 = time1
  dumbtime2 = time2
!-----------------------------------------------------------------
! Check to see if current time (timenow) has crossed past gdastime2,
! requiring that both gdastime2 be assigned to gdastime1 and a new
! gdastime2 be set 3 or 6 hours ahead of the current gdastime2.
!-----------------------------------------------------------------
  if ( timenow > gdasdrv%gdastime2 ) then
     movetime  = 1
     lis%f%findtime2 = 1
  end if
  
  if ( time1 > gdasdrv%griduptime1 .and. &
       time1 < gdasdrv%griduptime2 .and. & 
       gdasdrv%gridchange1 ) then 

     call lis_log_msg('MGS: getgdas -- changing grid to 2000-2002')
     
     gdasdrv%ncold = 512
     gdasdrv%nrold = 256
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 512
     gridDesci(3) = 256
     gridDesci(4) = 89.463
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.463
     gridDesci(8) = -0.703
     gridDesci(9) = 0.703
     gridDesci(10) = 128
     gridDesci(20) = 255

     if ( lis%f%interp == 1 ) then 
        call bilinear_interp_input(gridDesci,lis%d%gridDesc,&
             lis_nc_working*lis_nr_working)
     elseif ( lis%f%interp == 2 ) then 
        call bilinear_interp_input(gridDesci,lis%d%gridDesc,&
             lis_nc_working*lis_nr_working)
        call conserv_interp_input(gridDesci,lis%d%gridDesc,&
             lis_nc_working*lis_nr_working)
     endif

     gdasdrv%gridchange1 = .false.

  elseif ( time1 > gdasdrv%griduptime2 .and. & 
           time1 < gdasdrv%griduptime3 .and. &
           gdasdrv%gridchange2 ) then 

     call lis_log_msg('MGS: getgdas -- changing grid to 2002-2005')

     gdasdrv%ncold = 768
     gdasdrv%nrold = 384
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 768
     gridDesci(3) = 384
     gridDesci(4) = 89.642
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.642
     gridDesci(8) = -0.469
     gridDesci(9) = 0.469
     gridDesci(10) = 192
     gridDesci(20) = 255

     if ( lis%f%interp == 1 ) then 
        call bilinear_interp_input(gridDesci,lis%d%gridDesc, &
                                   lis_nc_working*lis_nr_working)
     elseif ( lis%f%interp == 2 ) then 
        call bilinear_interp_input(gridDesci,lis%d%gridDesc, &
                                   lis_nc_working*lis_nr_working)
        call conserv_interp_input(gridDesci,lis%d%gridDesc, &
                                  lis_nc_working*lis_nr_working)
     endif

     gdasdrv%gridchange2 = .false.

     if ( lis%f%ecor == 1 ) then 
        call update_gdas_elevdiff("2")
     endif
  elseif ( time1 > gdasdrv%griduptime3 .and. & 
           gdasdrv%gridchange3 ) then 

     call lis_log_msg('MGS: getgdas -- changing grid to 2005-')

     gdasdrv%ncold = 1152
     gdasdrv%nrold = 576
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

     if ( lis%f%interp == 1 ) then 
        call bilinear_interp_input(gridDesci,lis%d%gridDesc, &
                                   lis_nc_working*lis_nr_working)
     elseif ( lis%f%interp == 2 ) then 
        call bilinear_interp_input(gridDesci,lis%d%gridDesc, &
                                   lis_nc_working*lis_nr_working)
        call conserv_interp_input(gridDesci,lis%d%gridDesc, &
                                  lis_nc_working*lis_nr_working)
     endif

     gdasdrv%gridchange3 = .false.

     if ( lis%f%ecor == 1 ) then 
        call update_gdas_elevdiff("3")
     endif
  endif
!-----------------------------------------------------------------
! Establish gdastime1
!-----------------------------------------------------------------
  if ( lis%f%findtime1 == 1 ) then  !get new time1 from the past
     print *, 'Getting new time1 data'
     ferror = 0
     order = 1
     try = 0
     ts1 = -24*60*60
     
     do
        if ( ferror /= 0 ) exit
        try = try+1
        if ( lis%f%F06_flag == 0 ) then
           call gdasfile   ( name, gdasdrv%gdasdir, yr1, mo1, da1, hr1 )
        else
           call gdasfileF06( name, gdasdrv%gdasdir, yr1, mo1, da1, hr1 )
        end if
        print*, 'Reading GDAS file1 ',name
        call retgdas( order, name, nameF06, 0,ferror, try)
        
        if ( ferror == 1 ) then  
!-----------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------
           gdasdrv%gdastime1 = time1
        else  
!-----------------------------------------------------------------
! ferror still=0, so roll back one day & start again
!-----------------------------------------------------------------
           call tick( dumbtime1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
        end if
        if ( try > ndays ) then  
!-----------------------------------------------------------------
! data gap exceeds 10 days so stop
!-----------------------------------------------------------------
           print *, 'ERROR: GDAS data gap exceeds 10 days on file 1'
           call endrun
        end if
     end do  
  end if  
!-----------------------------------------------------------------
! Establish gdastime2
!-----------------------------------------------------------------
  if ( movetime == 1 ) then  
     gdasdrv%gdastime1 = gdasdrv%gdastime2
     lis%f%findtime2 = 1  
     do f = 1, nforce
        do c = 1, lis%d%ngrid
           glbdata1(f,c) = glbdata2(f,c)
        end do
     end do
  end if
  
  if ( lis%f%findtime2 == 1 ) then  
     print *, 'Getting new time2 data'
     lis%f%F00_flag = 0 
  
     ferror = 0
     order = 2
     try = 0
     ts2 = -24*60*60
!-----------------------------------------------------------------
! determine if both F00 & F06 files needed
!-----------------------------------------------------------------
     if ( lis%f%F06_flag == 0 ) then
        if ( modulo(hr2,2) > 0 ) then
           lis%f%F00_flag = 0 !odd hr, use F03 file
        else
           lis%f%F00_flag = 1 !even hr, use F00 & F06
        end if
     end if
     do
        if ( ferror /= 0 ) exit
        try = try+1
        if ( lis%f%F06_flag == 0 ) then
           if ( lis%f%F00_flag == 0 ) then
              call gdasfile( name,gdasdrv%gdasdir, yr2, mo2, da2, hr2 )
           else
              call gdasfile( name,gdasdrv%gdasdir, yr2, mo2, da2, hr2 )
              call gdasfileF06( nameF06, gdasdrv%gdasdir, yr2, mo2, da2, hr2 )
           end if
        else
           call gdasfileF06( name, gdasdrv%gdasdir, yr2, mo2, da2, hr2 )
        end if
        print*, 'Reading GDAS file2 ',name
        call retgdas( order, name, nameF06, &
             lis%f%F00_flag, ferror,try )
        if ( ferror == 1 ) then  
!-----------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------
           gdasdrv%gdastime2 = time2
        else  
!-----------------------------------------------------------------
! ferror still=0, so roll back one day & start again
!-----------------------------------------------------------------
           call tick( dumbtime2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )
        end if
        if ( try > ndays ) then  
!-----------------------------------------------------------------
! data gap exceeds 10 days so stop
!-----------------------------------------------------------------
           print *, 'ERROR: GDAS data gap exceeds 10 days on file 2'
           call endrun
        end if
     end do  
  end if
!-----------------------------------------------------------------
! loop through all forcing parameters & interpolate
!-----------------------------------------------------------------
  do f = 1, lis%f%nforce  
     if ( (f == 3) .or. (f == 4) ) then
        do c = 1, lis_nc_working
           do r = 1, lis_nr_working
              tindex = gindex(c,r)
              if(tindex .ne.-1) then
                 if ( (glbdata2(f,tindex) /= -9999.9) .and.  &
                      (glbdata2(f,tindex) < 0)) then
                    glbdata2(f,tindex) = (-1) * glbdata2(f,tindex)
                 end if
                 if ( (glbdata1(f,tindex) /= -9999.9) .and.  &
                      (glbdata1(f,tindex) < 0)) then
                    glbdata1(f,tindex) = (-1) * glbdata1(f,tindex)
                 end if
              end if
           end do  
        end do  
     end if  
  enddo
!EOC  
end subroutine getgdas


!BOP
! !ROUTINE: gdasfile
!
! !DESCRIPTION:
!   This subroutine puts together GDAS file name for 
!   3 hour file intervals
!
! !INTERFACE:
subroutine gdasfile( name, gdasdir, yr, mo, da, hr )

  implicit none

! !INPUT PARAMETERS:
  character(len=80), intent(in)    :: gdasdir
  integer, intent(in)              :: yr, mo, da, hr
! !OUTPUT PARAMETERS:
  character(len=80), intent (out)  :: name
!EOP
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: remainder
  integer :: doy
  real    :: gmt
  real*8  :: dumbtime
  character(len=2) :: initcode, fcstcode
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(21)
  character(LEN=100) :: temp
!=== End Variable Definition ===============

!=== formats for filename segments
!BOC
!-----------------------------------------------------------------
!  Make variables for the time used to create the file
!  We don't want these variables being passed out
!-----------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.
  remainder = modulo(uhr,2)  !if even, then remainder equals zero
                             !if odd, then remainder equals one
!-----------------------------------------------------------------
! if hour is   00 or 03, look for 00ZF00 or 00ZF03 
! if hour is   06 or 09, look for 06ZF00 or 06ZF03
! if hour is   12 or 15, look for 12ZF00 or 12ZF03
! if hour is   18 or 21, look for 18ZF00 or 18ZF03
!-----------------------------------------------------------------
  if ( uhr <= 3 ) then
     initcode = '00'
  else if (uhr <= 9 ) then
     initcode = '06'
  else if ( uhr <= 15 ) then
     initcode = '12'
  else if ( uhr <= 21 ) then
     initcode = '18'
  end if
  if ( remainder > 0 ) then
     fcstcode = '03'
  else
     fcstcode = '00'
  end if
  
  write(UNIT=temp, fmt='(a40)') gdasdir  
  read(UNIT=temp, fmt='(80a1)') (fbase(i), i=1,80)

  write(UNIT=temp, fmt='(a1, i4, i2, a1)') '/', uyr, umo, '/'
  read(UNIT=temp, fmt='(8a1)') fdir
  do i = 1, 8
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(UNIT=temp, fmt='(i4, i2, i2, a2)') uyr, umo, uda, initcode
  read(UNIT=temp, fmt='(10a1)') ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  write(UNIT=temp, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode, '.sg'
  read (UNIT=temp, fmt='(80a1)') (fsubs(i), i=1,21)

  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(UNIT=temp, fmt='(80a1)') (fbase(i), i=1,c), (fdir(i), i=1,8),  &
                       (ftime(i), i=1,10), (fsubs(i), i=1,21)

  read(UNIT=temp, fmt='(a80)') name

  return
!EOC
end subroutine gdasfile



!BOP
! !ROUTINE: gdasfilef06 
!
! !DESCRIPTION:
!   This subroutine puts together GDAS file name for 
!   6 hour forecast files (DAAC archive friendly)
!
! !INTERFACE:
subroutine gdasfileF06( name, gdasdir, yr, mo, da, hr )
! !USES:
  use time_manager

  implicit none

! !INPUT PARAMETERS:
  character(len=80) :: gdasdir
  integer :: yr, mo, da, hr
! !OUTPUT PARAMETERS:
  character(len=80) :: name
!EOP
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: doy
  real    :: gmt
  real*8    :: dumbtime
  character(len=2) :: initcode, fcstcode
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(21)
  character(LEN=100) :: temp
!=== End Variable Definition ===============

!BOC
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, a2)
95 format (10a1)
96 format (a40)
97 format (a16, a2, a3)
98 format (a1, i4, i2, a1)
99 format (8a1)
!-----------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!-----------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 6*(hr/6)  !hour needs to be a multiple of 6 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.
!-----------------------------------------------------------------
! if hour is 00Z look for 18ZF06 after rolling back one day
! if hour is 06Z look for 00ZF06
! if hour is 12Z look for 06ZF06
! if hour is 18Z look for 12ZF06
!-----------------------------------------------------------------
  fcstcode = '06'
  if ( uhr == 0 ) then
     call tick( dumbtime, doy, gmt, uyr, umo, uda, uhr, umn, uss, ts1 )
     initcode = '18'
  else if ( uhr == 6 ) then
     initcode = '00'
  else if ( uhr == 12 ) then
     initcode = '06'
  else if ( uhr == 18 ) then
     initcode = '12'
  end if

  write(UNIT=temp, fmt='(a40)') gdasdir  
  read(UNIT=temp, fmt='(80a1)') (fbase(i), i=1,80)

  write(UNIT=temp, fmt='(a1, i4, i2, a1)') '/', uyr, umo, '/'
  read(UNIT=temp, fmt='(8a1)') fdir
  do i = 1, 8
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(UNIT=temp, fmt='(i4, i2, i2, a2)') uyr, umo, uda, initcode
  read(UNIT=temp, fmt='(10a1)') ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do
  write(UNIT=temp, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', '06', '.sg'
  read (UNIT=temp, fmt='(80a1)') (fsubs(i), i=1,21)
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(UNIT=temp, fmt='(80a1)') (fbase(i), i=1,c), (fdir(i), i=1,8),  &
                       (ftime(i), i=1,10), (fsubs(i), i=1,21)

  read(UNIT=temp, fmt='(a80)') name
  return
!EOC
end subroutine gdasfileF06

subroutine get_gdasnew_diff()
   use lisdrv_module, only : lis, grid, gindex
!   use lisdrv_module, only : lis, tile
!   use tile_spmdMod
   use grid_spmdMod
   use spmdMod, only : iam, masterproc

   use lis_indices_module, only : lis_nc_data, lis_nr_data, lis_tnroffset

   implicit none

   real, allocatable, dimension(:,:) :: elevdiff
   integer :: ierr, i, j

   allocate(elevdiff(lis_nc_data,lis_nr_data), stat=ierr)
   call check_error(ierr,'Error allocating elev diff.',iam)

   elevdiff = 0.0

   call readelevdiff(lis%d%elev, elevdiff)
        
   call lis_log_msg('MSG: get_gdasnew_diff -- done reading elevation '// &
                    'difference file')

!   do i=1,lis%d%nch
!      tile(i)%elev = elevdiff(tile(i)%col, tile(i)%row-lis_tnroffset)
!   enddo
   do j = 1, lis_nr_data
      do i = 1, lis_nc_data
         if ( elevdiff(i,j) == -9999.0 ) then
            elevdiff(i,j) = 0.0
         endif
         if ( gindex(i,j) /= -1 ) then
            grid(gindex(i,j))%elev=elevdiff(i,j)
         endif
      enddo
   enddo

   deallocate(elevdiff,stat=ierr)
   call check_error(ierr,'Error deallocating elev diff.',iam)
 end subroutine get_gdasnew_diff


 subroutine update_gdas_elevdiff(period)

   use lisdrv_module, only : lis

   implicit none

   character(len=1), intent(in) :: period

   integer :: c
   character(len=40) :: elevfile, fpart1, fpart2

   lis%f%gridchange = 1

   elevfile = lis%p%elevfile
   c = index(elevfile,"gdas")
   fpart1 = elevfile(1:c+3)
   fpart2 = elevfile(c+5:40)
   lis%p%elevfile = trim(fpart1) // period // trim(fpart2)

   print*, 'Use newer elevation difference file: ', lis%p%elevfile
   print*, 'Transitioned to GDAS new grid dimensions.'

   call get_gdasnew_diff()
 end subroutine update_gdas_elevdiff
