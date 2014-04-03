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
! getgdas.f90: 
!
! DESCRIPTION:
!  Opens, reads, and interpolates NCEP-GDAS forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data 
!
!  The idea is to open either the 00 or 03 forecast file associated with
!  the most recent GDAS assimilation (available every 6 hours).
!  Precipitation rates and radiation fluxes will be taken from the F03 and 
!  F06 files, since averages are provided.
!  - if that fails, the strategy for missing data is to go backwards up to
!    10 days to get forcing at the same time of day.
!  When only archived data from DAAC is available, set F06_flag=1 
!  to read in the 06 forecast files, available every 6 hours.
!
! REVISION HISTORY:
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
!  22 Jan 2003: Urszula Jambor; Modified to accomodate new GDAS grid,
!               decision based on calendar date (12Z Oct 29, 2002).
!  09 May 2003: Urszula Jambor; Added elevation file name changes for GDAS.
!===============================================================================

subroutine getgdas(ldas, grid)

  use ldas_module      ! LDAS non-model-specific 1-D variables
  use grid_module      ! LDAS non-model-specific grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc, ldas%nr)	

	
!==== Local Variables=======================

  integer :: F06_flag=0  ! 0 indicates forcing file available every 3 hr
                         ! 1 indicates forcing file available every 6 hr
                         ! with ONLY F06 forecast files (not F00 or F03).
  integer :: c, r, f
  integer :: ferror
  integer :: try
  integer, parameter :: ndays = 10  ! # of days to look back for forcing data
  integer :: order     ! 1 indicates lesser interpolation boundary time
                       ! 2 indicates greater interpolation boundary time
  integer :: F00_flag=0! 1 indicates both F00 and corresponding F06 data needed
                       ! 0 indicates one file contains all parameter fields
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2
  integer :: idoy, iyr, imo, ida, ihr, imn, iss, its
  integer :: bdoy, byr, bmo, bda, bhr, bmn
  integer :: updoy, zdoy
  integer :: movetime  ! 1=move time2 into time1
  integer :: findtime1 ! 0=don't get a new file for 1st time (or 2nd)
  integer :: findtime2 ! 1=get a new file for 2nd time (or 1st)
  integer :: nforce    ! Number of forcing variables for model init. option

  real*8 :: griduptime,time1, time2
  real*8 :: btime, inittime, dumbtime1, dumbtime2
  real :: upgmt,gmt1, gmt2, igmt
  real :: wt1, wt2
  real :: zw1, zw2, czb, czm, cze  ! solar zenith weights & cos(zenith)

  character(len=80) :: name, nameF06="null", nameF03
  character*40 :: elevfile, fpart1, fpart2

!=== End Variable Definition =======================

!=== Determine the correct number of forcing variables
  if (ldas%tscount .eq. 0) then
    nforce = ldas%nmif
  else
    nforce = ldas%nf
  endif

!=== Assumption will be not to find or move any data
  findtime1 = 0
  findtime2 = 0
  movetime  = 0

!=== Determine required GDAS data times 
!=== (previous assimilation & future assimilation hours)
!=== If necessary, set grid upgrade date

  griduptime = 0.0
  if (ldas%gridchange==1) then
    yr1 = 2002     !grid update time
    mo1 = 10
    da1 = 29
    hr1 = 12
    mn1 = 0; ss1 = 0
    call date2time( griduptime,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
  endif

  yr1 = ldas%yr  !previous assimilation/forecast hour
  mo1 = ldas%mo
  da1 = ldas%da
  if ( F06_flag == 0 ) then
     hr1 = 3*(ldas%hr/3)
  else
     hr1 = 6*(ldas%hr/6)
  end if
  mn1 = 0
  ss1 = 0
  ts1 = 0
  call tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

  yr2 = ldas%yr  !next assimilation/forecast hour
  mo2 = ldas%mo
  da2 = ldas%da
  if ( F06_flag == 0 ) then
     hr2 = 3*(ldas%hr/3)
  else
     hr2 = 6*(ldas%hr/6)
  end if
  mn2 = 0
  ss2 = 0
  if ( F06_flag == 0 ) then
     ts2 = 3*60*60
  else
     ts2 = 6*60*60
  end if
  call tick( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!=== Use these if need to roll back time.
  dumbtime1 = time1
  dumbtime2 = time2

!=== Check to see if current time (ldas%time) has crossed past gdastime2,
!    requiring that both gdastime2 be assigned to gdastime1 and a new
!    gdastime2 be set 3 or 6 hours ahead of the current gdastime2.

  if ( ldas%time > ldas%gdastime2 ) then
     movetime  = 1
     findtime2 = 1
  end if

!=== What if current time happens to fall at the beginning of the run?
  if ( ldas%tscount == 1 .or. ldas%tscount == 0) then  ! beginning of the run.
     findtime1 = 1
     findtime2 = 1
     movetime  = 0
  end if

!=== Establish gdastime1
  if ( findtime1 == 1 ) then  !get new time1 from the past
     print *, 'Getting new time1 data'
     ferror = 0
     order = 1
     try = 0
     ts1 = -24*60*60

     do
        if ( ferror /= 0 ) exit
        try = try+1
        if (ldas%gridchange==1) then
          !If time1 >= griduptime, 2002OCT29H12, change GDAS resolution
          if (time1>=griduptime) then
             ldas%nrold = 384
             ldas%ncold = 768
          endif
        endif
        if ( F06_flag == 0 ) then
           call gdasfile   ( name, ldas%gdasdir, yr1, mo1, da1, hr1 )
        else
           call gdasfileF06( name, ldas%gdasdir, yr1, mo1, da1, hr1 )
        end if
	if (try>1) ferror = try+20
        call retgdas( order, ldas, grid, name, nameF06, 0, ferror )

        if ( ferror == 1 ) then  !successfully retrieved forcing data
           ldas%gdastime1 = time1
        else  !ferror still=0, so roll back one day & start again
           call tick( dumbtime1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
        end if
        if ( try > ndays ) then  !data gap exceeds 10 days so stop
           write(79,*) 'ERROR: GDAS data gap exceeds 10 days'
           print *, 'ERROR: GDAS data gap exceeds 10 days on file 1'
           STOP
        end if

     end do  !if no error (ferror NOT zero)
!     print *, 'End find time 1'
  end if  !need new time1


!=== Establish gdastime2
  if ( movetime == 1 ) then  !transfer time2 data to time1
     ldas%gdastime1 = ldas%gdastime2
     findtime2 = 1  !include to ensure getting new time2 data
     do f = 1, nforce
        do c = 1, ldas%nc
           do r = 1, ldas%nr
              grid(c,r)%glbdata1(f) = grid(c,r)%glbdata2(f)
           end do
        end do
     end do
  end if  !movetime1

  if ( findtime2 == 1 ) then  !get new time2 data
     print *, 'Getting new time2 data'
     ferror = 0
     order = 2
     try = 0
     ts2 = -24*60*60
     !===determine if both F00 & F06 files needed
     if ( F06_flag == 0 ) then
        if ( modulo(hr2,2) > 0 ) then
           F00_flag = 0             !odd hr, use F03 file
        else
           F00_flag = 1             !even hr, use F00 & F06
        end if
     end if

     do
        if ( ferror /= 0 ) exit
        try = try+1
        if (ldas%gridchange==1) then
          !If time2 is >= griduptime, 2002OCT29H12, change GDAS resolution
           if (time2>=griduptime) then
             ldas%nrold = 384
             ldas%ncold = 768
          endif
        endif
        if ( F06_flag == 0 ) then
           if ( F00_flag == 0 ) then
              call gdasfile   ( name,    ldas%gdasdir, yr2, mo2, da2, hr2 )
           else
              call gdasfile   ( name,    ldas%gdasdir, yr2, mo2, da2, hr2 )
              call gdasfileF06( nameF06, ldas%gdasdir, yr2, mo2, da2, hr2 )
           end if
        else
           call gdasfileF06( name, ldas%gdasdir, yr2, mo2, da2, hr2 )
        end if
	if (try>1) ferror = try+20
        call retgdas( order, ldas, grid, name, nameF06, F00_flag, ferror )
        if ( ferror == 1 ) then  !successfully retrieved forcing data
           ldas%gdastime2 = time2
        else  !ferror still=0, so roll back one day & start again
           call tick( dumbtime2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )
        end if
        if ( try > ndays ) then  !data gap exceeds 10 days so stop
           write(79,*) 'ERROR: GDAS data gap exceeds 10 days'
           print *, 'ERROR: GDAS data gap exceeds 10 days on file 2'
           STOP
        end if

     end do  !if no error (ferror NOT zero)
!     print *, 'End find time 2'
  end if  !need new time2

!=== Reset GMT times
  btime = ldas%gdastime1
  call time2date( btime, bdoy, gmt1, byr, bmo, bda, bhr, bmn )
  btime = ldas%gdastime2
  call time2date( btime, bdoy, gmt2, byr, bmo, bda, bhr, bmn )
  if ( (F06_flag==0) .and. (F00_flag==1) ) then
!=== Using a 6 hr mean, adjust time1=(time2-6hr)
     inittime = ldas%gdastime2
     call time2date( inittime, idoy, igmt, iyr, imo, ida, ihr, imn )
     its = -6*60*60
     call tick( inittime, idoy, igmt, iyr, imo, ida, ihr, imn, iss, its )
!     print *, 'Modified gmt (-3hr) ', hr1, igmt, hr2, gmt2
  end if

!=== INTERPOLATE DATA IN TIME
!=== - linearly interpolate instantaneous values.
!=== - zenith angle interpolation for SW
!=== - block interpolation for mean value parameters (precip,LW,albedo)

!=== Calculate weights
  wt1 = (ldas%gdastime2 - ldas%time) / (ldas%gdastime2 - ldas%gdastime1)
  wt2 = 1.0 - wt1

  do f = 1, nforce  !loop through all forcing parameters & interpolate

     !===if SW or LW radiation, if negative sign convention change to positive
     if ( (f == 3) .or. (f == 4) ) then
        do c = 1, ldas%nc
           do r = 1, ldas%nr

              if ( (grid(c,r)%glbdata2(f) /= -9999.9) .and.  &
                   (grid(c,r)%glbdata2(f) < 0)                ) then
                 grid(c,r)%glbdata2(f) = (-1) * grid(c,r)%glbdata2(f)
              end if
              if ( grid(c,r)%fimask == 0 ) then 
                 grid(c,r)%glbdata2(f) = ldas%udef
              end if
              if ( (grid(c,r)%glbdata1(f) /= -9999.9) .and.  &
                   (grid(c,r)%glbdata1(f) < 0)                ) then
                 grid(c,r)%glbdata1(f) = (-1) * grid(c,r)%glbdata1(f)
              end if
              if ( grid(c,r)%fimask == 0 ) then 
                 grid(c,r)%glbdata1(f) = ldas%udef
              end if
           end do  !r
        end do  !c
     end if  !f=3 or 4, SW or LW sign convention

     if ( f == 3 ) then  !shortwave

        do c = 1, ldas%nc
           do r = 1, ldas%nr

             !=== Process LAND-only points
             if (grid(c,r)%fimask > 0) then

               zdoy = ldas%doy
               if ( (F06_flag==0) .and. (F00_flag==1) ) then
                !=== Using a 6 hr mean, adjust time1=(time2-6hr)
                 call zterp( 0, grid(c,r)%lat, grid(c,r)%lon, igmt, gmt2, &
                      ldas%gmt, zdoy, zw1, zw2, czb, cze, czm, ldas, grid )
               else
                 call zterp( 0, grid(c,r)%lat, grid(c,r)%lon, gmt1, gmt2, &
                      ldas%gmt, zdoy, zw1, zw2, czb, cze, czm, ldas, grid )
               end if
               grid(c,r)%forcing(f) = zw1 * grid(c,r)%glbdata2(f)
               if (grid(c,r)%forcing(f) < 0) then
                 print *, '2 warning!!!  SW radiation is negative!!'
                 print *, 'sw=', grid(c,r)%forcing(f), '... negative'
                 print *, 'gdas2=', grid(c,r)%glbdata2(f)
                 print *, 'forcing mask=', grid(c,r)%fimask
                 STOP
               end if

	       if (grid(c,r)%forcing(f).gt.1367) then
                 grid(c,r)%forcing(f)=grid(c,r)%glbdata2(f)
               endif

             else
               grid(c,r)%forcing(f)=ldas%udef
             endif !=== Process LAND-only points

           end do  !r
        end do  !c

     else if ( (f==4) .or. (f==10) ) then  !mean longwave or albedo

        do c = 1, ldas%nc
           do r = 1, ldas%nr

             !=== Process LAND-only points
             if (grid(c,r)%fimask > 0) then

               if ( F00_flag==1 ) then !subtract 0-3hr mean from 0-6hr mean
                 grid(c,r)%forcing(f)=2*grid(c,r)%glbdata2(f) -grid(c,r)%glbdata1(f)
               else
                 grid(c,r)%forcing(f) = grid(c,r)%glbdata2(f)
               end if

             else
               grid(c,r)%forcing(f)=ldas%udef
             endif !=== Process LAND-only points

           end do
        end do

     else if ( (f==8) .or. (f==9 ) ) then  !mean precip rate,

        do c = 1, ldas%nc
           do r = 1, ldas%nr

             !=== Process LAND-only points
             if (grid(c,r)%fimask > 0) then

               if ( F00_flag == 1) then !subtract 0-3hr mean from 0-6hr mean
                 if (2*grid(c,r)%glbdata2(f) >= grid(c,r)%glbdata1(f)) then
                    grid(c,r)%forcing(f)=2*grid(c,r)%glbdata2(f) -grid(c,r)%glbdata1(f)
                 else
                    grid(c,r)%forcing(f) = 0.0
                 end if
              else
                 grid(c,r)%forcing(f) = grid(c,r)%glbdata2(f)
              end if

             else
               grid(c,r)%forcing(f)=ldas%udef
             endif !=== Process LAND-only points

           end do
        end do

     else  !linearly interpolate everything else

        do c = 1, ldas%nc
           do r = 1, ldas%nr

             !=== Process LAND-only points
             if (grid(c,r)%fimask > 0) then
               grid(c,r)%forcing(f) = wt1 * grid(c,r)%glbdata1(f) +  &
                                     wt2 * grid(c,r)%glbdata2(f)
             else
               grid(c,r)%forcing(f)=ldas%udef
             endif !=== Process LAND-only points

           end do
        end do

     end if  !a certain parameter f

  end do  !the f-parameter loop


  !=== If switched from T170 to T254 grid,
  !===  use appropriate elevation correction file
  if ((ldas%gridchange==1).and.(ldas%ncold==768)) then
     elevfile = ldas%elevfile
     fpart1 = elevfile(1:21)
     fpart2 = elevfile(23:40)
     ldas%elevfile = trim(fpart1) // "2" // trim(fpart2)
     print*, 'Use newer elevation difference file: ', ldas%elevfile
     write(79,*) 'Transitioned from GDAS T170 t0 T254 grid dimensions.'
     ldas%gridchange=0
  endif

end subroutine getgdas


!
!
!!!!!SSSSS  SUBROUTINES    SUBROUTINES    SUBROUTINES   SSSSS
!
!
!=======================================================
!
!  DESCRIPTION:
!   This subroutine puts together GDAS file name for 
!   3 hour file intervals
!
!=======================================================

subroutine gdasfile( name, gdasdir, yr, mo, da, hr )

  implicit none

!==== Local Variables=======================

  character(len=80) :: name, gdasdir
  integer :: yr, mo, da, hr
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: remainder
  integer :: doy
  real    :: gmt
  real*8  :: dumbtime
  character(len=2) :: initcode, fcstcode
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(21)

!=== End Variable Definition ===============

!=== formats for filename segments
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, a2)
95 format (10a1)
96 format (a40)
97 format (a16, a2, a3)
98 format (a1, i4, i2, a1)
99 format (8a1)

!=== Make variables for the time used to create the file
!=== We don't want these variables being passed out
  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.
  remainder = modulo(uhr,2)  !if even, then remainder equals zero
                             !if odd, then remainder equals one

!=== if hour is   00 or 03, look for 00ZF00 or 00ZF03 
!=== if hour is   06 or 09, look for 06ZF00 or 06ZF03
!=== if hour is   12 or 15, look for 12ZF00 or 12ZF03
!=== if hour is   18 or 21, look for 18ZF00 or 18ZF03

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
  
  open(unit=90, file='temp', form='formatted', access='direct', recl=80)
  write(90, 96, rec=1) gdasdir  !should be /GLDAS4/DATA/GDAS
  read(90, 92, rec=1) (fbase(i), i=1,80)

  write(90, 98, rec=1) '/', uyr, umo, '/'
  read(90, 99, rec=1) fdir
  do i = 1, 8
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(90, 94, rec=1) uyr, umo, uda, initcode
  read(90, 95, rec=1) ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  write(90, 97, rec=1) '.gdas1.sfluxgrbf', fcstcode, '.sg'
  read (90, 92, rec=1) (fsubs(i), i=1,21)

!sets c as the last character position of fbase
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(90, 92, rec=1) (fbase(i), i=1,c), (fdir(i), i=1,8),  &
                       (ftime(i), i=1,10), (fsubs(i), i=1,21)

  read(90, 93, rec=1) name

  close(90)
  return
end subroutine gdasfile



!=======================================================
!
!  DESCRIPTION:
!   This subroutine puts together GDAS file name for 
!   6 hour forecast files (DAAC archive friendly)
!
!=======================================================

subroutine gdasfileF06( name, gdasdir, yr, mo, da, hr )

  implicit none

!==== Local Variables=======================

  character(len=80) :: name, gdasdir
  integer :: yr, mo, da, hr
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: doy
  real    :: gmt
  real    :: dumbtime
  character(len=2) :: initcode, fcstcode
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(21)

!=== End Variable Definition ===============

!=== formats for filename segments
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, a2)
95 format (10a1)
96 format (a40)
97 format (a16, a2, a3)
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

!=== if hour is 00Z look for 18ZF06 after rolling back one day
!=== if hour is 06Z look for 00ZF06
!=== if hour is 12Z look for 06ZF06
!=== if hour is 18Z look for 12ZF06

  fcstcode = '06'
  if ( uhr == 0 ) then
     !=== ROLL BACK ONE DAY ===!
     call tick( dumbtime, doy, gmt, uyr, umo, uda, uhr, umn, uss, ts1 )
     initcode = '18'
  else if ( uhr == 6 ) then
     initcode = '00'
  else if ( uhr == 12 ) then
     initcode = '06'
  else if ( uhr == 18 ) then
     initcode = '12'
  end if

  open(unit=90, file='temp', form='formatted', access='direct', recl=80)
  write(90, 96, rec=1) gdasdir  !should be /GLDAS4/DATA/GDAS
  read(90, 92, rec=1) (fbase(i), i=1,80)

  write(90, 98, rec=1) '/', uyr, umo, '/'
  read(90, 99, rec=1) fdir
  do i = 1, 8
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(90, 94, rec=1) uyr, umo, uda, initcode
  read(90, 95, rec=1) ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  !Was having problems when using variable fcstcode, so filled in the value
  !  write(90, 97, rec=1) '.gdas1.sfluxgrbf', fcstcode, '.sg'
  write(90, 97, rec=1) '.gdas1.sfluxgrbf', '06', '.sg'
  read (90, 92, rec=1) (fsubs(i), i=1,21)

!sets c as the last character position of fbase
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(90, 92, rec=1) (fbase(i), i=1,c), (fdir(i), i=1,8),  &
                       (ftime(i), i=1,10), (fsubs(i), i=1,21)

  read(90, 93, rec=1) name

  close(90)
  return
end subroutine gdasfileF06
