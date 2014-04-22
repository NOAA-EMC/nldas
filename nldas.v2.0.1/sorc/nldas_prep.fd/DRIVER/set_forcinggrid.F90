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
! set_forcinggrid.F90
!
! DESCRIPTION:
!  Determines whether model baseline forcing grid dimensions change
!  during run.  Sets flag to watch for transition date during run. 
!  Updates ncold and nrold if start time of run is past transition date.
!  Also adjusts elevation file name if needed.
!  GDAS --> grid changes from T170 to T254 on 29OCT2002, 12Z
!           (256x512 --> 384x768)
!  GEOS --> grid changes from 1x1 to 1x1.25 deg on 01NOV2002, 00Z
!           (181x360 --> 181x288)
!
!===========================================================================
! REVISION HISTORY:
!  4 Mar 2003: Urszula Jambor; Initial code
!  9 May 2003: Urszula Jambor; Added elevation file name changes for GDAS.
!===========================================================================
subroutine set_forcinggrid(ldas)

  use ldas_module      ! LDAS non-model-specific 1-D variables
  IMPLICIT NONE
  type (ldasdec)                   ldas              

  !===  Begin declarations

  integer :: yr,mo,da,hr,mn,ss,ts,doy
  real :: gmt
  real*8 :: gridchangetime
  character*40 :: elevfile, fpart1, fpart2

  !=== End declarations

  !=== Calculate end time of run, since begin time is ldas%time
  call date2time(ldas%etime,ldas%edoy,ldas%egmt, &
       ldas%eyr,ldas%emo,ldas%eda,ldas%ehr,ldas%emn,ldas%ess)

  !=== Define time when model grid was changed
  gridchangetime = 0.0; doy=0; gmt=0.0
  if (ldas%fgdas==1) then
     yr=2002
     mo=10
     da=29
     hr=12
     mn=0; ss=0
     call date2time(gridchangetime,doy,gmt,yr,mo,da,hr,mn,ss)
  else if (ldas%fgeos==1) then
     yr=2002
     mo=11
     da=01
     hr=00; mn=0; ss=0
     call date2time(gridchangetime,doy,gmt,yr,mo,da,hr,mn,ss)
  endif

  !=== Establish over which time period LSM will run: 
  !=== AFTER...................reset ncold & nrold to new values
  !=== DURING .................transition occurs during run, flag=1
  !=== or BEFORE grid change...no action required, flag=0

  if (ldas%time >= gridchangetime) then
     !=== entire run performed in period after grid change
     if (ldas%fgdas==1) then
        ldas%nrold = 384
        ldas%ncold = 768
	elevfile = ldas%elevfile
    	fpart1 = elevfile(1:21)
        fpart2 = elevfile(23:40)
        ldas%elevfile = trim(fpart1) // "2" // trim(fpart2)
        print*, 'Use newer elevation difference file: ', ldas%elevfile
        write(79,*) 'Switched GDAS settings from T170 to T254 grid dimensions.'
     endif
     if (ldas%fgeos==1) then
        ldas%ncold = 288
        elevfile = ldas%elevfile
        fpart1 = elevfile(1:21)
        fpart2 = elevfile(23:40)
        ldas%elevfile = trim(fpart1) // "4" // trim(fpart2)
        print*, 'Use newer elevation difference file: ', ldas%elevfile
        write(79,*) 'Switched settings from GEOS3 to GEOS4 grid dimensions.'
     endif
  else if ( (ldas%time <= gridchangetime) &
       .and. (ldas%etime >= gridchangetime) ) then
     !=== grid transition occurs during run
     ldas%gridchange = 1
  else
     !=== entire run occurs before grid transition is an issue
     ldas%gridchange = 0
  endif

end subroutine set_forcinggrid
