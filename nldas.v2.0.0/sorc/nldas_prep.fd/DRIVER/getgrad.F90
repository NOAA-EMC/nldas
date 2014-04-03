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
! getgrad.F90: 
!
! DESCRIPTION:
!  Opens, reads, interpolates and overlays radiation forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!
! REVISION HISTORY:
!  28  Oct 1999: Brian Cosgrove; Initial code, see getrad.f
!  27  Apr 2000: Brian Cosgrove' Disabled zenith angle correction cutoff for 
!                cos(zen) less than .2
!  11  May 2000: Brian Cosgrove; Enabled correction cutoffs for cos(zen) less
!                than .1, stop model if computed value is greater than 1367 w/m2
!  08  Jan 2001: Brian Cosgrove; Added check to see if czb or czm is equal
!                to zero before trying to divide by czm or czb.  If it is
!                zero, set radiation value to zero
!  06  Mar 2001: Brian Cosgrove; Changed computation of WT1 and WT2 in cases
!                where the previous hour or the next hour of observed radiation
!                is not available.  Substituted TIME1 for LDAS%PINKTIME1 and
!                LDAS%NESTIME1 and TIME2 for LDAS%PINKTIME2 and LDAS%NESTIME2
!  15  Jun 2001: Urszula Jambor; Reworked algorithm for AGRMET data & GLDAS.
!  15  Oct 2001: Jesse Meng; Replace ldas%agrmet flag by ldas%agrmetsw and 
!                ldas%agrmetlw;  added call retagrlw() to calculate 
!                AGRMET LW;
!  25  Feb 2002: Urszula Jambor; check on both SW & LW file status before 
!                using 
!  04  Jun 2002: Urszula Jambor; allowed fall back to model SW.
!  10  Dec 2002: Urszula Jambor; replaced ldas%astat1,2 with local sstat1,2
!                Reorganized routine to mirror other get-routines,and 
!                corrected bug in file status check for initial time1. 
!=========================================================================
subroutine getgrad ( ldas, grid )

  use ldas_module      ! LDAS non-model-specific 1-D variables
  use grid_module      ! LDAS non-model-specific grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid (ldas%nc, ldas%nr)

!=== Local Variables =====================================================
  integer :: c, r, f, zdoy
  integer :: yr1, mo1, da1, hr1, mn1, ss1, doy1, ts1
  integer :: yr2, mo2, da2, hr2, mn2, ss2, doy2, ts2

  integer :: findtime1     ! 0=don't get new file for 1st/2nd time
  integer :: findtime2     ! 1=get a new file 1st/2nd time
  integer :: movetime      ! if 1=move time 2 data into time 1

  integer :: sstat1=0, sstat2=0	!AGRMET SW STATUS
  integer :: lstat1=0, lstat2=0	!AGRMET LW STATUS
  
  real*8 :: time1, time2
  real   :: wt1, wt2, zw1, zw2, czb, cze, czm, gmt1, gmt2
  real   :: obsw(ldas%nc,ldas%nr)

  character*80 :: nameNH, nameSH

!=== End Variable Definition =============================================

!=== Determine Required Observed Radiation Data Times 
  yr1 = ldas%yr    !Previous Hour
  mo1 = ldas%mo
  da1 = ldas%da
  hr1 = ldas%hr
  mn1 = 0
  ss1 = 0
  ts1 = 0
  call tick ( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

  yr2 = ldas%yr    !Next Hour
  mo2 = ldas%mo
  da2 = ldas%da
  hr2 = ldas%hr
  mn2 = 0
  ss2 = 0
  ts2 = 60*60
  call tick ( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!=== Assumption will be not to find or move any data
  findtime1=0
  findtime2=0
  movetime=0

  if(ldas%time.ge.ldas%agrmtime2) then  !crossed time2 bdry, now time1=time2
     movetime=1
     findtime2=1
  endif

  if(ldas%tscount.eq.0 .or. ldas%tscount.eq.1) then  !beginning of the run
     findtime1=1
     findtime2=1
     movetime=0
  endif

  if (findtime1==1) then  !need to get initial time1
     call agrSWfile ( nameNH, nameSH, ldas, yr1, mo1, da1, hr1 )
     sstat1 = 0
!     call retglbSW ( 1, ldas, grid, nameNH, nameSH, sstat1, 1 )
     if (ldas%AGRMETLW /= 0) then
        lstat1 = 0
        call retagrlw ( 1, ldas, grid, yr1, mo1, da1, hr1, lstat1, 1 )
        if ((sstat1 + lstat1) < 2) then
           sstat1 = 0
           lstat1 = 0
        end if
     end if
     if (sstat1 /= 0) ldas%agrmtime1 = time1
  endif ! need initial time1 data

  if(movetime.eq.1) then !transfer time2 data to time1
     ldas%agrmtime1 = ldas%agrmtime2
     sstat1 = sstat2
     lstat1 = lstat2
     do c=1, ldas%nc
        do r=1, ldas%nr
           grid(c,r)%obswdata1 = grid(c,r)%obswdata2
           grid(c,r)%oblwdata1 = grid(c,r)%oblwdata2
        end do !r
     end do !c
  endif ! need to transfer time2 to time1

  if(findtime2.eq.1) then ! need new time2 data
     call agrSWfile ( nameNH, nameSH, ldas, yr2, mo2, da2, hr2 )
     sstat2 = 0
!     call retglbSW ( 2, ldas, grid, nameNH, nameSH, sstat2, 1 )
     if (ldas%AGRMETLW /= 0) then
        lstat2 = 0
        call retagrlw ( 2, ldas, grid, yr2, mo2, da2, hr2, lstat2, 1 )
        if ((sstat2 + lstat2) < 2) then
           sstat2 = 0
           lstat2 = 0
        end if
     endif
     if (sstat2 /= 0) ldas%agrmtime2 = time2
  endif ! need new time2 data 

!=== Print out Status of data holdings
  if (ldas%time == time1) then
    if (sstat1==0) write(79,*) 'NO AGR SW USED',mo1,da1,yr1,hr1
    if (sstat1/=0) write(79,*) 'USED AGRMET SW',mo1,da1,yr1,hr1
    if (sstat2==0) write(79,*) 'NO AGR SW USED',mo2,da2,yr2,hr2
    if (sstat2/=0) write(79,*) 'USED AGRMET SW',mo2,da2,yr2,hr2

    if (ldas%agrmetlw /= 0) then
       if (lstat1==0) write(79,*) 'NO AGR LW USED',mo1,da1,yr1,hr1
       if (lstat1/=0) write(79,*) 'USED AGRMET LW',mo1,da1,yr1,hr1
       if (lstat2==0) write(79,*) 'NO AGR LW USED',mo2,da2,yr2,hr2
       if (lstat2/=0) write(79,*) 'USED AGRMET LW',mo2,da2,yr2,hr2
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Setup fsource array --> FSOURCE only used in NLDAS so far, see getrad.f
!
!  LDAS%FSOURCE(7)=0
!  LDAS%FSOURCE(8)=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

!== Loop through and replace data as possible with AGRMET data
!== This depends on options specified in ldas.crd as well as actual
!== data holdings.
!== **INSERT LOOPING STRUCTURE HERE IF MULTIPLE OBS DATA SETS AVAILABLE

!== AGRMET SW

  if ((sstat1 == 1) .and. (sstat2 == 1)) then
     !== Compute weights and zenith angle information.  Replace forcing 
     !== with AGRMET data.

     wt1 = (ldas%agrmtime2 - ldas%time) / (ldas%agrmtime2 - ldas%agrmtime1)
     wt2 = 1.0 - wt1
     do c=1, ldas%nc
        do r=1, ldas%nr

          !=== Process LAND-only points
          if (grid(c,r)%fimask > 0) then

           zdoy = ldas%doy           
           call zterp ( 1, grid(c,r)%lat, grid(c,r)%lon, gmt1, gmt2, &
                ldas%gmt, zdoy, zw1, zw2, czb, cze, czm, ldas, grid  )
           obsw(c,r) = ldas%udef
           if ((grid(c,r)%obswdata1>0.0) .and. (grid(c,r)%obswdata2>0.0)) then
              obsw(c,r) = grid(c,r)%obswdata1*zw1+grid(c,r)%obswdata2*zw2
              if ((obsw(c,r) > grid(c,r)%obswdata1  .and. &
                   obsw(c,r) > grid(c,r)%obswdata2) .and. &
                   (czb<0.1 .or. cze<0.1)                      ) then
                 obsw(c,r) = grid(c,r)%obswdata1*wt1+grid(c,r)%obswdata2*wt2
              end if
           end if
           if ((grid(c,r)%obswdata1 > 0.0) .and. &
                (grid(c,r)%obswdata2 <= 0.0)      ) then
              if (czb > 0.0) then
                 obsw(c,r) = grid(c,r)%obswdata1 * (czm/czb)
              else
                 obsw(c,r) = grid(c,r)%obswdata1 * (0.0)
              end if
              if ((obsw(c,r) > grid(c,r)%obswdata1 .and. &
                   obsw(c,r) > 0.0) .and. (czb<0.1 .or. cze<0.1)) then
                 obsw(c,r) = grid(c,r)%obswdata1*wt1 + 0.0*wt2
              end if
           end if
           if ((grid(c,r)%obswdata1<=0.0) .and. (grid(c,r)%obswdata2>0.0)) then
              if (cze > 0.0) then
                 obsw(c,r) = grid(c,r)%obswdata2 * (czm/cze)
              else
                 obsw(c,r) = grid(c,r)%obswdata2 * (0.0)
              end if
              if ((obsw(c,r)>0.0 .and. obsw(c,r)>grid(c,r)%obswdata2) &
                   .and. (czb < 0.1  .or.  cze < 0.1)) then
                 obsw(c,r) = 0.0*wt1 + grid(c,r)%obswdata2*wt2
              end if
           end if

           if (obsw(c,r) >= 0.0) then
             if (obsw(c,r) .gt. 1367) then
               print *,'warning, OBSERVED SW RADIATION HIGH'
               print *,'it is',grid(c,r)%forcing(3),' at',c,r
               print *,'agrm1=',grid(c,r)%obswdata1
               print *,'agrm2=',grid(c,r)%obswdata2
               print *,'wt1,wt2,czb,cze,czm,zw1,zw2'
               print *,wt1,wt2,czb,cze,czm,zw1,zw2
               print *,'FALLING BACK TO MODEL SW'
             else
               grid(c,r)%forcing(3) = obsw(c,r)
             end if
           end if

          else
            grid(c,r)%forcing(3)=ldas%udef
          endif !=== Process LAND-only points

        end do !r
     end do !c
  end if

  if ((sstat1 == 1) .and. (sstat2 /= 1)) then
     !== Compute weights and zenith angle information.  Replace forcing 
     !== with zenith extrapolated AGRMET data

     wt1 = (time2 - ldas%time) / (time2 - ldas%agrmtime1)
     wt2 = 1.0 - wt1
     do c=1, ldas%nc
        do r=1, ldas%nr


          !=== Process LAND-only points
          if (grid(c,r)%fimask > 0) then

           zdoy = ldas%doy
           call zterp ( 1, grid(c,r)%lat, grid(c,r)%lon, gmt1, gmt2, &
                ldas%gmt, zdoy, zw1, zw2, czb, cze, czm, ldas, grid )
           obsw(c,r) = ldas%udef
           if (grid(c,r)%obswdata1 > 0.0) then
              if (czb > 0.0) then
                 obsw(c,r) = grid(c,r)%obswdata1  * (czm/czb)
              else
                 obsw(c,r) = grid(c,r)%obswdata1 * (0.0)
              end if
              if ((obsw(c,r) > 400.0) .and. &
                   (czb < 0.1  .or.  cze < 0.1)  ) then
                 !--- Arbitrary cutoff value of 400 W/m2
                 obsw(c,r) = grid(c,r)%obswdata1*wt1
              end if

              if (obsw(c,r) >= 0.0) then
                if (obsw(c,r) .gt. 1367) then
                  print *,'warning, OBSERVED SW RADIATION HIGH'
                  print *,'it is',grid(c,r)%forcing(3),' at',c,r
                  print *,'agrm1=',grid(c,r)%obswdata1
                  print *,'agrm2=',grid(c,r)%obswdata2
                  print *,'wt1,wt2,czb,cze,czm,zw1,zw2'
                  print *,wt1,wt2,czb,cze,czm,zw1,zw2
                  print *,'FALLING BACK TO MODEL SW'
                else
                  grid(c,r)%forcing(3) = obsw(c,r)
                end if
              end if
           end if

          else
            grid(c,r)%forcing(3)=ldas%udef
          endif !=== Process LAND-only points

        end do !r
     end do !c
  end if

  if ((sstat1 /= 1) .and. (sstat2 == 1)) then
     !== Compute weights and zenith angle information.  Replace forcing 
     !== with zenith extrapolated AGRMET data

     wt1 = (ldas%agrmtime2 - ldas%time) / (ldas%agrmtime2 - time1)
     wt2 = 1.0 - wt1
     do c=1, ldas%nc
        do r=1, ldas%nr

          !=== Process LAND-only points
          if (grid(c,r)%fimask > 0) then

           zdoy = ldas%doy
           call zterp ( 1, grid(c,r)%lat, grid(c,r)%lon, gmt1, gmt2, &
                ldas%gmt, zdoy, zw1, zw2, czb, cze, czm, ldas, grid )
           obsw(c,r) = ldas%udef
           if (grid(c,r)%obswdata2 > 0.0) then
              if (cze > 0.0) then
                 obsw(c,r) = grid(c,r)%obswdata2 * (czm/cze)
              else
                 obsw(c,r) = grid(c,r)%obswdata2 * (0.0)
              end if
              if ((obsw(c,r) > 400.0) .and. &
                   (czb < 0.1  .or.  cze <0.1)   ) then
                 !--- Arbitrary cutoff value of 400 W/m2
                 obsw(c,r) = 0.0*wt1 + grid(c,r)%obswdata2*wt2
              end if

              if (obsw(c,r) >= 0.0) then
                if (obsw(c,r) .gt. 1367) then
                  print *,'warning, OBSERVED SW RADIATION HIGH'
                  print *,'it is',grid(c,r)%forcing(3),' at',c,r
                  print *,'agrm1=',grid(c,r)%obswdata1
                  print *,'agrm2=',grid(c,r)%obswdata2
                  print *,'wt1,wt2,czb,cze,czm,zw1,zw2'
                  print *,wt1,wt2,czb,cze,czm,zw1,zw2
                  print *,'FALLING BACK TO MODEL SW'
                else
                  grid(c,r)%forcing(3) = obsw(c,r)
                end if
              end if
           end if

          else
            grid(c,r)%forcing(3)=ldas%udef
          endif !=== Process LAND-only points

        end do !r
     end do !c
  end if

  if ((sstat1 /= 1) .and. (sstat2 /= 1)) then
     do c=1, ldas%nc
        do r=1, ldas%nr
           obsw(c,r) = ldas%udef
        end do
     end do
  end if

!== AGRMET LW

    if( LDAS%AGRMETLW == 1 ) then

!      print*, 'ldas%time = ', ldas%time
!      print*, 'agrmtime1 = ', ldas%agrmtime1
!      print*, 'agrmtime2 = ', ldas%agrmtime2

     if ((lstat1 == 1) .and. (lstat2 == 1)) then

      wt1 = (ldas%agrmtime2 - ldas%time) / (ldas%agrmtime2  - ldas%agrmtime1)
      wt2 = 1.0 - wt1

      do c=1, ldas%nc
         do r=1, ldas%nr

            if ( (grid(c,r)%OBLWDATA1 > 0.0) .AND. &
	         (grid(c,r)%OBLWDATA2 > 0.0)         ) then
	      grid(c,r)%forcing(4)=  grid(c,r)%oblwdata1*wt1 &
		                   + grid(c,r)%oblwdata2*wt2 
            end if
	   		
         end do
      end do

     end if !((lstat1 == 1) .and. (lstat2 == 1))
    end if !(LDAS%AGRMETLW = 1)

  return

end subroutine getgrad




!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - RAD - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! agrSWfile: 
!
! DESCRIPTION:
!  This subroutine puts together the radiation data filenames for both
!  global hemispheres.
!  filename should be like:
!  /GLDAS1/DATA/AGRMET/SWDN/200103/NH/swdn_2001030100n
!
! REVISION HISTORY:
!  28  Oct 1999: Brian Cosgrove; Initial code
!  18  Jun 2001: Urszula Jambor; Modified for AGRMET data use in GLDAS
!  24  Oct 2001: Jesse Meng; Modified for AGRMET directory structure
!=========================================================================
subroutine agrSWfile ( nameNH, nameSH, ldas, yr, mo, da, hr )

  use ldas_module         ! LDAS non-model-specific 1-D variables
  implicit none
  type (ldasdec) ldas

!=== Local Variables =====================================================

  character*80 :: nameNH, nameSH
  integer :: yr, mo, da, hr

  integer :: i, c

  character*1 :: fname(80), fbase(80), fsub(80)
  character*1 :: ftime(18), fdir(13), fhemi(1)

!=== End Variable Definition =============================================

!=== formats for filename segments
92 format (80a1)
93 format (a80)
94 format (a8, i4, i2, i2, i2)
95 format (18a1)
96 format (a40)
97 format (a1)
98 format (a6, i4, i2, a1)
99 format (13a1)

!=== Generate Northern Hemisphere filename
  open(unit=80, file='temp', form='formatted', access='direct', recl=80)
  write(80,96,rec=1) ldas%agrmdir  !should be /GLDAS?/DATA/AGRMET
  read(80,92,rec=1) (fbase(i), i=1,80)

  write(80,98,REC=1)'/SWDN/',yr,mo,'/'
  read(80,99,rec=1) fdir
  do i = 1, 13
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do
  write(80,94,rec=1) 'NH/swdn_', yr, mo, da, hr
  read(80,95,rec=1) ftime
  do i = 1, 18
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do
  write(80,97,rec=1) 'n'
  read (80,92,rec=1) fhemi(1)
!sets c as the last character position of fbase
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(80, 92, rec=1) (fbase(i),i=1,c), (fdir(i),i=1,13), &
       (ftime(i),i=1,18), fhemi(1)
  read(80, 93, rec=1) nameNH
!  write(*,*) '--- AGRMET ---'
!  write(*,*) nameNH
  close(80)

!=== Generate Southern Hemisphere filename
  open(unit=81, file='temp', form='formatted', access='direct', recl=80)
  write(81,96,rec=1) ldas%agrmdir  !should be /GLDAS?/DATA/AGRMET
  read(81,92,rec=1) (fbase(i), i=1,80)

  write(81,98,REC=1)'/SWDN/',yr,mo,'/'
  read(81,99,rec=1) fdir
  do i = 1, 13
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do
  write(81,94,rec=1) 'SH/swdn_', yr, mo, da, hr
  read(81,95,rec=1) ftime
  do i = 1, 18
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do
  write(81,97,rec=1) 's'
  read (81,92,rec=1) fhemi(1)
!sets c as the last character position of fbase
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(81, 92, rec=1) (fbase(i),i=1,c), (fdir(i),i=1,13), &
       (ftime(i),i=1,18), fhemi(1)
  read(81, 93, rec=1) nameSH
 ! write(*,*) nameSH
  close(81)

  return
end subroutine agrSWfile














