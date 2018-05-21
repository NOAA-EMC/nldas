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
! getreanlecmwf.F90: 
!
! DESCRIPTION:
!  Opens, reads, and interpolates 6-hrly, 1/2 degree Reanalysis 
!  ECMWF forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data 
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! REVISION HISTORY:
!  11 Apr 2002: Urszula Jambor; original code based on getgeos.f
!  22 Oct 2002: Urszula Jambor; Limited SW forcing processing to 
!               land-only grid points
!=========================================================================
subroutine getreanlecmwf(ldas,grid)

  use ldas_module     ! LDAS non-model-specific 1-D variables
  use grid_module     ! LDAS non-model-specific grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc,ldas%nr)


  !==== Local Variables=======================
  integer :: try, ferror
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer :: c,r,f,zdoy,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1,ts1,bdoy,byr,bmo
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2,ts2,bda,bhr,bmn
  real*8 :: time1,time2,dumbtime1,dumbtime2
  real*8 :: timenow,btime
  real*8 :: filetime1,filetime2     ! times that the files correspond to
  real :: wt1,wt2,zw1,zw2,czb,cze,czm,gmt1,gmt2
  integer :: findtime1     ! 0=don't get new file for 1st time (or 2nd)
  integer :: findtime2     ! 1=get a new file 1st time (or 2nd)
  integer :: movetime      ! 1=move time 2 data into time 1
  integer, parameter :: nforce=9  ! # forcing variables
	 
  !=== Assumption will be not to find or move any data
  findtime1=0
  findtime2=0
  movetime=0
	 
  !=== End Variable Definition =======================	 

  !=== Determine Required Data Times (The previous hour & the future hour)

  yr1=ldas%yr    !Time now
  mo1=ldas%mo
  da1=ldas%da
  hr1=ldas%hr
  mn1=ldas%mn
  ss1=0
  ts1=0        
  
  call tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  yr1=ldas%yr    !Previous Hour
  mo1=ldas%mo
  da1=ldas%da
  hr1=6*((ldas%hr)/6)
  mn1=0
  ss1=0
  ts1=0
  call tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2=ldas%yr    !Next Hour
  mo2=ldas%mo
  da2=ldas%da
  hr2=6*((ldas%hr)/6)
  mn2=0
  ss2=0
  ts2=6*60*60

  call tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)	 

  if(timenow.gt.ldas%fmodeltime2) then
     movetime=1
     findtime2=1
  endif

  if(ldas%tscount.eq.0 .or. ldas%tscount.eq.1) then  !beginning of the run
     findtime1=1
     findtime2=1
     movetime=0
  endif
      
  !=== Establish fmodeltime1
  if (findtime1==1) then  !need to get new time1 from the past
     order=1   !Get data for glbdata1
     ferror = 0
     try = 0
     ts1 = -24*60*60
     do
        if ( ferror /= 0 ) then
           exit
        end if
        try = try+1
        call ret_reanlecmwf(order,ldas,grid,yr1,mo1,da1,hr1,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           ldas%fmodeltime1=time1
        else  !ferror still=0, so roll back one day
           call tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        end if
        if ( try > ndays ) then 
           print *, 'ERROR: reanlECMWF data gap exceeds 10 days on file 1'
           STOP
        end if
     end do
  endif
  
  !  Repeat for time 2
  
  if(movetime.eq.1) then !transfer time2 data to time1
     ldas%fmodeltime1=ldas%fmodeltime2	
     findtime2=1 !include to ensure getting new time2 data
     
     do f=1,nforce
        do c=1,ldas%nc
           do r=1,ldas%nr
              grid(c,r)%glbdata1(f)=grid(c,r)%glbdata2(f)
           enddo
        enddo
     enddo
     
  endif  ! if movetime=1
  
  if(findtime2.eq.1) then ! need new time2 data
     order=2   !Get data for glbdata2
     ferror = 0
     try = 0
     ts2 = -24*60*60
     do
        if ( ferror /= 0 ) exit
        try = try+1
        call ret_reanlecmwf(order,ldas,grid,yr2,mo2,da2,hr2,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           ldas%fmodeltime2=time2
        else  !ferror still=0, so roll back one day
           call tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        end if
        if ( try > ndays ) then 
           print *, 'ERROR: reanlECMWF data gap exceeds 10 days on file 2'
           STOP
        end if
     end do
  endif
  
  btime=ldas%fmodeltime1
  call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=ldas%fmodeltime2
  call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  
  !===  Interpolate Data in Time
  wt1=(ldas%fmodeltime2-ldas%time)/(ldas%fmodeltime2-ldas%fmodeltime1)
  wt2=1.0-wt1
  do f=1,nforce
     
     if(f.eq.3) then     ! Time Averaged Shortwave       

        do c=1,ldas%nc
           do r=1,ldas%nr

            !=== Process LAND-only points
            if (grid(c,r)%fimask > 0) then

              zdoy=ldas%doy
              call zterp(0,grid(c,r)%lat,grid(c,r)%lon, &
                         gmt1,gmt2,ldas%gmt,zdoy, &
                         zw1,zw2,czb,cze,czm,ldas,grid)
              grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)*zw1
              
              if(grid(c,r)%fimask.eq.0) grid(c,r)%forcing(f)=ldas%udef
              
              if ((grid(c,r)%forcing(f).ne.ldas%udef).and.  &
                  (grid(c,r)%forcing(f).lt.0)                ) then
                 if (grid(c,r)%forcing(f) > -0.00001) then !arbitrary!!
                    grid(c,r)%forcing(f) = 0.0             !threshold!!
                 else
                    print*, 'Stopping because forcing not udef but lt 0, '
                    print*, f,c,r,grid(c,r)%forcing(f)
                    stop
                 end if
              endif
              
              if (grid(c,r)%forcing(f).gt.1367) then
                 grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)
              endif

            else
             grid(c,r)%forcing(f)=ldas%udef
            endif !=== Process LAND-only points

           enddo
        enddo

     else if (f.eq.4) then     ! Time Averaged Longwave, Block Interpolation
     
        do c=1,ldas%nc
           do r=1,ldas%nr
              grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)
           enddo
        enddo

     else if (f.eq.6) then     ! Set to Zero, f=5 is magnitude of wind
     
        do c=1,ldas%nc
           do r=1,ldas%nr
              grid(c,r)%forcing(f)=0.0
           enddo
        enddo

     else if(f.eq.8.or.f.eq.9) then    ! precip variable Block Interpolation

        do c=1,ldas%nc
           do r=1,ldas%nr
              grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)
           enddo
        enddo

     else     !Linearly interpolate everything else	

        do c=1,ldas%nc
           do r=1,ldas%nr
              grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)*wt1+ &
                                   grid(c,r)%glbdata2(f)*wt2
           enddo
        enddo

     endif
  
  enddo   !the f loop

end subroutine getreanlecmwf
	   






