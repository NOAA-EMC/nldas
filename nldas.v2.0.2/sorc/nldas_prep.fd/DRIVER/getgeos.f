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
! getgeos.f: 
!
! DESCRIPTION:
!  Opens, reads, and interpolates GEOS forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data 
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  25 Oct 1999: Jared Entin; Significant F90 Revision
!  11 Apr 2000: Brian Cosgrove; Fixed name construction error 
!               in Subroutine ETA6HRFILE 
!  27 Apr 2000: Brian Cosgrove; Added correction for use of old shortwave
!               data with opposite sign convention from recent shortwave data.
!               Added capability to use time averaged shortwave & longwave data
!               Altered times which are passed into ZTERP--used to be GMT1 
!               and GMT2, now they are LDAS%ETATIME1 and LDAS%ETATIME2
!  30 Nov 2000: Jon Radakovich; Initial code based on geteta.f
!  17 Apr 2001: Jon Gottschalck; A few changes to allow model init.  
!  13 Aug 2001: Urszula Jambor; Introduced missing data replacement.     
!   5 Nov 2001: Urszula Jambor; Reset tiny negative SW values to zero. 
!  17 Nov 2002: Urszula Jambor; Modified to accomodate both GEOS3 and new 
!               GEOS4 grids -- decision based on calendar date (Nov 1, 2002)
!  24 Feb 2003: Urszula Jambor; Elevation-difference file name switch from 
!               GEOS3 to GEOS4 is no longer hard-wired to 1/4deg.
!  04 Mar 2003: Urszula Jambor; began use of GRIDCHANGE flag to switch grids.
!               
!=========================================================================
      subroutine getgeos(ldas,grid)

      use ldas_module     ! LDAS non-model-specific 1-D variables
      use grid_module     ! LDAS non-model-specific grid variables
      implicit none
      type (ldasdec) ldas
      type (griddec) grid(ldas%nc,ldas%nr)
	
	
!==== Local Variables=======================
      integer :: try, ferror
      integer, parameter :: ndays = 10  ! # days to look back for forcing data
      integer :: c,r,f,zdoy,order,updoy
      integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1,ts1,bdoy,byr,bmo
      integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2,ts2,bda,bhr,bmn
      real*8 :: time1,time2,dumbtime1,dumbtime2
      real*8 :: griduptime,btime
      real*8 :: filetime1,filetime2     ! times that the files correspond to
      character*80 :: name
      character*80 :: prevname ! previous time period's file name
      character*40 :: elevfile, fpart1, fpart2
      real :: wt1,wt2,zw1,zw2,czb,cze,czm,gmt1,gmt2,upgmt
      integer :: findtime1     ! 0=don't get new file for 1st time (or 2nd)
      integer :: findtime2     ! 1=get a new file 1st time (or 2nd)
      integer :: movetime      ! 1=move time 2 data into time 1
      integer :: gt,nforce     ! GEOS forcing file time, # forcing variables
	 
!=== Assumption will be not to find or move any data
      findtime1=0
      findtime2=0
      movetime=0
	 
!=== End Variable Definition =======================	 

!=== Determine the correct number of forcing variables
      if (ldas%tscount .eq. 0) then
        nforce = ldas%nmif
      else
        nforce = ldas%nf
      endif

!=== Determine Required GEOS Data Times (The previous hour & the future hour)
!=== If necessary, set grid upgrade date

      griduptime = 0.0
      if (ldas%gridchange==1) then
         yr1 = 2002  !grid update time
         mo1 = 11
         da1 = 01
         hr1 = 0; mn1 = 0; ss1 = 0
         call date2time(griduptime,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
      endif

      yr1=ldas%yr    !Previous Hour
      mo1=ldas%mo
      da1=ldas%da
      hr1=3*((ldas%hr)/3)
      mn1=0
      ss1=0
      ts1=0
      call tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

      yr2=ldas%yr    !Next Hour
      mo2=ldas%mo
      da2=ldas%da
      hr2=3*((ldas%hr)/3)
      mn2=0
      ss2=0
      ts2=3*60*60

      call tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)	 

      if(ldas%time.gt.ldas%geostime2) then
       movetime=1
       findtime2=1
      endif

      if(ldas%tscount.eq.0 .or. ldas%tscount.eq.1) then  !beginning of the run
       findtime1=1
       findtime2=1
       movetime=0
      endif
      
!=== Set GEOS data flags
      ldas%shortflag=2            !Time averaged SW
      ldas%longflag=2             !Time averaged LW 

!=== Establish geostime1
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
          if (ldas%gridchange==1) then
             !=== If time1 >= griduptime, 2002NOV01, GEOS3->GEOS4
             if (time1>=griduptime) then
                ldas%ncold = 288
             endif
          endif
          call geosfile(name,ldas%geosdir,yr1,mo1,da1,hr1,ldas%ncold)
          call readgeos(order,name,ldas,grid,ferror)
          if ( ferror == 1 ) then !successfully retrieved forcing data
             ldas%geostime1=time1
          else  !ferror still=0, so roll back one day
             call tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
          end if
          if ( try > ndays ) then 
             write(79,*) 'ERROR: GEOS data gap exceeds 10 days'
             print *, 'ERROR: GEOS data gap exceeds 10 days on file 1'
             STOP
          end if
       end do
      endif

!  Repeat for time 2

      if(movetime.eq.1) then !transfer time2 data to time1
       ldas%geostime1=ldas%geostime2	
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
          if (ldas%gridchange==1) then
             !=== If time2 >= griduptime, 2002NOV01, GEOS3->GEOS4
             if (time2>=griduptime) then
                ldas%ncold = 288
             endif
          endif
          call geosfile(name,ldas%geosdir,yr2,mo2,da2,hr2,ldas%ncold)
          call readgeos(order,name,ldas,grid,ferror)
          if ( ferror == 1 ) then !successfully retrieved forcing data
             ldas%geostime2=time2
          else  !ferror still=0, so roll back one day
             call tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
          end if
          if ( try > ndays ) then 
             write(79,*) 'ERROR: GEOS data gap exceeds 10 days'
             print *, 'ERROR: GEOS data gap exceeds 10 days on file 2'
             STOP
          end if
       end do
      endif

      btime=ldas%geostime1
      call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
      btime=ldas%geostime2
      call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)

!===  Interpolate Data in Time
      wt1=(ldas%geostime2-ldas%time)/(ldas%geostime2-ldas%geostime1)
      wt2=1.0-wt1
      do f=1,nforce

       if(f.eq.3) then     !Shortwave
	if (ldas%shortflag.eq.1) then    !Got Instantaneous SW
         do c=1,ldas%nc
          do r=1,ldas%nr

            !=== Process LAND-only points
            if (grid(c,r)%fimask > 0) then

              zdoy=ldas%doy
	      call zterp(1,grid(c,r)%lat,grid(c,r)%lon,
     &             gmt1,gmt2,ldas%gmt,zdoy,
     &             zw1,zw2,czb,cze,czm,ldas,grid)
              grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)*zw1+
     &             grid(c,r)%glbdata2(f)*zw2

c     In cases of small cos(zenith) angles, use linear weighting
C     to avoid overly large weights
                   
              if ( (grid(c,r)%forcing(f).gt.grid(c,r)%glbdata1(f))
     &             .and. (grid(c,r)%forcing(f).gt.grid(c,r)%glbdata2(f))
     &             .and. (czb.lt.0.1.or.cze.lt.0.1))then
                 grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)*wt1+
     &                grid(c,r)%glbdata2(f)*wt2
              endif
              if (grid(c,r)%forcing(f).gt.1367) then
                 grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)*wt1+
     &                grid(c,r)%glbdata2(f)*wt2
              endif

              if ((grid(c,r)%forcing(f).ne.ldas%udef).and.
     &             (grid(c,r)%forcing(f).lt.0) ) then
                 if (grid(c,r)%forcing(f) > -0.00001) then !arbitrary!!
                    grid(c,r)%forcing(f) = 0.0 !threshold!!
                 else
                    print*, 'Stop because SW forcing not udef but < 0, '
                    print*, f,c,r,grid(c,r)%forcing(f)
                    stop
                 end if
              endif

            else
              grid(c,r)%forcing(f)=ldas%udef
            endif !=== Process LAND-only points

          enddo
         enddo
        endif

        if (ldas%shortflag.eq.2) then !Got Time Averaged SW
         do c=1,ldas%nc
          do r=1,ldas%nr

            !=== Process LAND-only points
            if (grid(c,r)%fimask > 0) then

             zdoy=ldas%doy
             call zterp(0,grid(c,r)%lat,grid(c,r)%lon,
     &                  gmt1,gmt2,ldas%gmt,zdoy,
     &                  zw1,zw2,czb,cze,czm,ldas,grid)
             grid(c,r)%forcing(f)=grid(c,r)%glbdata2(f)*zw1

             if ((grid(c,r)%forcing(f).ne.ldas%udef).and.
     &           (grid(c,r)%forcing(f).lt.0) ) then
               if (grid(c,r)%forcing(f) > -0.00001) then !arbitrary!!
                 grid(c,r)%forcing(f) = 0.0              !threshold!!
               else
                 print*, 'Stopping because forcing not udef but lt 0, '
                 print*, f,c,r,grid(c,r)%forcing(f)
                 stop
               end if
	     endif

	     if (grid(c,r)%forcing(f).gt.1367) then
               grid(c,r)%forcing(f)=grid(c,r)%glbdata2(f)
             endif

            else
             grid(c,r)%forcing(f)=ldas%udef
            endif !=== Process LAND-only points

          enddo
         enddo

        endif

	else if(f.eq.8.or.f.eq.9) then    ! precip variable Block Interpolation
         do c=1,ldas%nc
          do r=1,ldas%nr

            if (grid(c,r)%fimask > 0) then
              grid(c,r)%forcing(f)=grid(c,r)%glbdata2(f)
            else
              grid(c,r)%forcing(f)=ldas%udef
            endif

          enddo
         enddo
	    
        else if (f.eq.4) then     !Longwave
         if (ldas%longflag.eq.1) then    !Got Instantaneous LW
          do c=1,ldas%nc
           do r=1,ldas%nr
            grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)*wt1+
     &		                 grid(c,r)%glbdata2(f)*wt2 
           enddo
          enddo
         endif    

         if (ldas%longflag.eq.2) then    !Got Time Averaged LW
          do c=1,ldas%nc
           do r=1,ldas%nr

             if (grid(c,r)%fimask > 0) then
               grid(c,r)%forcing(f)=grid(c,r)%glbdata2(f)
             else
               grid(c,r)%forcing(f)=ldas%udef
             endif

           enddo
          enddo

	 endif
	      
         else     !Linearly interpolate everything else	
          do c=1,ldas%nc
           do r=1,ldas%nr

             if (grid(c,r)%fimask > 0) then

               grid(c,r)%forcing(f)=grid(c,r)%glbdata1(f)*wt1+
     &               grid(c,r)%glbdata2(f)*wt2

             else
               grid(c,r)%forcing(f)=ldas%udef
             endif

           enddo
          enddo  
	    
        endif

      enddo   !the f loop

84       format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)
      if(findtime2.eq.1) then
       write(83,*) 'nxt-name: ',name
       write(83,84) yr1,mo1,da1
      endif

!=== If switched to GEOS4 forcing, use appropriate elevation correction file
      if ((ldas%gridchange==1).and.(ldas%ncold==288)) then
         elevfile = ldas%elevfile
         fpart1 = elevfile(1:21)
         fpart2 = elevfile(23:40)
         ldas%elevfile = trim(fpart1) // "4" // trim(fpart2)
         print*, 'Use newer elevation difference file: ', ldas%elevfile
         write(79,*) 'Transitioned from GEOS3 to GEOS4 grid dimensions.'
         ldas%gridchange=0
      endif

      return 
      end
	   

C
C
!!!!!SSSSS  SUBROUTINES    SUBROUTINES    SUBROUTINES   SSSSS
C
C
!=======================================================
!
!  DESCRIPTION:
!   This subroutine puts together GEOS file name
!
!=======================================================
      subroutine geosfile(name,geosdir,yr,mo,da,hr,ncold)
	   
      implicit none
	   
!=== Local Variables
      
      character*80 name
      character*40 geosdir
      integer yr,mo,da,hr,ncold
      integer i,c
      integer ii,jj
  
      integer uyr,umo,uda,uhr
      character(len=2) :: initcode
      character*1 fbase(80),fsubs(80)
      character*1 ftime(10),fdir(8)
	    
!==== End Variable Definition

      ii = ncold
      jj = 181
     
!=== Put together filename
 91   format(a4,i3,a11,i3)
 92   format(80a1)
 93   format(a80)
 94   format(i4,i2,i2,a2)
 95   format(10a1)
 96   format(a40)
 87   format(a6)
 97   format(a8)
 98   format(a1,i4,i2,a1)
 99   format(8a1)
C
!==== Make variables for the time used to create the file
!====  We don't want these variables being passed out
      uyr=yr
      umo=mo
      uda=da
      uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours

!=== Determine initcode for the hour of the forecast file

C  If the time is 12 or later the file is time stamped
C  with the next day.  So check for that first

      if(uhr<3)then
       initcode = '00'   
      elseif(uhr<6)then
       initcode = '03'
      elseif(uhr<9)then
       initcode = '06'
      elseif(uhr<12)then
       initcode = '09'
      elseif(uhr<15)then
       initcode = '12'
      elseif(uhr<18)then
       initcode = '15'
      elseif(uhr<21)then
       initcode = '18'
      elseif(uhr<24)then
       initcode = '21'
      endif

      open(90,file='temp',form='formatted',
     &	   access='direct',recl=80)

      write(90,96,rec=1) geosdir  !should be some ../DATA/GEOS/BEST_LK
      read(90,92,rec=1) (fbase(i),i=1,80)

 
      write(90,98,rec=1)'/',uyr,umo,'/'
      read(90,99,rec=1) fdir
      do i=1,8
      if(fdir(i).eq.(' ')) fdir(i)='0'
      enddo
	   
      write(90,94,rec=1) uyr,umo,uda,initcode
      read(90,95,rec=1) ftime
      do i=1,10
      if(ftime(i).eq.(' ')) ftime(i)='0'	     
      enddo	   

      if (ncold==360) then
         write(90,97,rec=1)'.GEOS323'
         read(90,92,rec=1) (fsubs(i),i=1,8)
      else
         write(90,87,rec=1)'.GEOS4'
         read(90,92,rec=1) (fsubs(i),i=1,6)
      endif
   
      c=0
      do i=1,80
      if(fbase(i).eq.(' ').and.c.eq.0) c=i-1 
      enddo
	
      if (ncold==360) then
         write(90,92,rec=1) (fbase(i),i=1,c),(fdir(i),i=1,8),
     &        (ftime(i),i=1,10),(fsubs(i),i=1,8)
      else
         write(90,92,rec=1) (fbase(i),i=1,c),(fdir(i),i=1,8),
     &        (ftime(i),i=1,10),(fsubs(i),i=1,6)
      endif

      read(90,93,rec=1) name
      write(679,*)hr,name
      close(90)
      return
      end     
         	    	   
