      subroutine grib2_wrt_g2func(bmap,hour1,yesterday,
     &           today,fld,iptable,lugb,ierr)
!-----------------------------------------------------------------
! ABSTRACT: This routine is to write out a new grib2 file
!   March 2013:     J.Wang
!   23 September 2013: Revised by Youlong Xia for NLDAS output   
!-----------------------------------------------------------------
!
      implicit none
      integer iptable(13) ! ipds table
      integer year, mon, day, hour1 
      integer year0, mon0, day0, hour0
!
      integer, parameter   :: max_bytes=20000000
      integer, parameter   :: nx=464,ny=224
!
      integer listsec0(2)
      integer listsec1(13)
      integer igds(5)
      integer igdstmpllen
      integer ipdstmpllen
      integer idrstmpllen
      integer idrsnum,ibmap,numcoord,ipdsnum,idefnum    
 
      integer,dimension(100) :: igdstmpl
      integer,dimension(100) :: ipdstmpl
      integer,dimension(100) :: idrstmpl
!
      integer ideflist(1)
      real(4) coordlist(1)
!
      character*1 cgrib(max_bytes)
      character*8 yesterday, today

      logical*1 bmap(nx,ny)
!
      real(4),dimension(nx*ny) :: fld
      integer ifilw,i,j,lengrib,lonstt,lonlst,latstt,latlst,ierr
      integer yy,mm,dd,hh,mn,sc,lugb
      real(4) :: dxval
!
! code start
!-----------------------------------------------------------------
!
!-- set file unit
      ifilw=lugb
!
!-----------------------------------------------------------------
      if(hour1.ge.1.and.hour1.le.23) then
      read(today,'(i4,i2,i2)') year0,mon0,day0
      year=year0
      mon=mon0
      day=day0
      if(iptable(2).eq.8) then
      hour0=hour1-1
      else
      hour0=hour1 
      endif
      endif  
   
      if (hour1.eq.0) then
      if(iptable(2).eq.8) then
      hour0=23
      read(yesterday,'(i4,i2,i2)') year0,mon0,day0     
      read(today,'(i4,i2,i2)') year,mon,day
      else 
      read(today,'(i4,i2,i2)') year,mon,day
      hour0=hour1
      year0=year
      mon0=mon
      day0=day
      endif
      endif

!-- Open GRIB2 file
      cgrib=''
!
!-- section 0: indicator section 
      listsec0(1)=iptable(1) !Discipline: table 0.0 
!---- (0:Meteorological;1: Hydrlogical; 2:Land) ----------------------
      listsec0(2)=2         ! grib edition number (2:grib2)
!
!-- section 1: identification section
      listsec1(1)=7  ! Identification of orginating center (Table 0)  (7:ncep)
      listsec1(2)=4  ! Identification of orginating subcenter (ON388-Table C) (4:emc)
      listsec1(3)=2         ! GRIB master tables version number (Table 1.0)  (2: Nov 2003)
      listsec1(4)=1         ! Version number of GRIB local tables used to augment Master Tables (Table 1.1)
      listsec1(5)=1         ! Significance of reference time (Table 1.2) (0:ana 1:start of fcst 2:vrfy 3:obs)
      listsec1(6)=year0     ! Reference time - Year (4 digits)
      listsec1(7)=mon0      ! Reference time - Month
      listsec1(8)=day0      ! Reference time - Day
      listsec1(9)=hour0     ! Reference time - Hour
      listsec1(10)=0        ! Reference time - Minute
      listsec1(11)=0        ! Reference time - Second
      listsec1(12)=0        ! Production status of data (Table 1.3) (0:opn products 1:opn test products)
      listsec1(13)=1        ! Type of processed data (Table 1.4) (0:ana products 1:fcst products 2:ana & fcst 3: cntl fcst)

       call gribcreate(cgrib,max_bytes,listsec0,listsec1,ierr)
!       print*,'gribcreate status=',ierr
!
!-- section 3: grid definition section
      igds(1)=0             ! Source of grid definition (Table 3.0) (0:specified in the code)
      igds(2)=nx*ny         ! Number of grid points in the defined grid
      igds(3)=0             ! Number of octets for optional list of numbers defining number of points 
      igds(4)=0             ! Interpretation of list of numbers defining number of points 
!-- example: Lat/lon
      igds(5)=0            ! Grid definition template number (Table 3.1) (0:Lat/lon)
      igdstmpl=0
      if( igds(5)==0) then
      igdstmpllen=19
!
!-- set up grid definition template 3.0
        igdstmpl=0
        igdstmpl(1)=6       ! Shape of the Earth (Table 3.2) (6:Shape of the Earth = 6,371,229.0 m)
        igdstmpl(8)=nx      ! Ni . number of points along a paralell 
        igdstmpl(9)=ny      ! Nj . number of points along a meridian 
        igdstmpl(10)=0      ! Basic angle of the initial production domain 
        igdstmpl(11)=0      ! Subdivisions of basic angle used to define extreme longitudes and latitudes, and direction increments 
        latstt=25063000
        lonstt=235062000
        latlst=52938000
        lonlst=292937000
        dxval=125000
        igdstmpl(12)=latstt ! La1 - latitude of first grid point
        igdstmpl(13)=lonstt ! Lo1 - longitude of first grid point 
        igdstmpl(14)=48     ! Resolution and component flags (Table 3.3, bits order reversed)
        igdstmpl(15)=latlst ! La2 - latitude of last grid point
        igdstmpl(16)=lonlst ! Lo2 - longitude of last grid point 
        igdstmpl(17)=dxval  ! i direction increment
        igdstmpl(18)=dxval   ! j direction increment
        igdstmpl(19)=64      ! Scanning mode (Table 3.4) (+i,+j,i consecutive,row scan same direction)
      endif 
!
      idefnum=1             ! Used if igds(3) .ne. 0. The number of entries in array ideflist
      ideflist=0            ! Used if igds(3) .ne. 0. number of grid points contained in each row ( or column ), Dummy array otherwise
      call addgrid(cgrib,max_bytes,igds,igdstmpl,igdstmpllen,ideflist,idefnum,ierr)
!      print*,'addgrid status=',ierr
!
!-- section 4: product definition section
      ipdstmpl=0
! ------------ product definition for simultaneous variable ----------
      if(iptable(2).eq.0) then       ! used for simultaneous variables
      ipdsnum=iptable(2)      ! Product Definition Template Number (Table 4.0) 
!(0: at a point in time; 8 for average or accumulation) 
      ipdstmpllen=iptable(3)  ! pdt template length
      ipdstmpl(1)=iptable(4)  ! catogory
      ipdstmpl(2)=iptable(5)  ! parameter
      ipdstmpl(3)=2         ! Type of generating process (Table 4.3) (0:ana, 1:ic, 2:fcst)
      ipdstmpl(4)=0            ! Background generating process identifier 
      ipdstmpl(5)=141          ! Land Data Assimilation and Forecast System identified (ON388TableA) 
      ipdstmpl(6)=0            ! Hours of observational data cutoff after reference time
      ipdstmpl(7)=0            ! Minutes of observational data cutoff after reference time
      ipdstmpl(8)=1            ! Indicator of unit of time range (Table 4.4) (0:minute, 1:hour 2:day)
      ipdstmpl(9)=0            ! Forecast time in units defined by ipdstmpl(8)
      ipdstmpl(10)=iptable(6)  ! Type of first fixed surface (see Code table 4.5) (100:isobaric level)
      ipdstmpl(11)=iptable(7)  ! Scale factor of first fixed surface
      ipdstmpl(12)=iptable(8)  ! Scaled value of first fixed surface
      ipdstmpl(13)=iptable(9)  ! Type of first second surface (see Code table 4.5) (100:isobaric level)
      ipdstmpl(14)=iptable(10) ! Scale factor of second fixed surface
      ipdstmpl(15)=iptable(11) ! Scaled value of second fixed surface
      endif
! ----- product difinition for average or accumulation variable -----
      if(iptable(2).eq.8) then   
      ipdsnum=iptable(2)    ! Product Definition Template
!     Number (Table 4.8) (0: at a point in time; 8 for average or accumulation)
      ipdstmpllen=iptable(3)  ! pdt template length
      ipdstmpl(1)=iptable(4)  ! catogory
      ipdstmpl(2)=iptable(5)  ! parameter
      ipdstmpl(3)=2         ! Type of generating process (Table 4.3)(0:ana, 1:ic, 2:fcst)
      ipdstmpl(4)=0            ! Background generating process identifier
      ipdstmpl(5)=141          ! Land Data Assimilation and Forecast System identified (ON388TableA)
      ipdstmpl(6)=0            ! Hours of observational data cutoff after reference time
      ipdstmpl(7)=0            ! Minutes of observational data cutoff after reference time
      ipdstmpl(8)=1            ! Indicator of unit of time range (Table4.4) (0:minute, 1:hour 2:day)
      ipdstmpl(9)=0            ! Forecast time in units defined by ipdstmpl(8)
      ipdstmpl(10)=iptable(6)  ! Type of first fixed surface (see Code table 4.5) (100:isobaric level)
      ipdstmpl(11)=iptable(7)  ! Scale factor of first fixed surface
      ipdstmpl(12)=iptable(8)  ! Scaled value of first fixed surface
      ipdstmpl(13)=iptable(9)  ! Type of first second surface (see Code table 4.5)
      ipdstmpl(14)=iptable(10) ! Scale factor of second fixed surface
      ipdstmpl(15)=iptable(11) ! Scaled value of second fixed surface
      ipdstmpl(16)=year        ! End year of overall time interval 
      ipdstmpl(17)=mon         ! End month
      ipdstmpl(18)=day         ! End day
      ipdstmpl(19)=hour1       ! End hour
      ipdstmpl(20)=0           ! End minute
      ipdstmpl(21)=0           ! End second
      ipdstmpl(22)=1           ! number of time ranges
      ipdstmpl(23)=255         ! total number of data values missing
      ipdstmpl(24)=iptable(12) ! average or accumulation
      ipdstmpl(25)=1           ! 1: analysis, 2: forecast from table 4.11 
      ipdstmpl(26)=1
      ipdstmpl(27)=1
      ipdstmpl(28)=255
      ipdstmpl(29)=0
      ipdstmpl(30)=0
      endif

      numcoord=0               ! Number of coordinate values after template 
      coordlist=0.             ! Optional list of coordinate values
! ----------- end of Section 4 -----------------------------------------    
!-- section 5: Data Representation Section
      idrstmpl=0
      idrsnum=3             ! Data representation section template number (Table 5.0) (3:Grid Point Data - Complex Packing and Spatial Differencing)
      idrstmpllen=18            ! Length of Data representation section
      idrstmpl(2)=0             ! Binary scale factor
      idrstmpl(3)=iptable(13)   ! Decimal scale factor
      idrstmpl(7)=0             ! Missing value management used (see Code Table 5.5)
      idrstmpl(8)=0             ! Primary missing value substitute
      idrstmpl(9)=0             ! Secondary missing value substitute 
      idrstmpl(17)=2            ! Order of spatial difference (see Code Table 5.6) 
!
!-- section 6:       
      ibmap=0  ! Bit-map indicator (Table 6.0) (0:A bit map applies, 255:A bit map doesn't apply)
!
      call addfield(cgrib,max_bytes,ipdsnum,ipdstmpl,ipdstmpllen,coordlist,
     &numcoord,idrsnum,idrstmpl,idrstmpllen,fld,nx*ny,ibmap,bmap,ierr)
!      print*,'addfield status=',ierr

!-- finalize  GRIB message after all section
!-- adds the End Section ( "7777" )
      call gribend(cgrib,max_bytes,lengrib,ierr)
!      print*,'gribend status=',ierr
!      print*,'length of the final GRIB2 message in octets =',lengrib
!
      call wryte(ifilw, lengrib, cgrib)

      end subroutine grib2_wrt_g2func 

