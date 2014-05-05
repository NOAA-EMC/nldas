C-----------------------------------------------------------------------
C postproc: 	LDAS post-processor for VIC model
C		Converts VIC binary output data into Grib format
C Author:	Justin Sheffield
C Org:		Princeton University, Dept. Civil and Env. Eng.
C Email:	justin@princeton.edu
C Date:		Jan-2000
C-----------------------------------------------------------------------

      program postproc
      
      implicit none

      include 'postproc.h'

C general variables
      integer i, j, rec, rec_yester, hour 
      integer funit, status
      parameter(funit=50)

C space and time variables
      integer nx, ny
      parameter (nx=464, ny=224)
      character*8 startdate, yesterday
      character*2 chour

C data arrays: energy balance components
      real    nswrs(nx,ny),nlwrs(nx,ny),lhtfl(nx,ny),shtfl(nx,ny) 
      real    gflux(nx,ny),snohf(nx,ny),dswrf(nx,ny),dlwrf(nx,ny)

c data arrays: water balance components
      real    asnow(nx,ny),arain(nx,ny),evp(nx,ny),ssrun(nx,ny)
      real    bgrun(nx,ny),snom(nx,ny)

c data arrays: surface state variables
      real    snowt(nx,ny),avsft(nx,ny),radt(nx,ny),albdo(nx,ny)
      real    weasd(nx,ny),cwat(nx,ny)

C data arrays: subsurface state variables
      real    soilt1(nx,ny),soilt2(nx,ny),soilt3(nx,ny)
      real    soilm(nx,ny),soilmr(nx,ny),soilm1m(nx,ny)
      real    soilm1(nx,ny),soilm2(nx,ny),soilm3(nx,ny)
      real    lsoil1(nx,ny),lsoil2(nx,ny),lsoil3(nx,ny)
      real    mstav(nx,ny),mstavr(nx,ny)

C data arrays: evaporation components
      real    evcw(nx,ny),trans(nx,ny),evbs(nx,ny),sbsno(nx,ny)
      real    acond(nx,ny),lai(nx,ny)
     
C data arrays: cold season processes
      real    snod(nx,ny),snoc(nx,ny),salbd(nx,ny)

C file names
      character*256 ctrl_file, bin_file

C general VIC variables 
      character*256 vic_output
      integer nrecs, dt, statehour
      integer startyear, startmonth, startday, starthour
      integer endyear, endmonth, endday, endhour
      
c postproc options
      character*256 mask_file, pmask_file, pp_output,
     &              lats_table, vars_table
      integer pmask_value
      logical*1 mask(nx,ny)

C local variables mask(nx,ny)
      integer pmask(nx,ny)
   
C date stuff
      integer dmy_hour(MAXRECS), dmy_day(MAXRECS), dmy_month(MAXRECS)
      integer dmy_year(MAXRECS), dmy_day_in_year(MAXRECS)

C output variable stuff
      character*40 var_name(MAXVARS)
      real var_mult(MAXVARS), var_offset(MAXVARS)
      logical*1 var_rate(MAXVARS)
      integer nvars, ivar

C-----------------------------------------------------------------------

C off we go
      if(VERBOSE) then
        write(*,'(a)') 'LDAS Post-Processor for the VIC Model'
        write(*,*)
      endif
        
C get the command line arguments
      call getarg(1, ctrl_file)
      if(ctrl_file .eq. '') call usage()

C read the post-processor options
      call get_options(ctrl_file, nrecs, dt, vic_output,
     &           startyear, startmonth, startday, starthour,
     &           endyear, endmonth, endday, endhour, 
     &           statehour, lats_table, vars_table,
     &           mask_file, pmask_file, pmask_value, pp_output)
      
C make the dates
      call make_dates(nrecs, 
     &   startyear, startmonth, startday, starthour,
     &   endyear, endmonth, endday, endhour, dt,
     &   dmy_hour, dmy_day, dmy_month, dmy_year, dmy_day_in_year)

c set the date strings
      write(startdate, '(i4,i2.2,i2.2)') 
     &              dmy_year(1), dmy_month(1), dmy_day(1)

C read in the ldas mask
      call get_mask(nx, ny, mask_file, pmask_file, 
     &              mask, pmask, pmask_value)

C get the variables
      call get_vars(vars_table, var_name, var_mult, 
     &              var_offset, var_rate, nvars)

C open the VIC binary data files
      do ivar=1,nvars
        write(chour, '(i2.2)') starthour
        bin_file = trim(vic_output)//trim(var_name(ivar))//
     &             '_'//startdate//chour//char(0)
c        open(unit=ivar+funit, file=bin_file, form="formatted", 
c        open(unit=ivar+funit, file=bin_file, form="binary",
        open(unit=ivar+funit, file=bin_file, form="unformatted",
     &       access="direct", recl=1, status="old", iostat=status)
        if( status .gt. 0 ) then
          write(*,*) 'ERROR: Can''t open VIC file: ',trim(bin_file)
          stop
        endif
        if(VERBOSE) then
          write(*,*) 'Opened VIC output file: ',trim(bin_file)
        endif
      enddo

c loop through the records
      do rec=1,nrecs
     
        hour = rec
        if(VERBOSE) then
          write(*,'(a,i3,a,i4.4,i2.2,i2.2,i2.2)') 'rec ',rec,': ',
     &    dmy_year(rec),dmy_month(rec),dmy_day(rec),dmy_hour(rec) 
        endif

C read in the VIC binary output data
        call read_bin(1+funit,nx,ny,mask,pmask,pmask_value,rec,nswrs,
     &                var_mult(1),var_offset(1),var_rate(1),dt,nrecs)
        call read_bin(2+funit,nx,ny,mask,pmask,pmask_value,rec,nlwrs,
     &                var_mult(2),var_offset(2),var_rate(2),dt,nrecs)
        call read_bin(3+funit,nx,ny,mask,pmask,pmask_value,rec,lhtfl,
     &                var_mult(3),var_offset(3),var_rate(3),dt,nrecs)
        call read_bin(4+funit,nx,ny,mask,pmask,pmask_value,rec,shtfl,
     &                var_mult(4),var_offset(4),var_rate(4),dt,nrecs)
        call read_bin(5+funit,nx,ny,mask,pmask,pmask_value,rec,gflux,
     &                var_mult(5),var_offset(5),var_rate(5),dt,nrecs)
        call read_bin(6+funit,nx,ny,mask,pmask,pmask_value,rec,snohf,
     &                var_mult(6),var_offset(6),var_rate(6),dt,nrecs)
        call read_bin(7+funit,nx,ny,mask,pmask,pmask_value,rec,dswrf,
     &                var_mult(7),var_offset(7),var_rate(7),dt,nrecs)
        call read_bin(8+funit,nx,ny,mask,pmask,pmask_value,rec,dlwrf,
     &                var_mult(8),var_offset(8),var_rate(8),dt,nrecs)

        call read_bin(9+funit,nx,ny,mask,pmask,pmask_value,rec,asnow,
     &                var_mult(9),var_offset(9),var_rate(9),dt,nrecs)
        call read_bin(10+funit,nx,ny,mask,pmask,pmask_value,rec,arain,
     &                var_mult(10),var_offset(10),var_rate(10),dt,nrecs)
        call read_bin(11+funit,nx,ny,mask,pmask,pmask_value,rec,evp,
     &                var_mult(11),var_offset(11),var_rate(11),dt,nrecs)
        call read_bin(12+funit,nx,ny,mask,pmask,pmask_value,rec,ssrun,
     &                var_mult(12),var_offset(12),var_rate(12),dt,nrecs)
        call read_bin(13+funit,nx,ny,mask,pmask,pmask_value,rec,bgrun,
     &                var_mult(13),var_offset(13),var_rate(13),dt,nrecs)
        call read_bin(14+funit,nx,ny,mask,pmask,pmask_value,rec,snom,
     &                var_mult(14),var_offset(14),var_rate(14),dt,nrecs)

        call read_bin(15+funit,nx,ny,mask,pmask,pmask_value,rec,snowt,
     &                var_mult(15),var_offset(15),var_rate(15),dt,nrecs)
        call read_bin(16+funit,nx,ny,mask,pmask,pmask_value,rec,avsft,
     &                var_mult(16),var_offset(16),var_rate(16),dt,nrecs)
        call read_bin(17+funit,nx,ny,mask,pmask,pmask_value,rec,radt,
     &                var_mult(17),var_offset(17),var_rate(17),dt,nrecs)
        call read_bin(18+funit,nx,ny,mask,pmask,pmask_value,rec,albdo,
     &                var_mult(18),var_offset(18),var_rate(18),dt,nrecs)
        call read_bin(19+funit,nx,ny,mask,pmask,pmask_value,rec,weasd,
     &                var_mult(19),var_offset(19),var_rate(19),dt,nrecs)
        call read_bin(20+funit,nx,ny,mask,pmask,pmask_value,rec,cwat,
     &                var_mult(20),var_offset(20),var_rate(20),dt,nrecs)

        call read_bin(21+funit,nx,ny,mask,pmask,pmask_value,rec,soilt1,
     &                var_mult(21),var_offset(21),var_rate(21),dt,nrecs)
        call read_bin(22+funit,nx,ny,mask,pmask,pmask_value,rec,soilt2,
     &                var_mult(22),var_offset(22),var_rate(22),dt,nrecs)
        call read_bin(23+funit,nx,ny,mask,pmask,pmask_value,rec,soilt3,
     &                var_mult(23),var_offset(23),var_rate(23),dt,nrecs)
        call read_bin(24+funit,nx,ny,mask,pmask,pmask_value,rec,soilm,
     &                var_mult(24),var_offset(24),var_rate(24),dt,nrecs)
        call read_bin(25+funit,nx,ny,mask,pmask,pmask_value,rec,soilmr,
     &                var_mult(25),var_offset(25),var_rate(25),dt,nrecs)
        call read_bin(26+funit,nx,ny,mask,pmask,pmask_value,rec,soilm1m,
     &                var_mult(26),var_offset(26),var_rate(26),dt,nrecs)
        call read_bin(27+funit,nx,ny,mask,pmask,pmask_value,rec,soilm1,
     &                var_mult(27),var_offset(27),var_rate(27),dt,nrecs)
        call read_bin(28+funit,nx,ny,mask,pmask,pmask_value,rec,soilm2,
     &                var_mult(28),var_offset(28),var_rate(28),dt,nrecs)
        call read_bin(29+funit,nx,ny,mask,pmask,pmask_value,rec,soilm3,
     &                var_mult(29),var_offset(29),var_rate(29),dt,nrecs)
        call read_bin(30+funit,nx,ny,mask,pmask,pmask_value,rec,lsoil1,
     &                var_mult(30),var_offset(30),var_rate(30),dt,nrecs)
        call read_bin(31+funit,nx,ny,mask,pmask,pmask_value,rec,lsoil2,
     &                var_mult(31),var_offset(31),var_rate(31),dt,nrecs)
        call read_bin(32+funit,nx,ny,mask,pmask,pmask_value,rec,lsoil3,
     &                var_mult(32),var_offset(32),var_rate(32),dt,nrecs)
        call read_bin(33+funit,nx,ny,mask,pmask,pmask_value,rec,mstav,
     &                var_mult(33),var_offset(33),var_rate(33),dt,nrecs)
        call read_bin(34+funit,nx,ny,mask,pmask,pmask_value,rec,mstavr,
     &                var_mult(34),var_offset(34),var_rate(34),dt,nrecs)

        call read_bin(35+funit,nx,ny,mask,pmask,pmask_value,rec,evcw,
     &                var_mult(35),var_offset(35),var_rate(35),dt,nrecs)
        call read_bin(36+funit,nx,ny,mask,pmask,pmask_value,rec,trans,
     &                var_mult(36),var_offset(36),var_rate(36),dt,nrecs)
        call read_bin(37+funit,nx,ny,mask,pmask,pmask_value,rec,evbs,
     &                var_mult(37),var_offset(37),var_rate(37),dt,nrecs)
        call read_bin(38+funit,nx,ny,mask,pmask,pmask_value,rec,sbsno,
     &                var_mult(38),var_offset(38),var_rate(38),dt,nrecs)
        call read_bin(39+funit,nx,ny,mask,pmask,pmask_value,rec,acond,
     &                var_mult(39),var_offset(39),var_rate(39),dt,nrecs)
        call read_bin(40+funit,nx,ny,mask,pmask,pmask_value,rec,lai,
     &                var_mult(40),var_offset(40),var_rate(40),dt,nrecs)

        call read_bin(41+funit,nx,ny,mask,pmask,pmask_value,rec,snod,
     &                var_mult(41),var_offset(41),var_rate(41),dt,nrecs)
        call read_bin(42+funit,nx,ny,mask,pmask,pmask_value,rec,snoc,
     &                var_mult(42),var_offset(42),var_rate(42),dt,nrecs)
        call read_bin(43+funit,nx,ny,mask,pmask,pmask_value,rec,salbd,
     &                var_mult(43),var_offset(43),var_rate(43),dt,nrecs)

C write the data to Grib format
        call gribout(nx,ny,pp_output,lats_table,mask,
     &        dmy_day(rec),dmy_month(rec),dmy_year(rec),dmy_hour(rec),  
     &        nswrs,nlwrs,lhtfl,shtfl,gflux,snohf,dswrf,dlwrf,
     &        asnow,arain,evp,ssrun,bgrun,snom,
     &        snowt,avsft,radt,albdo,weasd,cwat,
     &        soilt1,soilt2,soilt3,
     &        soilm,soilmr,soilm1m,
     &        soilm1,soilm2,soilm3,
     &        lsoil1,lsoil2,lsoil3,
     &        mstav,mstavr,
     &        evcw,trans,evbs,sbsno,acond,lai,
     &        snod,snoc,salbd) 

      enddo 

      stop
      end
      



