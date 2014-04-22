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
! !ROUTINE: vic_singleout.F90
!
! !DESCRIPTION:  
!  LIS VIC data writer:  Binary and stat files in ALMA convention
!
! !REVISION HISTORY:
! 08 Apr 2004: James Geiger; ALMA version 
! 
! !INTERFACE:
subroutine vic_singleout(index, tmp)
! !USES:
!  use netcdf
  use lisdrv_module, only : lis,tile
  use vic_varder      ! VIC-specific variables
  use drv_output_mod, only : drv_writevar_bin
  
  implicit none 
!EOP

  real,intent(in)    :: tmp(lis%d%glbnch)
  integer,intent(in) :: index

  integer :: t,c,r,m,i,n,iret,ftn,length
  character*80 mkfyrmo,filenmt,filenmg,cdir,namet,nameg,filengb
  character*80 mkfyrmo2
  character*1  fname(80),fbase(40),fmkdir(80)
  character*1  ftime(8),fcd(3),frm(3),flats(13),ftimec(4)
  character*1  fyrmodir(26),fsubft(80)
  character*1  fsubfg(80),ftimeb(10),fsubgb(9)
  
  character (len=100) :: fbinname
  character (len=100) :: temp1

  integer,parameter :: nvarsg=38
  character*80 :: vname(nvarsg)
  CHARACTER*80 :: vname1(nvarsg)
  
  character*40 file
  character*80 name

  real, dimension(:,:), allocatable :: g2tmp


  data vname1 / "swnet",        &
                "lwnet",        &
                "qle",          &
                "qh",           &
                "qg",           &
                "rainf",        &
                "snowf",        &
                "evap",         &
                "qs",           &
                "qsb",          &
                "qfz",          &
                "snowt",        &
                "avgsurft",     &
                "radt",         &
                "albedo",       &
                "soiltemp1",    &
                "soiltemp2",    &
                "soiltemp3",    &
                "soilmoist1",   &
                "soilmoist2",   &
                "soilmoist3",   &
                "tveg",         &
                "esoil",        &
                "soilwet",      &
                "rootmoist",    &
                "swe",          &
                "qsm",          &
                "delsoilmoist", &
                "delswe",       &
                "acond",        &
                "wind",         &
                "rainfforc",    &
                "snowfforc",    &
                "tair",         &
                "qair",         &
                "psurf",        &
                "swdown",       &
                "lwdown" /

  data vname / "SWnet(W/m2)",         &
               "LWnet(W/m2)",         &
               "Qle(W/m2)",           &
               "Qh(W/m2)",            &
               "Qg(W/m2)",            &
               "Rainf(kg/m2s)",       &
               "Snowf(kg/m2s)",       &
               "Evap(kg/m2s)",        &
               "Qs(kg/m2s)",          &
               "Qsb(kg/m2s)",         &
               "Qfz(kg/m2s)",         &
               "SnowT(K)",            &
               "AvgSurfT(K)",         &
               "RadT(K)",             &
               "Albedo(-)",           &
               "SoilTemp1(K)",        &
               "SoilTemp2(K)",        &
               "SoilTemp3(K)",        &
               "SoilTemp1(kg/m2)",    &
               "SoilMoist2(kg/m2)",   &
               "SoilMoist3(kg/m2)",   &
               "TVeg(kg/m2s)",        &
               "ESoil(kg/m2s)",       &
               "SoilWet(-)",          &
               "RootMoist(kg/m2)",    &
               "SWE(kg/m2)",          &
               "Qsm(kg/m2s)",         &
               "DelSoilMoist(kg/m2)", &
               "DelSWE(kg/m2)",       &
               "ACond(m/s)",          &
               "Wind(m/s)",           &
               "Rainfforc(kg/m2s)",   &
               "Snowfforc(kg/m2s)",   &
               "Tair(K)",             &
               "Qair(kg/kg)",         &
               "Psurf(Pa)",           &
               "SWdown(W/m2)",        &
               "LWdown(W/m2)" /
!BOC
!-------------------------------------------------------------------------
! Test to see if output writing interval has been reached
!-------------------------------------------------------------------------
  if ( mod(lis%t%gmt,vicdrv%writeintvic) == 0 ) then

     write(unit=temp1,fmt='(i4,i2,i2)')lis%t%yr,lis%t%mo,lis%t%da
     read(unit=temp1,fmt='(8a1)') ftime
     do i=1,8
        if ( ftime(i) == (' ') ) then
           ftime(i)='0'
        endif
     enddo
     write(unit=temp1,fmt='(i4)')lis%t%yr
     read(unit=temp1,fmt='(8a1)')ftimec
     do i=1,4
        if ( ftimec(i) == (' ') ) then
           ftimec(i)='0'
        endif
     enddo
     write(unit=temp1,fmt='(a40)') lis%o%odir
     read(unit=temp1,fmt='(40a1)') (fbase(i),i=1,40)
     c=0
     do i=1,40
        if ( fbase(i) == (' ') .and. c == 0 ) then
           c=i-1
        endif
     enddo
     
     write(unit=temp1,fmt='(a4,i3,a5,i4,a1,i4,i2,i2)')'/EXP', & 
          lis%o%expcode,'/VIC/', & 
          lis%t%yr,'/',lis%t%yr,lis%t%mo,lis%t%da
     read(unit=temp1,fmt='(80a1)') (fyrmodir(i),i=1,25)
     do i=1,25
        if ( fyrmodir(i) == (' ') ) then
           fyrmodir(i)='0'
        endif
     enddo
     
     write(unit=temp1,fmt='(a9)')'mkdir -p '
     read(unit=temp1,fmt='(80a1)')(fmkdir(i),i=1,9)
     
     write(unit=temp1,fmt='(80a1)')(fmkdir(i),i=1,9),(fbase(i),i=1,c), & 
          (fyrmodir(i),i=1,26)
     read(unit=temp1,fmt='(a80)')mkfyrmo
     call system(mkfyrmo)
!----------------------------------------------------------------------
! Generate file name for BINARY output
!----------------------------------------------------------------------
     write(unit=fbinname, fmt='(i4,i2,i2,i2)') lis%t%yr,lis%t%mo, & 
          lis%t%da,lis%t%hr
     read(unit=fbinname,fmt='(10a1)') ftimeb
     do i=1,10
        if ( ftimeb(i) == (' ') ) then
           ftimeb(i)='0'
        endif
     enddo
     if ( lis%o%wout == 1 ) then
        if ( lis%o%wtil == 1 ) then 
           write(unit=fbinname,fmt='(a8)') '.ls4r   '
        else
#if ( defined FARMER_DOG_BONES )
           write(unit=fbinname,fmt='(a8)') '.gd4r   '
#else
           write(unit=fbinname,fmt='(a8)') '.gs4r   '
#endif
        endif
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,8)
     elseif ( lis%o%wout == 2 ) then 
        write(unit=fbinname,fmt='(a8)') '.VIC.grb'
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,8)
     elseif ( lis%o%wout == 3 ) then         
        write(unit=fbinname,fmt='(a8)') '.VIC.nc '
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,8)
     endif

     length = len(trim(vname1(index)))
     write(unit=temp1, fmt='(a10)') vname1(index)
     read(unit=temp1,fmt='(10a1)') (fname(i), i=1,length)

     write(unit=fbinname,fmt='(80a1)')(fbase(i),i=1,c), & 
          (fyrmodir(i),i=1,25), '/',& 
          (ftimeb(i),i=1,10), & 
          (fname(i),i=1,length),&
          (fsubgb(i),i=1,8)
     read(unit=fbinname,fmt='(a80)')filengb

     ftn = 73
     if ( lis%o%wout == 1 ) then
#if ( defined FARMER_DOG_BONES )
        allocate(g2tmp(lis%d%lnc,lis%d%lnr))
        g2tmp = lis%d%UDEF
        do i = 1, lis%d%glbnch
           g2tmp(tile(i)%col, tile(i)%row) = tmp(i)
        enddo

        open(ftn,file=filengb,form='unformatted',access='direct', &
             recl=lis%d%lnc * lis%d%lnr * 4)
        write(ftn, rec=1) g2tmp
        close(ftn)
        deallocate(g2tmp)
#else
        open(ftn,file=trim(filengb),form='unformatted')
        call drv_writevar_bin(ftn,tmp)
        close(ftn)
#endif
     endif

!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
#if 0
     if ( vicdrv%vicopen == 0 ) then
        file='VICstats.dat'
        call openfile(name,lis%o%odir,lis%o%expcode,file)
        if ( lis%o%startcode == 1 ) then
           open(65,file=name,form='formatted',status='unknown', & 
                position='append')
        else
           open(65,file=name,form='formatted',status='replace')       
        endif
        vicdrv%vicopen=1
     endif
          
     write(65,996)'       Statistical Summary of Vic output for:  ', & 
          lis%t%mo,'/',lis%t%da,'/',lis%t%yr,lis%t%hr,':',lis%t%mn,':',lis%t%ss
996    format(a47,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
       write(65,*)
       write(65,997)
997    format(t27,'Mean',t41,'Stdev',t56,'Min',t70,'Max')
       call vic_writestats(65)
       write(65,*)
       write(65,*)
#endif
    endif
!EOC
  end subroutine vic_singleout
  
