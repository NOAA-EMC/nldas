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
! !ROUTINE: retecmwf 
!
! !DESCRIPTION:
!  ECMWF model output variables used to force LDAS fall into 2 categories:
!  --> inst3, Instantaneous values, available every 3 hours
!  --> accum, Time integrated values (accumulations), updated every 3 hours
!
!  NOTE-1: be aware that ECMWF outputs large-scale and convective precipitation
!  separately.  For total precipitation, need to sum the two fields, LSP+CP=TP.
!  Reversed traditional order of model-forcing precip fields for easier
!  processing -- interchange at end of routine to preserve traditional order.
!  NOTE 2: only time2 SW flux accumulations used in interpolation
!  NOTE-3: read in of Albedo is currently suppressed.  This field is 
!  instantaneous and available every 6 hours.  At this time, all LDAS LSMs 
!  replace this parameter.
!
!  ECMWF FORCING VARIABLES:
!  1. T         inst3, near-surfcae temperature, ~10 metres [K]
!  2. Q         inst3, near-surface specific humidity, ~10 metres[kg/kg]
!  3. SSRD      accum, surface solar radiation downwards [W m**-2 s]
!  4. STRD      accum, surface thermal radiation downwards [W m**-2 s]
!  5. U         inst3, zonal wind,~10 metres [m/s]
!  6. V         inst3, meridional wind,~10 metres[m/s]
!  7. SP        inst3, surface pressure [Pa]
!  8. CP        accum, convective precipitation [m]
!  9. LSP       accum, large scale precipitation [m]
!
! !REVISION HISTORY:
!  18 Jun 2003: Urszula Jambor; original code
!
! !INTERFACE:
subroutine retecmwf( order, yr, mon, da, hr, ferror )
! !USES:
  use lisdrv_module, only : lis, grid, gindex
  use baseforcing_module, only : glbdata1, glbdata2
  use time_manager
  use ecmwfdomain_module, only : ecmwfdrv
!EOP  
  implicit none
  
  integer, intent(in)    :: order  ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: yr,mon,da,hr      ! data and hour (multiple of 3)
  integer, intent(inout) :: ferror           ! set to zero if there's an error

  character(200) :: dir
  integer :: i,j,r,c
  integer :: iv,nforce,count
  integer :: errorflag,endloop
  integer :: necmwf,nlis
  integer, dimension(ecmwfdrv%nmif) :: pds5
  integer :: gridDesc(200)
  integer :: fret,iret
  integer :: mo
  logical*1 :: lb(ecmwfdrv%ncold*ecmwfdrv%nrold)

  real, dimension(ecmwfdrv%ncold*ecmwfdrv%nrold) :: f,tp
  real, dimension(lis%d%lnc*lis%d%lnr) :: lis_1d
  real, dimension(lis%d%lnc,lis%d%lnr) :: varfield

  !=== End Defining Local Variables =====================================

  !=== set GRIB parameter specifier
  pds5 = (/ 130,133,169,175,131,132,134,143,142 /)
  mo = lis%d%lnc*lis%d%lnr

  ferror = 1  !if there's a problem then ferror is set to zero
  errorflag = 0
  dir = ecmwfdrv%ecmwfdir
  necmwf = ecmwfdrv%ncold*ecmwfdrv%nrold
  nlis = lis%d%lnc*lis%d%lnr
  nforce = ecmwfdrv%nmif
  endloop = 0

  do iv=1,nforce
     if (endloop == 1) exit

     f = 0.0
     varfield = 0.0
     lb = .true.
     fret = 0
     select case (iv)
     case(1,2,5,6,7)
!        print*, iv,', inst3'
        call ret_inst3(dir,yr,mon,da,hr,necmwf,f,pds5(iv),gridDesc,lb,fret)
        if (fret == 0) then
           call interp_iv(iv,lis%d%gridDesc,necmwf,nlis,lb,f,lis_1d,mo,iret)
        endif
     case(3,4,8,9)
!        print*, iv,', accum'
        if (order==2) then
           call ret_accum(dir,yr,mon,da,hr,necmwf,f,pds5(iv),gridDesc,lb,fret)
        endif
        if ((order==2).and.(fret==0)) then
           if (iv == 8) tp = f
           if (iv == 9) then ! CP + LSP = total precipitation
              tp = tp+f
              call interp_iv(iv,lis%d%gridDesc,necmwf,nlis,lb,tp,lis_1d,mo,iret)
           else
              call interp_iv(iv,lis%d%gridDesc,necmwf,nlis,lb,f,lis_1d,mo,iret)
           endif
        endif
        if (order==1) lis_1d = lis%d%udef
     end select

     if (fret == 0) then
        !=== Create 2D array for main program
        count = 0
        do r=1,lis%d%lnr
           do c=1,lis%d%lnc
              varfield(c,r) = lis_1d(c+count)
           end do
           count = count+lis%d%lnc
        end do
     else
        !=== Problem with retrieving file
        errorflag = 1
     end if

     if (errorflag == 1) then
        endloop = 1
        ferror = 0
        print*, 'Could not find forcing for parameter', iv,yr,'/',mon,'/',da,hr
     else
        do c=1,lis%d%lnc
          do r=1,lis%d%lnr
	     if (gindex(c,r) .ne. -1) then
              if (order == 1) then
                 glbdata1(iv,gindex(c,r)) = varfield(c,r)
              else
                 glbdata2(iv,gindex(c,r)) = varfield(c,r)
              end if
	     end if
           end do
        end do
     end if !errorflag check

  end do !iv

  if (order==2) then
     do r=1,lis%d%lnr
        do c=1,lis%d%lnc
	   if (gindex(c,r) .ne. -1) then
            varfield(c,r) = glbdata2(9,gindex(c,r))
            glbdata2(9,gindex(c,r)) = glbdata2(8,gindex(c,r))
            glbdata2(8,gindex(c,r)) = varfield(c,r)
	   endif
        enddo
     enddo
  endif

end subroutine retecmwf
  
!BOP
! !ROUTINE: ret_inst3
!
! !DESCRIPTION: 
! This routine opens the corresponding ECMWF data file to extract
! the specified variable, which represents an instantaneous value.
! Should be used for near-surface temperature, specific humidity,
! winds, and surface pressure.
!
! INPUT:
!   dir         -- directory path name for ECMWF data
!   id          -- index indicates forcing variable (T,Q,U,V,SP) 
!   yr,mo,da,hr -- date information
!   necwmf      -- number of elements of input grid
! OUTPUT:
!   f           -- 1D forcing field, necmwf elements
!   gridDesc        -- GDS array describing grid
!   lb          -- unpacked bitmap
!   fret        -- integer error return flag
! !INTERFACE:
subroutine ret_inst3(dir,yr,mo,da,hr,necmwf,f,id,gridDesc,lb,fret)
!EOP
  implicit none

  character(200), intent(in) :: dir
  integer, intent(in)        :: yr,mo,da,hr,necmwf, id
  integer, intent(out)       :: gridDesc(200)
  integer :: fret
  real, intent(out)          :: f(necmwf)
  logical*1, intent(out)     :: lb(necmwf)

  integer :: remainder
  integer :: iret,gbret,jret
  integer :: i,j,lugb,jpds(200),jgds(200)
  integer :: lubi,kf,k,kpds(200)

  character :: cyr*4,cmo*2,cda*2,chr*2,fhr*2
  character(200) :: filename,fullpath

  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mo
  write(cda, '(i2.2)') da
  write(chr, '(i2.2)') hr
  remainder = modulo(hr,2)
  if (remainder==0) then
     !=== Use analysis field, rather than forecast
     fhr = chr
  else if ((hr==03).or.(hr==15)) then
     write(fhr, '(i2.2)') hr-3
  else if ((hr==09).or.(hr==21)) then
     write(fhr, '(i2.2)') hr-9
  end if
  filename = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
  fullpath = trim(dir)//'/'//cyr//cmo//'/'//trim(filename)

  !=== Set up to open file and retrieve specified field 
  lugb = id+da  ! must alternate file unit # 
  fret=0;  iret=0;  gbret=0;  jret=0
  j = 0
  jpds = -1
  jpds(5) = id
  call baopen(lugb,fullpath,iret)
  if (iret==0) then
     call getgb(lugb,lubi,necmwf,j,jpds,jgds,kf,k,kpds,gridDesc,lb,f,gbret)
!     print*,'avg ',sum(f)/size(f), ' max ',maxval(f),' min ', minval(f)
  else
     gbret = 99
  endif
  call baclose(lugb,jret)
  fret = iret+gbret+jret

end subroutine ret_inst3


!==========================================================================
! ret_accum
!
! This routine opens the corresponding ECMWF data file to extract
! the specified variable, which represents an accumulated value.
! Should be used for radiation & precipitation.
! Returns field valid for 3-hour interval prior to given time, except
! for downward sw flux, where interval length varies.
! Returns forcing in units appropriate to LDAS: 
!      [W m**-2 s] ---> [W m**-2]
!      [m] -----------> [mm/s]
!
! INPUT:
!   dir         -- directory path name for ECMWF data
!   id          -- index indicates forcing variable (SSRD,STRD,LSP,CP) 
!   yr,mo,da,hr -- date information
!   necwmf      -- number of elements of input grid
! OUTPUT:
!   f           -- 1D forcing field, necmwf elements
!   gridDesc        -- GDS array describing grid
!   lb          -- unpacked bitmap
!   fret        -- integer error return flag
!==========================================================================

subroutine ret_accum(dir,yr,mo,da,hr,necmwf,f,id,gridDesc,lb,fret)

  use time_manager, only : tick

  implicit none

  !=== Begin variable declarations
  character(200) :: dir
  integer :: yr,mo,da,hr
  integer :: necmwf,id,fret
  integer :: gridDesc(200) 
  real :: f(necmwf)
  logical*1 :: lb(necmwf)

  integer :: iyr,imo,ida,ihr,imn,iss,ts,idoy
  real*8 :: itime
  real :: igmt
  real :: f1(necmwf),f2(necmwf)

  integer :: iret,gbret,jret,c
  integer :: j,lugb,jpds(200),jgds(200)
  integer :: lubi,kf,k,kpds(200)

  character :: cyr*4,cmo*2,cda*2,chr*2,fhr*2,chr3*2
  character :: ciyr*4,cimo*2,cida*2
  character(200) :: file1,full1,file2,full2

  !=== End variable declarations

  !=== Establish file name
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mo
  write(cda, '(i2.2)') da
  write(chr, '(i2.2)') hr
  write(chr3,'(i2.2)') hr-3
  if (hr==0) then
     !=== roll back one day for forecast initialization time
     iyr=yr;  imo=mo;  ida=da
     ihr=hr;  imn=0;   iss=0
     ts = -24*60*60
     call tick(itime,idoy,igmt,iyr,imo,ida,ihr,imn,iss,ts)
     write(ciyr, '(i4.4)') iyr
     write(cimo, '(i2.2)') imo
     write(cida, '(i2.2)') ida
     file2 = 'ecmwf.'//ciyr//cimo//cida//'12.'//cmo//cda//'00.1_4'
     file1 = 'ecmwf.'//ciyr//cimo//cida//'12.'//cimo//cida//'21.1_4'
  else if ((hr==03).or.(hr==15)) then
     write(fhr, '(i2.2)') hr-3
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'none'
  else if ((hr==06).or.(hr==18)) then
     write(fhr, '(i2.2)') hr-6
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  else if ((hr==09).or.(hr==21)) then
     write(fhr, '(i2.2)') hr-9
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  else if (hr==12) then
     write(fhr, '(i2.2)') hr-12
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  end if
  if (hr==0) then
     full2 = trim(dir)//'/'//ciyr//cimo//'/'//trim(file2)
     full1 = trim(dir)//'/'//ciyr//cimo//'/'//trim(file1)
  else
     full2 = trim(dir)//'/'//cyr//cmo//'/'//trim(file2)
     full1 = trim(dir)//'/'//cyr//cmo//'/'//trim(file1)
  end if
  print*, trim(full2)
  print*, trim(full1)

  !=== Set up to open file(s) and retrieve specified field 
  !=== file 2
  lugb = id+da  ! must alternate file unit # 
  fret=0;  iret=0;  gbret=0;  jret=0
  j = 0
  jpds = -1
  jpds(5) = id
  call baopen(lugb,full2,iret)
  if (iret==0) then
     call getgb(lugb,lubi,necmwf,j,jpds,jgds,kf,k,kpds,gridDesc,lb,f2,gbret)
  else
     gbret = 99
  endif
  call baclose(lugb,jret)
  fret = iret+gbret+jret
  !=== get file 1, if id is not 169 (SW field)
  !=== or if not hr 03 or hr 15
  
  if (id /= 169) then
     if (.NOT.((hr == 03).or.(hr == 15))) then
        lugb = id+da+5  ! must alternate file unit # 
        fret=0;  iret=0;  gbret=0;  jret=0
        j = 0
        jpds = -1
        jpds(5) = id
        call baopen(lugb,full1,iret)
        if (iret==0) then
           call getgb(lugb,lubi,necmwf,j,jpds,jgds,kf,k,kpds,gridDesc,lb,f1,gbret)
        else
           gbret = 99
        endif
        call baclose(lugb,jret)
        fret = fret+iret+gbret+jret
     else
        f1 = 0.0
     endif ! issues if hr is 03 or 15
  end if !any field but SW
  
  !=== Need to address issue of units and time period of interest
  !=== Downward shortwave flux is converted from an accumulation
  !=== to a mean based on a 3, 6, 9, or 12 hour interval.
  !=== Downward longwave flux, large-scale and convective precipitation
  !=== are extracted from 2 adjacent accumulations and converted to a
  !=== 3-hr mean flux or rate.
  if (id == 169) then
     if ((hr==3).or.(hr==15)) then
        f = f2 / (3.0*60*60)
     else if ((hr==6).or.(hr==18)) then
        f = f2 / (6.0*60*60)
     else if ((hr==9).or.(hr==21)) then
        f = f2 / (9.0*60*60)
     else if ((hr==0).or.(hr==12)) then
        f = f2 / (12.0*60*60)
     endif
  else
     f = (f2 - f1) / (3.0*60*60)     
     if ((id == 142).or.(id == 143)) then
!=== Added because sometimes the f2 data is less than the f1 data
!=== creating negative precip. This ocurs only for small values of precip
!=== and the differences are on the order of 10E-4 so for intents and purposes 
!=== is not a substantial problem. This is an issue at the raw processing data 
!=== level issue that will need to be looked at more closely. Reprocessing 
!=== of entire dataset may be necessary if it is found that the problem can
!== be fixed. Will revisit. JG 3/10/2004.
        do c=1,necmwf
        if ( (f2(c)-f1(c)) < 0.0 ) then
	  f(c) = 0.0
	else
	  f(c) = f(c) * 1000.0
	endif
	enddo
!=== Original code line commented out	
!!!     f = f * 1000.0  !convert m/s to mm/s
     endif
  endif

end subroutine ret_accum


!=========================================================================
!===
!=== interp_iv
!===
!=== Interpolate forcing field to LDAS grid
!===
!=== INPUT:
!===   iv             -- field identifier
!===   gridDesc           -- gds of data to be interpolated
!===   ldas_gridDesc      -- gds of ldas grid to be interpolated to
!===   necmwf, nlis -- number of elements in old and new grids
!===   lb             -- input bitmap
!===   f              -- original data array to be interpolated
!=== OUTPUT:
!===   lis_1d         -- interpolated data array
!===   iret           -- integer return value
!=========================================================================

subroutine interp_iv (iv,lis_gridDesc,necmwf,nlis,lb,f,lis_1d,mo,iret)

  use bilinear_interpMod, only : w110,w120,w210,w220,n110,n120,n210,n220,&
         rlat0,rlon0
  use ecmwfdomain_module, only : mi
 
  implicit none

  !=== Begin variable declarations
  integer :: iv, c
  real :: lis_gridDesc(50)
  integer :: necmwf, nlis, iret
  real :: f(necmwf)
  real, dimension(nlis) :: lis_1d,rlat,rlon
  integer :: ip,ipopt(20),ibi,km
  integer :: no,ibo
  integer :: mo
  logical*1 :: lb(necmwf),lo(nlis)
  !=== End variable declarations
  !=== Set ipolates options:
  !=== ip=0, bilinear === ip=3, budget bilinear, for precip
  !=== km=1, one parameter === ibi=1, use undefined bitmap
  if ((iv == 8).or.(iv == 9)) then
     ip = 3
     ipopt = 0
     ipopt(1) = -1
     ipopt(2) = -1
  else
     ip = 0
     ipopt = 0
  endif
  km = 1
  ibi = 1
  lo = .true.  ! initialize output bitmap
!!!  call ipolates (ip,ipopt,gridDesc,ldas_gridDesc,necmwf,nlis, &
!!!       km,ibi,lb,f,no,rlat,rlon,ibo,lo,lis_1d,iret)

! do c=1,mo
!  print*, c, rlat0(c),rlon0(c)
!    print*, w110(c),w120(c),w210(c),w220(c)
! 55  format(f10.5,1x,f10.5)
! enddo

! stop
 call bilinear_interp(lis_gridDesc,ibi,lb,f,ibo,lo,lis_1d,necmwf,mo,&
        rlat0,rlon0,w110,w120,w210,w220,n110,n120,n210,n220,iret)

end subroutine interp_iv


