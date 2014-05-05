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
! !ROUTINE: retberg.F90
!
! !DESCRIPTION:
!  For the given time, reads in 8 parameters from 1/2 degree
!  Reanalysis ECMWF data, transforms into 9 GLDAS forcing 
!  parameters and interpolates to the LDAS domain.
!
! Reanal. ECMWF FORCING VARIABLES available every 6 hours: 
!  1. 2T      2 metre air temperature [K], instantaneous
!  2. 2D      2 metre dew point temperature [K], instantaneous
!  3. SSRD    Downward shortwave flux, surface [W/m^2 s], accumulation
!  4. STRD    Downward longwave radiation, surface [W/m^2 s], accumulation
!  5. WIND    Calculated absolute 10 metre wind speed [m/s], instantaneous
!  6. P       Calculate surface pressure [Pa], instantaneous
!  7. LSP+CP  Total precipitation [m], accumulation
!  8. CP      Convective precipatation [m], accumulation
!
! !REVISION HISTORY:
!  09 Apr 2002: Urszula Jambor; Original code, based on readgeos.f
!  20 Dec 2002: Matt Rodell; Don't ipolate if running 1/2 deg, also fix
!               error in setting v-wind to zero
!  25 Mar 2003: Urszula Jambor; Modified argument list passed to GEOGFILL2.
!  24 Nov 2003; Sujay Kumar; Included in LIS
!  15 Mar 2004; Matt Rodell; Fix test for beginning of restart run
!
! !INTERFACE:
subroutine retberg(order, yr, mon, da, hr, ferror)
! USES:
  use lisdrv_module, only : lis, grid,gindex
  use baseforcing_module, only: glbdata1, glbdata2
  use bergdomain_module, only : bergdrv,mi
  use bilinear_interpMod, only : w110,w120,w210,w220,n110,n120,n210,n220, & 
       rlat0,rlon0
  use spmdMod, only : iam
!EOP
  implicit none
  !Subroutine inputs
  integer :: order            ! lower(1) or upper(2) time interval bdry
  integer :: yr,mon,da,hr      ! data and hour (multiple of 6)
  integer :: ferror           ! set to zero if there's an error
  integer :: iii
  integer :: mo
  ! Specify Reanalysis BERG forcing parameters & file paths
  real, parameter :: no_data_value = -9999.0
  integer, parameter :: N_EF=8, NF=9  ! # ECMWF & GLDAS forcing variables
  real, parameter :: epsln = 0.622        !Constants' values taken from 
  real, parameter :: A = 2.53E8, B=5.42E3 !Roger&Yau, A Short Course
                                          !in Cloud Physics, pp.12-17
  character(25), dimension(N_EF), parameter :: ecmwf_fv = (/  &
       'TEMP-AIR   ',    &
       'TEMP-DEW   ',    & 
       'RAD-SW-DWN ',    &
       'RAD-LW-DWN ',    &
       'WIND       ',    &
       'PRES-SRF   ',    &
       'PREC-TOTL  ',    &
       'PREC-CONV  '     /)

  character(100), dimension(N_EF), parameter :: v_inpath = (/ &
       '/X5RAID/sujay/GRID/     ',  &
       '/X5RAID/sujay/GRID/     ',  &
       '/X5RAID/sujay/GRID/     ',  &
       '/X5RAID/sujay/GRID/     ',  &
       '/X5RAID/sujay/GRID/     ',  &
       '/X5RAID/sujay/GRID/     ',  &
       '/X5RAID/sujay/GRID/     ',  &
       '/X5RAID/sujay/GRID/     '   /)

  integer, dimension(N_EF), parameter :: fnum = (/ &
       31, 33, 34, 35, 36, 37, 38, 39 /) 

  character(200) :: infile

  real :: tempecmwf(bergdrv%ncold,bergdrv%nrold,N_EF)
  real :: temp2ecmwf(bergdrv%ncold,bergdrv%nrold,NF)
  real,allocatable :: templdas(:, :)
  real :: p_kPa(bergdrv%ncold,bergdrv%nrold)
  character :: cyr*4, cmo*2 
  integer :: timestep
  integer :: i, j, t, v, ii, r, kk, eindex, x, y
  integer :: ios,ioerror      ! set to non-zero if there's an error
  integer :: istat

  integer :: gldas,necmwf     ! Size of I/O 1D fields
  real    :: f(bergdrv%ncold*bergdrv%nrold),go(lis%d%lnc*lis%d%lnr)  ! 1D I/O fields
  real    :: tg(lis%d%lnc,lis%d%lnr)           ! Interpolated 2D data field
  integer :: ibi,no,ibo,ipopt(20)
  integer :: iret,km,ip,c
  real :: gridDesco(50)       ! Input,output grid info arrays
  logical*1 :: lb(bergdrv%ncold*bergdrv%nrold)
  logical*1 :: lo(lis%d%lnc*lis%d%lnr)      ! Input and output bitmaps
  logical*1 :: geogmask(lis%d%lnc,lis%d%lnr)! 2D output bitmap
  logical*1 :: tmask(lis%d%lnc,Lis%d%lnr)   ! 2d valid temperature bitmap
  logical :: file_open

  ! if a problem, ferror is set to zero
  ferror = 1
  mo = lis%d%lnc*lis%d%lnr
  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour correspond to
  ! and which file to open.
  !----------------------------------------------------------------

  timestep = 4*(da - 1) + (1 + hr/6)
  print*, 'order and month-timestep', order, timestep
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mon
  allocate(templdas(lis%d%ngrid,bergdrv%nmif), stat=ios)
  if(ios .ne.0) then 
     print*, 'Error allocating templdas.',iam
     stop 344
  endif
  !=== Open ECMWF forcing files

  if ((lis%f%findtime1==1) .and. (order==1)) then  !beginning of run
     do v = 1, 7 !N_EF

        ! File name for data variable(v)/year/yearmo
        infile=trim(v_inpath(v))//trim(ecmwf_fv(v))//'/'//cyr//'/'//cyr//cmo
	inquire(unit=fnum(v), opened=file_open, name=infile)
	if (file_open) then
	  write(*,*) 'RETBERG(1): unit ', fnum(v), '= ', infile
	  stop
	end if
        open(fnum(v), file=infile, form='unformatted', iostat=istat)
        if (istat/=0) then
           print *, 'Problem opening file: ', trim(infile)
           print *, 'Stopping...'
           stop
        else
           print*, 'Opened file: ', trim(infile)
        end if
        ! Fast-forward to desired location within data file
        if (timestep>1) then
           do t = 1, (timestep-1)
              do i = 1, bergdrv%nrold
                 read(fnum(v)) (tempecmwf(j,i,v), j=1, bergdrv%ncold)
              end do 
           end do
        end if
     end do !v
     
  else if ( timestep==1 ) then !close previous month files & open current
     do v = 1, 7!N_EF

        close(fnum(v))
        ! New file name for data variable(v)/year/yearmo
        infile=trim(v_inpath(v))//trim(ecmwf_fv(v))//'/'//cyr//'/'//cyr//cmo
	inquire(unit=fnum(v), opened=file_open, name=infile)
	if (file_open) then
	  write(*,*) 'RETBERG(2): unit ', fnum(v), '= ', infile
	  stop
	end if
        open(fnum(v), file=infile, form='unformatted', iostat=istat)
        if (istat/=0) then
           print *, 'Problem opening file: ', trim(infile)
           print *, 'Stopping...'
           stop
        else
           print*, 'Opened file: ', trim(infile)
        end if

     end do !v

  end if
  do v = 1, 7!N_EF
     !----------------------------------------------------------------
     ! Extract grid array for chosen time
     !----------------------------------------------------------------
     do i = 1, bergdrv%nrold
        read(fnum(v)) (tempecmwf(j,i,v), j=1, bergdrv%ncold)
     end do

     !----------------------------------------------------------------
     ! Change data from ECMWF grid convention to GLDAS one 
     ! Swap latitude band and shift longitudes by 180deg. 
     !----------------------------------------------------------------
     call berggrid_2_gldasgrid(bergdrv%ncold,bergdrv%nrold,tempecmwf(:,:,v))
  end do !v
  
  !-----------------------------------------------------------------
  ! Calculate specific humidity from Dew point Temp. & Sfc Pressure.
  ! Approximate: q~epsilon*e/p, p~p_sfc, e=e_s(Td), e_s=A*exp(-B/T)
  !-----------------------------------------------------------------

  where (tempecmwf(:,:,6)>=0.0) &
       p_kPa(:,:) = tempecmwf(:,:,6) / 1000.0
  where (tempecmwf(:,:,6)< 0.0) &
       p_kPa(:,:) = no_data_value
  where (tempecmwf(:,:,2)>no_data_value) &
       tempecmwf(:,:,2) = epsln*(A * EXP(-B/tempecmwf(:,:,2))) / p_kPa(:,:) 
  !-----------------------------------------------------------------
  ! Filter out any unrealistic forcing values.
  !-----------------------------------------------------------------

  ! Shortwave
  where (tempecmwf(:,:,3) < 0.0001) &
         tempecmwf(:,:,3) = 0.0001
  ! Wind
  where (tempecmwf(:,:,5) < 0.0001) &
         tempecmwf(:,:,5) = 0.0001
  ! Total precipitation
  where (tempecmwf(:,:,7) < 0.0) &
         tempecmwf(:,:,7) = 0.0
  ! Convective precipitation
  where (tempecmwf(:,:,8) < 0.0) &
         tempecmwf(:,:,8) = 0.0

  !----------------------------------------------------------------
  ! Transfer ECMWF forcing fields to GLDAS format
  !----------------------------------------------------------------

  temp2ecmwf(:,:,1) = tempecmwf(:,:,1)
  temp2ecmwf(:,:,2) = tempecmwf(:,:,2)
  temp2ecmwf(:,:,3) = tempecmwf(:,:,3)
  temp2ecmwf(:,:,4) = tempecmwf(:,:,4)
  temp2ecmwf(:,:,5) = tempecmwf(:,:,5)    !Since absolute wind speed 
  temp2ecmwf(:,:,6) = 0.0                 ! already provided as WIND,
  temp2ecmwf(:,:,7) = tempecmwf(:,:,6)    ! let U=WIND and V=0.0
  temp2ecmwf(:,:,8) = tempecmwf(:,:,7)/(6*3600)
  temp2ecmwf(:,:,9) = tempecmwf(:,:,8)/(6*3600)

!  open(112,file="forcing.bin",form='unformatted')
!  write(112) temp2ecmwf(:,:,3)
!  print*, temp2ecmwf(:,:,3)
  !-----------------------------------------------------------------
  ! Interpolate each forcing variable to GLDAS domain
  !-----------------------------------------------------------------

  !=== Initialize input & output grid arrays
  gridDesco = 0
           
  !=== Set input & output grid array values (reanlECMWF to GLDAS)
        
  gridDesco = lis%d%gridDesc

  !=== Define input & output data bitmaps
  necmwf = bergdrv%ncold*bergdrv%nrold
  gldas  = lis%d%lnc*lis%d%lnr
  do i=1,necmwf
     if ( bergdrv%remask1d(i) > 0 ) then
        lb(i)=.true.
     else
        lb(i)=.false.
     end if
  enddo
!<debug>
!  tempecmwf(:,:,1) = -9999.0
!  v = 0
!  do r=1,bergdrv%nrold
!     do c=1,bergdrv%ncold
!        if(lb(c+v)) then 
!           tempecmwf(c,r,1) = 1
!        endif
!     enddo
!     v= v + bergdrv%ncold
!  enddo
!  open(112,file="forcing.bin",form='unformatted')
!  write(112) tempecmwf(:,:,1)
!  stop
!</debug>
  lo = .false.
  tmask = .false.

  do v=1,NF
     if (v .ne. 6) then ! not the v-wind component, which is set to zero.
     
      !=== Transferring current data to 1-D array for interpolation
      c=0
      do i=1,bergdrv%nrold
          do j=1,bergdrv%ncold
             c = c + 1
             f(c) = temp2ecmwf(j,i,v)
!             if(v.eq.1 .and. f(c).ne.-9999.0) print*, c,f(c)
          enddo
      enddo

      !=== Interpolate if forcing and model grids are not both half deg.
!      if (lis%d%domain .ne. 3) then
      if (lis%d%gridDesc(9).ne.0.5) then
          ip = 0
          do ii=1,20
             ipopt(ii) = 0
          enddo
          km = 1
          ibi = 1
          call bilinear_interp(gridDesco,ibi,lb,f,ibo,lo,go,mi,mo, & 
               rlat0,rlon0,w110,w120,w210,w220, & 
               n110,n120,n210,n220,iret)
      else ! forcing and model grids both half degree
	kk = 0
	do r=1,lis%d%lnr
	  y = r + (lis%d%gridDesc(4) + 89.750) / 0.500
	  do c=1,lis%d%lnc
	    x = c + (lis%d%gridDesc(5) + 179.750) / 0.500
	    kk = kk + 1
            eindex = ((y - 1) * 720) + x
	    go(kk) = f(eindex)
	    lo(kk) = lb(eindex)
	  end do
	end do

      end if ! lis%d%domain 
	    
       !=== Convert data to original 3D array & a 2D array to 
       !=== fill in of missing points due to geography difference  
!      tg = -9999.0
      c = 0
      do j = 1, lis%d%lnr
         do i = 1, lis%d%lnc
!            print*, i,j,gindex(i,j),lo(i+c)
            if(gindex(i,j).ne.-1) then 
               geogmask(i,j) = lo(i+c)
               tg(i,j) = go(i+c)
            endif
         enddo
         c = c + lis%d%lnc
      enddo
!      if(v==3) then 
!         open(112,file="forcing.bin",form='unformatted')
!         write(112) tg
!         stop
!      endif
      if(v==2) then 
!         do j=1,lis%d%lnr
!            do i=1, lis%d%lnc
!               print*, 'before ',tg(14,1)
!            enddo
!         enddo
      endif

      call geogfill2(lis%d%lnc,lis%d%lnr,geogmask,tg,v,tmask)
      !       print*, gindex(21,3),tg(21,3),v
#if 0 
      if(v==3) then 
         do j=1,lis%d%lnr
            do i=1, lis%d%lnc
               print*, 'after ',i,j, tg(i,j)
            enddo
         enddo
      endif
#endif
      do j = 1, lis%d%lnr
         do i = 1, lis%d%lnc
            if(gindex(i,j).ne.-1) then 
               if (tg(i,j) .eq. -1.0) then
                  print*, 'No nearest neighbours, v, i, j',v,i,j
                  stop
               endif
               templdas(gindex(i,j),v) = tg(i,j)
            endif
         end do !c
       enddo ! r

    else ! v==6, v-wind component, always zero
       templdas(:,6) = 0.0
    endif
    
 enddo !v

  do v= 1,NF
    !=== Fill in undefined and ocean points
      do c = 1,lis%d%ngrid
        if(order.eq.1)then
           glbdata1(v,c)=templdas(c,v)
        else
           glbdata2(v,c)=templdas(c,v)
        endif
    enddo !r
  enddo !v
  deallocate(templdas)

end subroutine retberg


!BOP
! 
! !ROUTINE: berggrid_2_gldasgrid
! 
! !DESCRIPTION:
! Changes grid_data from ECMWF data convention to GLDAS convention
!
! ECMWF: North-to-South around Greenwich Meridian
! Global grid. Data are written as flat binary from "upper left to 
! lower right" starting at 0.5-degree grid point center coordinates: 
! 0.25E,89.75N and going to 0.25W,89.75S. Here is the write statement:
!
! do i = 1,360   
!   write(14) (val(j,i),j=1,720)
! end do
!
! GLDAS: South-to-North around Date Line
! Full global grid.  Starts at the southernmost latitude and date line, 
! going east and then north.
!
! !REVISION HISTORY:
!  10 Apr 2002: Urszula Jambor;  Code adapted from 
!               ecmwfgrid_2_grid2catgrid, by R. Reichle
! !INTERFACE:
subroutine berggrid_2_gldasgrid( nx, ny, grid_data )
!EOP  
  implicit none
  
  integer, intent(in) :: nx, ny
  real, intent(inout), dimension(nx,ny) :: grid_data
  
  integer :: i, j, m, n
  real :: tmp, tmp_data1(nx), tmp_data2(nx)
  
  ! ------------------------------------------------------------------
  ! some checks
  
  if ((nx /= 720) .or. (ny /= 360)) then
     write (*,*) 'Recmwfgrid_2_gldasgrid(): This routine has only been'
     write (*,*) 'checked for nx=720 and ny=360. Make sure you know'
     write (*,*) 'what you are doing. STOPPING.'
     stop
  end if  
  if ((mod(nx,2) /= 0) .or. (mod(ny,2) /= 0)) then
     write (*,*) 'Recmwfgrid_2_gldasgrid(): This routine can only work'
     write (*,*) 'for even nx and ny. Make sure you know'
     write (*,*) 'what you are doing. STOPPING.'
     stop
  end if
  
  !-------------------------------------------------------------------
  
  do j=1,ny/2
     
     ! swap latitude bands (North-to-South becomes South-to-North)
     n = ny-j+1
     tmp_data1      = grid_data(:,j)
     tmp_data2      = grid_data(:,n)
     
     do i=1,nx/2

        ! shift longitudes (wrapping around Greenwhich Meridian becomes
        !  wrapping around Date Line)
        m = i + nx/2
        tmp          = tmp_data1(i)
        tmp_data1(i) = tmp_data1(m)
        tmp_data1(m) = tmp
        
        tmp          = tmp_data2(i)
        tmp_data2(i) = tmp_data2(m)
        tmp_data2(m) = tmp
        
     end do
     
     grid_data(:,j) = tmp_data2
     grid_data(:,n) = tmp_data1
     
  end do

end subroutine Berggrid_2_gldasgrid
!EOC

!BOP
! 
! !ROUTINE: fillgaps
! 
! !DESCRIPTION:
! Fills in values for NSIPP tilespace land points
! where no ECMWF reanalysis data is available via GEOGFILL 
! by assigning most appropriate land-point value along
! the latitudinal circle of original tilespace point.
! Developed manually with 17 points in mind.
!
! !REVISION HISTORY:
!  23 Jul 2002: Urszula Jambor
!  03 Sep 2002: Urszula Jambor, revised points to reflect 
!               NSIPP land mask correction.
! !INTERFACE:
subroutine fillgaps( nc, nr, v, arr )
!EOP
  implicit none

  integer :: nc, nr, v
  real :: arr(nc,nr)
!BOC
  arr( 55, 1) = arr( 59, 4)
  arr( 62, 2) = arr( 59, 4)
  arr( 88, 8) = arr( 46, 8)
  arr( 94, 8) = arr( 46, 8)
  arr(  2, 9) = arr(142, 9)
  arr(  3,20) = arr(133,20)
  arr(139,20) = arr(133,20)
  arr(140,20) = arr(133,20)
  arr( 95,29) = arr( 89,29)
  arr( 37,30) = arr( 41,30)
  arr( 36,31) = arr( 41,31)
  arr( 37,31) = arr( 41,31)
  arr(136,34) = arr(123,34)
  arr( 62,50) = arr( 69,50)
  arr( 63,50) = arr( 69,50)
  arr(142,57) = arr(144,57)
  arr( 70,67) = arr( 83,67)

  if (v == 3) then
    arr( 88, 8) = arr( 84,14)
    arr( 94, 8) = arr( 84,14)
  endif

end subroutine fillgaps
!EOP


