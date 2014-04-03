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
! ret_reanlecmwf.F90:
!
! DESCRIPTION:
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
! REVISION HISTORY:
!  09 Apr 2002: Urszula Jambor; Original code, based on readgeos.f
!  20 Dec 2002: Matt Rodell; Don't ipolate if running 1/2 deg, also fix
!               error in setting v-wind to zero
!  25 Mar 2003: Urszula Jambor; Modified argument list passed to GEOGFILL2.
!=========================================================================

subroutine ret_reanlecmwf(order, ldas, grid, yr, mo, da, hr, ferror)

  use ldas_module      ! LDAS non-model-specific 1-D variables
  use grid_module      ! LDAS non-model-specific grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc,ldas%nr)

  !=== Local Variables ==================================================

  !Subroutine inputs
  integer :: order            ! lower(1) or upper(2) time interval bdry
  integer :: yr,mo,da,hr      ! data and hour (multiple of 6)
  integer :: ferror           ! set to zero if there's an error

  ! Specify Reanalysis ECMWF forcing parameters & file paths
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
       '/NILE/ECMWF/GRID/     ',  &
       '/GLDAS8/ECMWF/GRID/   ',  &
       '/NSIPP/ECMWF/GRID/    ',  &
       '/NSIPP/ECMWF/GRID/    ',  &
       '/GLDAS8/ECMWF/GRID/   ',  &
       '/GLDAS8/ECMWF/GRID/   ',  &
       '/GLDAS8/ECMWF/GRID/   ',  &
       '/GLDAS8/ECMWF/GRID/   '   /)

  integer, dimension(N_EF), parameter :: fnum = (/ &
       31, 33, 34, 35, 36, 37, 38, 39 /) 

  character(200) :: infile

  real :: tempecmwf(ldas%ncold,ldas%nrold,N_EF)
  real :: temp2ecmwf(ldas%ncold,ldas%nrold,NF)
  real :: templdas(ldas%nc,ldas%nr,NF)
  real :: p_kPa(ldas%ncold,ldas%nrold)
  real :: outdata(ldas%ncold,ldas%nrold)
  real :: moredata(ldas%nc,ldas%nr)
  character :: cyr*4, cmo*2 
  integer :: timestep
  integer :: i, j, t, v, ii, jj
  integer :: ios,ioerror      ! set to non-zero if there's an error
  integer :: istat

  integer :: gldas,necmwf     ! Size of I/O 1D fields
  integer :: lmask(ldas%nc*ldas%nr)                     ! ID land mask
  real :: f(ldas%ncold*ldas%nrold),go(ldas%nc*ldas%nr)  ! 1D I/O fields
  real :: rlat(ldas%nc*ldas%nr)
  real :: rlon(ldas%nc*ldas%nr)         ! Output lat & lon
  real :: tg(ldas%nc,ldas%nr)           ! Interpolated 2D data field
  integer :: ibi,no,ibo,ipopt(20)
  integer :: iret,km,ip,c
  integer :: kgds(200),kgdso(200)       ! Input,output grid info arrays
  logical*1 :: lb(ldas%ncold*ldas%nrold)
  logical*1 :: lo(ldas%nc*ldas%nr)      ! Input and output bitmaps
  logical*1 :: geogmask(ldas%nc,ldas%nr)! 2D output bitmap
  logical*1 :: tmask(LDAS%NC,LDAS%NR)   ! 2d valid temperature bitmap

  !=== End Defining Local Variables =====================================

  ! if a problem, ferror is set to zero
  ferror = 1

  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour correspond to
  ! and which file to open.
  !----------------------------------------------------------------

  timestep = 4*(da - 1) + (1 + hr/6)
  print*, 'order and month-timestep', order, timestep
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mo

  !=== Open ECMWF forcing files

  if ( (ldas%tscount==1).and.(order==1) ) then  !beginning of run

     do v = 1, N_EF

        ! File name for data variable(v)/year/yearmo
        infile=trim(v_inpath(v))//trim(ecmwf_fv(v))//'/'//cyr//'/'//cyr//cmo
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
              do i = 1, ldas%nrold
                 read(fnum(v)) (tempecmwf(j,i,v), j=1, ldas%ncold)
              end do
           end do
        end if

     end do !v

  else if ( timestep==1 ) then !close previous month files & open current

     do v = 1, N_EF

        close(fnum(v))
        ! New file name for data variable(v)/year/yearmo
        infile=trim(v_inpath(v))//trim(ecmwf_fv(v))//'/'//cyr//'/'//cyr//cmo
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

  do v = 1, N_EF

     !----------------------------------------------------------------
     ! Extract grid array for chosen time
     !----------------------------------------------------------------
     do i = 1, ldas%nrold
        read(fnum(v)) (tempecmwf(j,i,v), j=1, ldas%ncold)
     end do

     !----------------------------------------------------------------
     ! Change data from ECMWF grid convention to GLDAS one 
     ! Swap latitude band and shift longitudes by 180deg. 
     !----------------------------------------------------------------

     call Recmwfgrid_2_gldasgrid(ldas%ncold,ldas%nrold,tempecmwf(:,:,v))

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

  outdata = temp2ecmwf(:,:,2)
  open(22, file='aaronTd.bin', form='unformatted')
  write(22) outdata
  close(22) 

  !-----------------------------------------------------------------
  ! Interpolate each forcing variable to GLDAS domain
  !-----------------------------------------------------------------

  !=== Initialize input & output grid arrays
  kgds  = 0
  kgdso = 0
           
  !=== Set input & output grid array values (reanlECMWF to GLDAS)
  kgds(1) = 0
  kgds(2) = ldas%ncold
  kgds(3) = ldas%nrold
  kgds(4) = -89750
  kgds(5) = -179750
  kgds(6) = 128
  kgds(7) = 89750
  kgds(8) = 179750
  kgds(9) = 500
  kgds(10) = 500
  kgds(20) = 255
        
  kgdso = ldas%ldas_kgds

  !=== Define input & output data bitmaps
  necmwf = ldas%ncold*ldas%nrold
  gldas  = ldas%nc*ldas%nr

  do i=1,necmwf
     if ( ldas%remask1d(i) > 0 ) then
        lb(i)=.true.
     else
        lb(i)=.false.
     end if
  enddo
  lo = .true.
  tmask = .false.

  do v=1,NF

     if (v .ne. 6) then ! not the v-wind component, which is set to zero.
     
      !=== Interpolate if forcing and model grids are not both half deg.
      if (ldas%domain .ne. 5) then
     
       !=== Transferring current data to 1-D array for interpolation
       c=0
       do i=1,ldas%nrold
          do j=1,ldas%ncold
             c = c + 1
             f(c) = temp2ecmwf(j,i,v)
          enddo
       enddo

       !=== Setting interpolation options 
       !=== (ip=0, bilinear),(iopt=0, no options), 
       !=== (km=1, one variable),(ibi=1, use bitmap)
       !=== Adjust to budget-bilinear for precip forcing fields
       if (v .eq. 8 .or. v .eq. 9) then
          ip = 3
          ipopt(1) = -1
          ipopt(2) = -1
          km = 1
          ibi = 1
       else                 
          ip = 0
          do ii=1,20
             ipopt(ii) = 0
          enddo
          km = 1
          ibi = 1
       endif

       !=== Interpolate data from ECMWF grid to GLDAS grid
       CALL IPOLATES(ip,ipopt,kgds,kgdso,necmwf,gldas,km,ibi,lb, &
                     f,no,rlat,rlon,ibo,lo,go,iret)
     
       !=== Convert data to original 3D array & a 2D array to 
       !=== fill in of missing points due to geography difference  
       c = 0
       do j = 1, ldas%nr
          do i = 1, ldas%nc
             geogmask(i,j) = lo(i+c)
             tg(i,j) = go(i+c)
          enddo
         c = c + ldas%nc
       enddo

       CALL GEOGFILL2(ldas%nc,ldas%nr,grid%fimask,geogmask,tg,v,tmask)
       if (ldas%koster==1) call fillgaps(ldas%nc,ldas%nr,v,tg)

       do j = 1, ldas%nr
         do i = 1, ldas%nc
            if (tg(i,j) .eq. -1.0) then
		print*, 'No nearest neighbours, v, i, j',v,i,j
                stop
            endif
            templdas(i,j,v) = tg(i,j)
         end do !c
       enddo ! r

      else ! forcing and model grids both half degree
       templdas(:,:,v) = temp2ecmwf(:,:,v)
      end if ! ldas%domain 

     else ! v==6, v-wind component, always zero
!        templdas(i,j,6) = 0.0  THIS WAS AN ERROR, FIXED 12/20/02
        templdas(:,:,6) = 0.0
     endif

   enddo !v

  do v= 1,NF
    !=== Fill in undefined and ocean points
    do j = 1, ldas%nr
      do i = 1, ldas%nc

        if (grid(i,j)%mask < 1.0) then
            templdas(i,j,v) = ldas%udef
        endif

        if(order.eq.1)then
           grid(i,j)%glbdata1(v)=templdas(i,j,v)
        else
           grid(i,j)%glbdata2(v)=templdas(i,j,v)
        endif

      enddo !c
    enddo !r
  enddo !v

!  moredata = 0.0
!  moredata = templdas(:,:,3)
!  open(99, file='currentP.mask', form='unformatted')
!  write(99) moredata
!  close(99)

end subroutine ret_reanlecmwf


!=================================================================
! S U B R O U T I N E S
!=================================================================

!=================================================================
! DESCRIPTION:
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
!  REVISION HISTORY:
!  10 Apr 2002: Urszula Jambor;  Code adapted from 
!               ecmwfgrid_2_grid2catgrid, by R. Reichle
!=================================================================

subroutine Recmwfgrid_2_gldasgrid( nx, ny, grid_data )
  
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

end subroutine Recmwfgrid_2_gldasgrid



!=================================================================
! DESCRIPTION:
! Fills in values for NSIPP tilespace land points
! where no ECMWF reanalysis data is available via GEOGFILL 
! by assigning most appropriate land-point value along
! the latitudinal circle of original tilespace point.
! Developed manually with 17 points in mind.
!
!  REVISION HISTORY:
!  23 Jul 2002: Urszula Jambor
!  03 Sep 2002: Urszula Jambor, revised points to reflect 
!               NSIPP land mask correction.
!=================================================================


subroutine fillgaps( nc, nr, v, arr )

  implicit none

  integer :: nc, nr, v
  real :: arr(nc,nr)

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


!========EOF=========================================================

