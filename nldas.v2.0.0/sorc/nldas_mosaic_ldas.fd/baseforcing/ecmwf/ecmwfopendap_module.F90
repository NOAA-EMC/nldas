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
! !MODULE: ecmwfopendap_module.F90
! 
! !DESCRIPTION: 
!
!  This module contains routines needed to initialize and control variables
! required for the execution of GDS-based I/O specfific to ECMWF forcing
! routines
!
! !REVISION HISTORY:
! 10 Jul 2003; James Geiger  Initial Specification 
! 22 Dec 2003; Sujay Kumar Separted ECMWF specific routines from the main
!              module
!
! !INTERFACE:
module ecmwfopendap_module
#if ( defined OPENDAP )
! TODO: Stubs to allow the code to compile
! TODO: ECMWF does not yet work in GDS-mode
   contains
   subroutine opendap_ecmwf_init()
   end subroutine opendap_ecmwf_init
   subroutine def_gridDesc(gridDesci)
      real :: gridDesci(50)
   end subroutine def_gridDesc
#endif
#if ( defined OPENDAP_KLUGE )
! TODO: Correct this section
! !USES:
  use ecmwfdomain_module, only : ecmwfdrv
  use ESMF_TimeMgmtMod
  use lisdrv_module, only : lis, grid, gindex
  use grid_spmdMod,   only : gdi, gdisp
  use tile_spmdMod,   only : di_array
  use spmdMod

  implicit none
   type(esmf_date), save :: ecmwf_ref_date !Reference date for ECMWF forcing data
   integer               :: ecmwf_slat     !Southern latitude boundary for each subdomain (for ECMWF)
   integer               :: ecmwf_nlat     !Northern latitude boundary for each subdomain (for ECMWF)
   integer               :: ecmwf_nc       !Number of columns for each subdomain (for ECMWF)
   integer               :: ecmwf_nr       !Number of rows for each subdomain (for ECMWF)
   integer               :: input_slat     !Southern latitude boundary for the native domain
   integer               :: input_nlat     !Northern latitude boundary for the native domain
   character*4           :: cecmwf_slat    !character representation of ecmwf$_-$slat
   character*4           :: cecmwf_nlat    !character representation of ecmwf$_-$nlat
   character*3           :: ciam           !character representation of processor id
   integer               :: grid_offset    !Global to local grid mapping offset
   integer               :: fnroffset

#if ( defined ABSOFT )
   character*15 :: opendap_ecmwf_home = '/home/jim/'
#else
   character*2  :: opendap_ecmwf_home = './'
#endif
   contains
!BOP
! !ROUTINE: opendap_ecmwf_init
! 
! !DESCRIPTION:
!  Initializes the ECMWF-GDS variables
!
! !INTERFACE:
   subroutine opendap_ecmwf_init()
!EOP
     implicit none
!BOC       
     call init_ecmwf_vars()
     call reset_ecmwf_filepaths()
!EOC       
   end subroutine opendap_ecmwf_init

!BOP
! !ROUTINE: reset_ecmwf_filepaths
! 
! !DESCRIPTION:
!  Resets input data filenames for ECMWF forcing for execution through GDS
!
! !INTERFACE:
     subroutine reset_ecmwf_filepaths()
!EOPI
       implicit none
!BOC
       ecmwfdrv%ecmwfdir = trim(opendap_ecmwf_home)//trim(adjustl(ciam))//'/'//ecmwfdrv%ecmwfdir
!EOC
     end subroutine reset_ecmwf_filepaths
!BOP
! !ROUTINE: init_ecmwf_vars
! 
! !DESCRIPTION:
! 
! Computes domain decomposition for native as well as 
! interpolated domains for input ECMWF forcing
! 
! !INTERFACE:
    subroutine init_ecmwf_vars()
!EOPI
      use lisdrv_module, only: tile
      implicit none

      integer :: res, center, bottom, domain
      character*5 :: line
!BOC
      if(lis%d%gridDesc(9) .eq. 0.01) then 
         domain = 8 
      elseif(lis%d%gridDesc(9).eq.0.05) then 
         domain = 7
      elseif(lis%d%gridDesc(9) .eq. 0.125) then 
         domain = 6
      elseif(lis%d%gridDesc(9) .eq. 0.25) then 
         domain = 5
      elseif(lis%d%gridDesc(9) .eq. 0.50) then 
         domain = 4
      elseif(lis%d%gridDesc(9) .eq. 1.0) then      
         domain = 3
      elseif((lis%d%gridDesc(9) .eq. 2) .and.  &
           (lis%d%gridDesc(10) .eq. 2.5)) then 
         domain = 2
      endif
      select case (domain)
      case(1)
         print*,'Error!  Cannot handle ndas.'
         stop 999
      case(2)
         res = 250
         center = 875
         bottom = -59875
      case(3)
         print*,'Error!  Cannot handle 2x2.5.'
         stop 999
      case(4)
         res = 1000
         center = 500
         bottom = -59500
      case(5)
         res = 500
         center = 750
         bottom = -59750
      case(6)
         res = 50
         center = 975
         bottom = -59975
      case DEFAULT
         print*, "Select domain size (1,2,3,4,5,6)"
         stop 999
      end select
      
      call set_ecmwf_lat(ecmwf_slat,ecmwf_nlat,input_slat,input_nlat)
      
      lis%d%ngrid = gdi(iam)
      lis%d%nch   = di_array(iam)
      
      ecmwf_nc    = 1440
      ecmwf_nr    = ( ecmwf_nlat - ecmwf_slat + 1 )

      fnroffset = ecmwf_slat - 1
      grid_offset = tile(1)%index-1

      write(ciam, '(i3)') iam
      write(cecmwf_slat, '(i4)') ecmwf_slat
      write(cecmwf_nlat, '(i4)') ecmwf_nlat

      print*,'DBG: ecmwf_init -- ecmwf_slat', ecmwf_slat, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- ecmwf_nlat', ecmwf_nlat, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- ngrid', lis%d%ngrid, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- glbngrid', lis%d%glbngrid, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- ncold', ecmwfdrv%ncold, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- nrold', ecmwfdrv%nrold, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- lnc', lis%d%lnc, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- lnr', lis%d%lnr, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- fnroffset', fnroffset, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- grid_offset', grid_offset, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- cecmwf_slat ', cecmwf_slat, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- cecmwf_nlat ', cecmwf_nlat, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- ciam ', ciam, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- ecmwf_nc', ecmwf_nc, ' (', iam, ')'
      print*,'DBG: ecmwf_init -- ecmwf_nr', ecmwf_nr, ' (', iam, ')'
      line = '"'//ciam//'"'
      print*,'DBG: ecmwf_init -- line', line, ' (', iam, ')'
!EOC
    end subroutine init_ecmwf_vars
!BOP
! !ROUTINE: def_gridDesc
! 
! !DESCRIPTION:  
! Initializes the gridDesc array for ECMWF-GDS runs. 
! 
! !INTERFACE:
    subroutine def_gridDesc(gridDesci)
!EOP      
      real :: gridDesci(50)
!BOC
      gridDesci(1) = 0
      gridDesci(2) = ecmwf_nc
      gridDesci(3) = ecmwf_nr
      gridDesci(4) = input_slat
      gridDesci(7) = input_nlat
      gridDesci(5) = -180.000
      gridDesci(6) = 128
      gridDesci(8) = 179.750
      gridDesci(9) = .250
      gridDesci(10) = .250
      gridDesci(20) = 255
!EOC
    end subroutine def_gridDesc
!BOP
! !ROUTINE: init_ecmwf_ref_date
! 
! !DESCRIPTION:  
! Initializes the reference date for ECMWF forcing data
! 
! !INTERFACE:
   subroutine init_ecmwf_ref_date()
!EOP
      implicit none

      integer :: ref_ymd = 20001219 ! 0z19dec2000
      integer :: ref_tod = 0
      integer :: rc
!BOC
      ecmwf_ref_date = esmf_dateinit(esmf_no_leap, ref_ymd, ref_tod, rc)
!EOC
    end subroutine init_ecmwf_ref_date

!BOP
! !ROUTINE: get_ecmwf_index
!
! !DESCRIPTION: 
! Computes the time-based index for ECMWF forcing data
! 
! !INTERFACE:
   function get_ecmwf_index(offset)
     implicit none
! !INPUT PARAMETERS:
     integer, intent(in) :: offset   ! offset from current date in hours
!EOP
     integer :: get_ecmwf_index
     type(esmf_date) :: current_date
     type(esmf_time) :: diff
     logical :: islater
     integer :: rc
     integer :: ymd, tod
     integer :: ndays, nsecs
     logical, save :: ref_data_uninit = .true.
!BOC     
     if ( ref_data_uninit ) then
        print*,'DBG: get_ecmwf_index -- initializing ref date',' (', iam, ')'
        call init_ecmwf_ref_date()
        ref_data_uninit = .false.
     endif
     ymd = ( lis%t%yr * 10000 ) + ( lis%t%mo * 100 ) + lis%t%da
     tod = ( lis%t%hr *  3600 ) + ( lis%t%mn *  60 ) + lis%t%ss
     current_date = esmf_dateinit(esmf_no_leap, ymd, tod, rc)
     
     diff = esmf_timeinit()
     call esmf_datediff(current_date, ecmwf_ref_date, diff, islater, rc)
     call esmf_timeget(diff, ndays, nsecs, rc)
     
     get_ecmwf_index = ( (ndays * 24) + (nsecs / 3600) + offset ) / 3 + 1
     print*,'DBG: get_ecmwf_index -- get_ecmwf_index = ', get_ecmwf_index, &
          ' (', iam, ')'
!EOC     
   end function get_ecmwf_index
!BOP
! !ROUTINE: set_ecmwf_lat
! 
! !DESCRIPTION:
!
! Computes the latitudes of the decomposed domain for ECMWF forcing data
!
! !INTERFACE:
   subroutine set_ecmwf_lat(slat, nlat, islat, inlat)
      implicit none
! !OUTPUT PARAMETERS:
      integer, intent(out) :: slat, nlat, islat, inlat
!EOP
      real :: lat
!BOC
!--------------------------------------------------------------------
! Set southern latitude index
!--------------------------------------------------------------------
      print*, 'DBG: set_ecmwf_lat -- grid(1)%lat', grid(1)%lat,' (', iam, ')'
      lat = grid(1)%lat
      if ( lat < 0.0 ) then
         slat = int(lat) - 1 ! E.g. map -54.5 \to -55
      else
         slat = int(lat)     ! E.g. map  40.5 \to  40
      endif
      if ( slat < -90 ) then
         print*, 'ERR: set_ecmwf_lat -- Setting slat = -90', &
                 ' (', iam, ')'
         slat = -90
      endif
!--------------------------------------------------------------------
! Set southern latitude boundary
!--------------------------------------------------------------------
      islat = 250 * slat
!--------------------------------------------------------------------
! Set northern latitude index
!--------------------------------------------------------------------
      print*, 'DBG: set_ecmwf_lat -- grid(gdi(iam))', &
              gdi(iam), grid(gdi(iam))%lat,' (', iam, ')'
      lat = grid(gdi(iam))%lat
      if ( lat < 0.0 ) then
         nlat = int(lat)     ! E.g. map -10.5 \to -10
      else
         nlat = int(lat) + 1 ! E.g. map  10.5 \to  11
      endif
      if ( nlat > 90 ) then
         print*, 'ERR: set_ecmwf_lat -- Setting nlat = 90', &
                 ' (', iam, ')'
         nlat = 90
      endif
!--------------------------------------------------------------------
! Set northern latitude boundary
!--------------------------------------------------------------------
      inlat = 250 * nlat

      slat = slat + 90 + 0.25 ! map [-90,90] \to [1,181]
      nlat = nlat + 90 + 0.25 ! map [-90,90] \to [1,181]
      slat = 1
      nlat = 601
      islat = -90000
      inlat = 90000
!EOC
    end subroutine set_ecmwf_lat
#endif

  end module ecmwfopendap_module
