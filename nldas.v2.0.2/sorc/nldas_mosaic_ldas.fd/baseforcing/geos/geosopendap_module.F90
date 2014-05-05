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
! !MODULE: geosopendap_module.F90
! 
! !DESCRIPTION: 
!
!  This module contains routines needed to initialize and control variables
! required for the execution of GDS-based I/O specfific to GEOS forcing
! routines
!
! !REVISION HISTORY:
! 10 Jul 2003; James Geiger  Initial Specification 
! 22 Dec 2003; Sujay Kumar Separted geos specific routines from the main
!              module
!
! !INTERFACE:
module geosopendap_module

#if ( defined OPENDAP )
! !USES:
  use geosdomain_module, only : geosdrv
#if ( defined ESMF_TIMEMANAGER_KLUGE )
  use ESMF_TimeMgmtMod
#else
  use time_manager
#endif
  use lisdrv_module,     only : lis, grid, gindex
  use grid_spmdMod,      only : gdi, gdisp
  use tile_spmdMod,      only : di_array
  use spmdMod
  use opendap_module,    only : opendap_data_prefix, ciam

  implicit none

#if ( defined ESMF_TIMEMANAGER_KLUGE )

!   type(esmf_date), save :: geos_ref_date !Reference date for GEOS forcing data

#else

   type date_struct
      integer :: year
      integer :: month
      integer :: day
      integer :: hours
      integer :: minutes
      integer :: seconds
   end type date_struct

   type(date_struct), save :: geos_ref_date !Reference date for GEOS forcing data
#endif
   integer               :: geos_slat     !Southern latitude boundary for each subdomain (for GEOS)
   integer               :: geos_nlat     !Northern latitude boundary for each subdomain (for GEOS)
   integer               :: geos_nc       !Number of columns for each subdomain (for GEOS)
   integer               :: geos_nr       !Number of rows for each subdomain (for GEOS)
   real                  :: input_slat    !Southern latitude boundary for the native domain
   real                  :: input_nlat    !Northern latitude boundary for the native domain
   character*4           :: cgeos_slat    !character representation of geos$_-$slat
   character*4           :: cgeos_nlat    !character representation of geos$_-$nlat
   integer               :: grid_offset   !Global to local grid mapping offset
   integer               :: fnroffset

   contains
!BOP
! !ROUTINE: opendap_geos_init
! 
! !DESCRIPTION:
!  Initializes the GEOS-GDS variables
!
! !INTERFACE:
   subroutine opendap_geos_init()
!EOP
     implicit none
!BOC       
     call init_geos_vars()
!     call reset_geos_filepaths()
!EOC       
   end subroutine opendap_geos_init

!BOP
! !ROUTINE: reset_geos_filepaths
! 
! !DESCRIPTION:
!  Resets input data filenames for GEOS forcing for execution through GDS
!
! !INTERFACE:
     subroutine reset_geos_filepaths()
!EOP
       implicit none
!BOC
       geosdrv%geosdir      = trim(opendap_data_prefix)//'/'// &
                              trim(adjustl(ciam))//'/'//geosdrv%geosdir
!EOC
     end subroutine reset_geos_filepaths
!BOP
! !ROUTINE: init_geos_vars
! 
! !DESCRIPTION:
! 
! Computes domain decomposition for native as well as 
! interpolated domains for input GEOS forcing
! 
! !INTERFACE:
    subroutine init_geos_vars()
!EOP
      use lisdrv_module, only: tile
      implicit none

      integer :: res, center, bottom, domain
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
         print*,'Error!  Cannot handle nldas.'
         stop 999
      case(2)
         print*,'Error!  Cannot handle 2x2.5.'
         stop 999
      case(3) ! 1 deg
         res = 1000
         center = 500
         bottom = -59500
      case(4) ! 1/2 deg
         res = 500
         center = 750
         bottom = -59750
      case(5) ! 1/4 deg
         res = 250
         center = 875
         bottom = -59875
      case(6) ! 5km
         res = 50
         center = 975
         bottom = -59975
      case(7) ! 5km
         res = 50
         center = 975
         bottom = -59975
      case(8) ! 1km
         res = 10
         center = 995
         bottom = -59995
      case DEFAULT
         print*, "Select domain size (1,2,3,4,5,6,7,8)"
         stop 999
      end select
      
      call set_geos_lat(geos_slat,geos_nlat,input_slat,input_nlat)
      
      lis%d%ngrid = gdi(iam)
      lis%d%nch   = di_array(iam)
      
      geos_nc    = 360
      geos_nr    = ( geos_nlat - geos_slat + 1 )

      fnroffset = geos_slat - 1
      grid_offset = tile(1)%index-1

      write(cgeos_slat, '(i4)') geos_slat
      write(cgeos_nlat, '(i4)') geos_nlat

      print*,'DBG: geos_init -- geos_slat', geos_slat, ' (', iam, ')'
      print*,'DBG: geos_init -- geos_nlat', geos_nlat, ' (', iam, ')'
      print*,'DBG: geos_init -- ngrid', lis%d%ngrid, ' (', iam, ')'
      print*,'DBG: geos_init -- glbngrid', lis%d%glbngrid, ' (', iam, ')'
      print*,'DBG: geos_init -- ncold', geosdrv%ncold, ' (', iam, ')'
      print*,'DBG: geos_init -- nrold', geosdrv%nrold, ' (', iam, ')'
      print*,'DBG: geos_init -- lnc', lis%d%lnc, ' (', iam, ')'
      print*,'DBG: geos_init -- lnr', lis%d%lnr, ' (', iam, ')'
      print*,'DBG: geos_init -- fnroffset', fnroffset, ' (', iam, ')'
      print*,'DBG: geos_init -- grid_offset', grid_offset, ' (', iam, ')'
      print*,'DBG: geos_init -- cgeos_slat ', cgeos_slat, ' (', iam, ')'
      print*,'DBG: geos_init -- cgeos_nlat ', cgeos_nlat, ' (', iam, ')'
      print*,'DBG: geos_init -- geos_nc', geos_nc, ' (', iam, ')'
      print*,'DBG: geos_init -- geos_nr', geos_nr, ' (', iam, ')'
!EOC
    end subroutine init_geos_vars

!BOP
! !ROUTINE: def_gridDesc
! 
! !DESCRIPTION:  
! Initializes the grid description array for GEOS-GDS runs. 
! 
! !INTERFACE:
    subroutine def_gridDesc(gridDesci)
!EOP      
      real, intent(out) :: gridDesci(50)
!BOC
       gridDesci(1) = 0
       gridDesci(2) = geos_nc    ! geosdrv%ncold
       gridDesci(3) = geos_nr    ! geosdrv%nrold
       gridDesci(4) = input_slat ! -90.000
       gridDesci(5) = -180.000
       gridDesci(7) = input_nlat ! 90.000
       gridDesci(8) = 179.000
       gridDesci(6) = 128
       gridDesci(9) = 1.000
       gridDesci(10) = 1.000
       gridDesci(20) = 255

!EOC
    end subroutine def_gridDesc

!BOP
! !ROUTINE: init_geos_ref_date
! 
! !DESCRIPTION:  
! Initializes the reference date for GEOS forcing data
! 
! !INTERFACE:
   subroutine init_geos_ref_date()
!EOP
      implicit none

#if ( defined ESMF_TIMEMANAGER_KLUGE )
      integer :: ref_ymd = 20001219 ! 0z19dec2000
      integer :: ref_tod = 0
      integer :: rc
#endif
!BOC
#if ( defined ESMF_TIMEMANAGER_KLUGE )
      geos_ref_date = esmf_dateinit(esmf_no_leap, ref_ymd, ref_tod, rc)
#else
      geos_ref_date%year    = 2000
      geos_ref_date%month   = 12
      geos_ref_date%day     = 19
      geos_ref_date%hours   = 0
      geos_ref_date%minutes = 0
      geos_ref_date%seconds = 0
#endif
!EOC
    end subroutine init_geos_ref_date

!BOP
! !ROUTINE: get_geos_index
!
! !DESCRIPTION: 
! Computes the time-based index for GEOS forcing data
! 
! !INTERFACE:
   function get_geos_index(offset)
     implicit none
! !INPUT PARAMETERS:
     integer, intent(in) :: offset   ! offset from current date in hours
!EOP
     integer :: get_geos_index

     integer :: ndays, nsecs
     logical, save :: ref_data_uninit = .true.

#if ( defined ESMF_TIMEMANAGER_KLUGE )
     type(esmf_date) :: current_date
     type(esmf_time) :: diff
     logical :: islater
     integer :: rc
     integer :: ymd, tod
#endif
!BOC     
     if ( ref_data_uninit ) then
        print*,'DBG: get_geos_index -- initializing ref date',' (', iam, ')'
        call init_geos_ref_date()
        ref_data_uninit = .false.
     endif

#if ( defined ESMF_TIMEMANAGER_KLUGE )
     ymd = ( lis%t%yr * 10000 ) + ( lis%t%mo * 100 ) + lis%t%da
     tod = ( lis%t%hr *  3600 ) + ( lis%t%mn *  60 ) + lis%t%ss
     current_date = esmf_dateinit(esmf_no_leap, ymd, tod, rc)
     
     diff = esmf_timeinit()
     call esmf_datediff(current_date, geos_ref_date, diff, islater, rc)
     call esmf_timeget(diff, ndays, nsecs, rc)

#else

     call diff_date(lis%t%yr, lis%t%mo, lis%t%da, &
                    lis%t%hr, lis%t%mn, lis%t%ss, &
                    geos_ref_date%year,           &
                    geos_ref_date%month,          &
                    geos_ref_date%day,            &
                    geos_ref_date%hours,          &
                    geos_ref_date%minutes,        &
                    geos_ref_date%seconds,        &
                    ndays, nsecs)
#endif

     
     get_geos_index = ( (ndays * 24) + (nsecs / 3600) + offset ) / 3 + 1
     print*,'DBG: get_geos_index -- get_geos_index = ', get_geos_index, &
          ' (', iam, ')'
!EOC     

   end function get_geos_index
!BOP
! !ROUTINE: set_geos_lat
! 
! !DESCRIPTION:
!
! Computes the latitudes of the decomposed domain for GEOS forcing data
!
! !INTERFACE:
   subroutine set_geos_lat(slat, nlat, islat, inlat)
      implicit none
! !OUTPUT PARAMETERS:
      integer, intent(out) :: slat, nlat
      real, intent(out)    :: islat, inlat
!EOP
      real :: lat
!BOC
!--------------------------------------------------------------------
! Set southern latitude index
!--------------------------------------------------------------------
      print*, 'DBG: set_geos_lat -- grid(1)%lat', grid(1)%lat,' (', iam, ')'
      lat = grid(1)%lat
      if ( lat < 0.0 ) then
         slat = int(lat) - 1 ! E.g. map -54.5 \to -55
      else
         slat = int(lat)     ! E.g. map  40.5 \to  40
      endif
      if ( slat < -90 ) then
         print*, 'ERR: set_geos_lat -- Setting slat = -90', &
                 ' (', iam, ')'
         slat = -90
      endif
!--------------------------------------------------------------------
! Set southern latitude boundary
!--------------------------------------------------------------------
      islat = 1000 * slat
!--------------------------------------------------------------------
! Set northern latitude index
!--------------------------------------------------------------------
      print*, 'DBG: set_geos_lat -- grid(gdi(iam))', &
              gdi(iam), grid(gdi(iam))%lat,' (', iam, ')'
      lat = grid(gdi(iam))%lat
      if ( lat < 0.0 ) then
         nlat = int(lat)     ! E.g. map -10.5 \to -10
      else
         nlat = int(lat) + 1 ! E.g. map  10.5 \to  11
      endif
      if ( nlat > 90 ) then
         print*, 'ERR: set_geos_lat -- Setting nlat = 90', &
                 ' (', iam, ')'
         nlat = 90
      endif
!--------------------------------------------------------------------
! Set northern latitude boundary
!--------------------------------------------------------------------
      inlat = 1000 * nlat

      slat = slat + 90 + 1 ! map [-90,90] \to [1,181]
      nlat = nlat + 90 + 1 ! map [-90,90] \to [1,181]
      slat = 1
      nlat = 181
      !islat = -90000
      !inlat = 90000
      islat = -90.000
      inlat = 90.000
!EOC
    end subroutine set_geos_lat
#endif

  end module geosopendap_module
