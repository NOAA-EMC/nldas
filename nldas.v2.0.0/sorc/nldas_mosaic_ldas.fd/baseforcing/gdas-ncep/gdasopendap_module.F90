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
! !MODULE: gdasopendap_module.F90
! 
! !DESCRIPTION: 
!
!  This module contains routines needed to initialize and control variables
! required for the execution of GDS-based I/O specfific to GDAS forcing
! routines
!
! !REVISION HISTORY:
! 05 Feb 2004; James Geiger  Initial Specification 
!
! !INTERFACE:
module gdasopendap_module
#if ( defined OPENDAP )
! !USES:
  use lisdrv_module, only : lis, grid
  use grid_spmdMod,   only : gdi, gdisp
  use tile_spmdMod,   only : di_array
  use spmdMod

  implicit none
   integer               :: gdas_slat     !Southern latitude boundary for each subdomain (for GDAS)
   integer               :: gdas_nlat     !Northern latitude boundary for each subdomain (for GDAS)
   integer               :: gdas_wlon     !Western longitude boundary for each subdomain (for GDAS)
   integer               :: gdas_elon     !Eastern longitude boundary for each subdomain (for GDAS)
   integer               :: gdas_nc       !Number of columns for each subdomain (for GDAS)
   integer               :: gdas_nr       !Number of rows for each subdomain (for GDAS)
   integer               :: input_slat    !Southern latitude boundary for the native domain
   integer               :: input_nlat    !Northern latitude boundary for the native domain
   character*4           :: cgdas_slat    !character representation of gdas$_-$slat
   character*4           :: cgdas_nlat    !character representation of gdas$_-$nlat
   character*4           :: cgdas_wlon    !character representation of gdas$_-$wlon
   character*4           :: cgdas_elon    !character representation of gdas$_-$elon
   integer               :: grid_offset   !Global to local grid mapping offset
   integer               :: fnroffset

   contains
!BOP
! !ROUTINE: opendap_gdas_init
! 
! !DESCRIPTION:
!  Initializes the GDAS-GDS variables
!
! !INTERFACE:
   subroutine opendap_gdas_init(gdasdrv)
!EOP
     use gdasdrv_module
     implicit none
     type(gdasdrvdec)      :: gdasdrv
!BOC       
     call init_gdas_vars()
     call reset_gdas_filepaths(gdasdrv)
!EOC       
   end subroutine opendap_gdas_init

!BOP
! !ROUTINE: reset_gdas_filepaths
! 
! !DESCRIPTION:
!  Resets input data filenames for GDAS forcing for execution through GDS
!
! !INTERFACE:
     subroutine reset_gdas_filepaths(gdasdrv)
!EOPI
       use gdasdrv_module
       use opendap_module, only : opendap_data_prefix, ciam
       implicit none
       type(gdasdrvdec)      :: gdasdrv
!BOC
       gdasdrv%gdasdir      = trim(opendap_data_prefix)//'/'// &
                              trim(adjustl(ciam))//'/'//gdasdrv%gdasdir
!EOC
     end subroutine reset_gdas_filepaths
!BOP
! !ROUTINE: init_gdas_vars
! 
! !DESCRIPTION:
! 
! Computes domain decomposition for native as well as 
! interpolated domains for input GDAS forcing
! 
! !INTERFACE:
    subroutine init_gdas_vars()
!EOPI
      use lisdrv_module, only: tile
      implicit none

!BOC
      call set_gdas_lat(gdas_slat,gdas_nlat,gdas_wlon,gdas_elon)
      
      lis%d%ngrid = gdi(iam)
      lis%d%nch   = di_array(iam)
      
      gdas_nc    = ( gdas_elon - gdas_wlon + 1 )
      gdas_nr    = ( gdas_nlat - gdas_slat + 1 )

      fnroffset = gdas_slat - 1
      grid_offset = tile(1)%index-1

      write(cgdas_slat, '(i4)') gdas_slat
      write(cgdas_nlat, '(i4)') gdas_nlat
      write(cgdas_wlon, '(i4)') gdas_wlon
      write(cgdas_elon, '(i4)') gdas_elon

      print*,'DBG: gdas_init -- gdas_slat', gdas_slat, ' (', iam, ')'
      print*,'DBG: gdas_init -- gdas_nlat', gdas_nlat, ' (', iam, ')'
      print*,'DBG: gdas_init -- gdas_wlon', gdas_wlon, ' (', iam, ')'
      print*,'DBG: gdas_init -- gdas_elon', gdas_elon, ' (', iam, ')'
      print*,'DBG: gdas_init -- ngrid', lis%d%ngrid, ' (', iam, ')'
      print*,'DBG: gdas_init -- glbngrid', lis%d%glbngrid, ' (', iam, ')'
      !print*,'DBG: gdas_init -- ncold', gdasdrv%ncold, ' (', iam, ')'
      !print*,'DBG: gdas_init -- nrold', gdasdrv%nrold, ' (', iam, ')'
      print*,'DBG: gdas_init -- lnc', lis%d%lnc, ' (', iam, ')'
      print*,'DBG: gdas_init -- lnr', lis%d%lnr, ' (', iam, ')'
      print*,'DBG: gdas_init -- fnroffset', fnroffset, ' (', iam, ')'
      print*,'DBG: gdas_init -- grid_offset', grid_offset, ' (', iam, ')'
      print*,'DBG: gdas_init -- cgdas_slat ', cgdas_slat, ' (', iam, ')'
      print*,'DBG: gdas_init -- cgdas_nlat ', cgdas_nlat, ' (', iam, ')'
      print*,'DBG: gdas_init -- cgdas_wlon ', cgdas_wlon, ' (', iam, ')'
      print*,'DBG: gdas_init -- cgdas_elon ', cgdas_elon, ' (', iam, ')'
      print*,'DBG: gdas_init -- gdas_nc', gdas_nc, ' (', iam, ')'
      print*,'DBG: gdas_init -- gdas_nr', gdas_nr, ' (', iam, ')'
!EOC
    end subroutine init_gdas_vars
!BOP
! !ROUTINE: def_kgds
! 
! !DESCRIPTION:  
! Initializes the kgds array for GDAS-GDS runs. 
! 
! !INTERFACE:
    subroutine def_kgds(kgdsi)
!EOP      
      integer :: kgdsi(200)
!BOC
      kgdsi(1) = 4
      kgdsi(2) = 512
      kgdsi(3) = 256
      kgdsi(4) = 89463
      kgdsi(5) = 0
      kgdsi(6) = 128
      kgdsi(7) = -89463
      kgdsi(8) = -703
      kgdsi(9) = 703
      kgdsi(10) = 128
      kgdsi(20) = 255
!EOC
    end subroutine def_kgds

!BOP
! !ROUTINE: set_gdas_lat
! 
! !DESCRIPTION:
!
! Computes the latitudes of the decomposed domain for GDAS forcing data
!
! !INTERFACE:
   subroutine set_gdas_lat(slat, nlat, wlon, elon)
      implicit none
! !OUTPUT PARAMETERS:
      integer, intent(out) :: slat, nlat, wlon, elon
!EOP
!BOC
   ! Return global domain
   gdas_slat = 1
   gdas_nlat = 256
   gdas_wlon = 1
   gdas_elon = 512
!EOC
    end subroutine set_gdas_lat

!BOP
! !ROUTINE: twotoone
! 
! !DESCRIPTION:
!
! Remaps a 2-dimensional array of GDAS forcing values into a 
! 1-dimensional array.
!
! !INTERFACE:
subroutine twotoone(f2d,f,nc,nr)
   implicit none

! !INPUT PARAMETERS:
   integer, intent(in) :: nc, nr
   real, dimension(nc,nr), intent(in) :: f2d
! !OUTPUT PARAMETERS:
   real, dimension(nc*nr), intent(out) :: f
!EOP
!BOC
   integer :: i, j, k

   k = 1
   do j = nr, 1, -1
      do i = 1, nc
         f(k) = f2d(i,j)
         k = k + 1
      enddo
   enddo
end subroutine twotoone
!EOC

!BOP
! !ROUTINE: get_kpds
! 
! !DESCRIPTION:
!
! Sets the local kpds array for use in interp\_gdas.
!
! !INTERFACE:
subroutine get_kpds(kpds,iv)
   implicit none
! !INPUT PARAMETERS:
   integer, intent(in) :: iv
! !OUTPUT PARAMETERS:
   integer, intent(out) :: kpds(200)
!EOP

!BOC
   kpds = 0

   select case ( iv )
   case ( 1 )
      kpds(1:25) = (/7,82,255,128,11,105,2,1,6,10,18,0,1,3,0,&
                     10,0,1,2,0,21,1,0,0,32/)
   case ( 2 )
      kpds(1:25) = (/7,82,255,128,51,105,2,1,6,10,18,0,1,3,0,&
                     10,0,1,2,0,21,4,0,0,32/)
   case ( 3 )
      kpds(1:25) = (/7,82,255,128,204,1,0,1,6,10,18,0,1,0,3,3,&
                     0,1,2,0,21,0,0,0,32/)
   case ( 4 )
      kpds(1:25) = (/7,82,255,128,205,1,0,1,6,10,18,0,1,0,3,3,0,&
                     1,2,0,21,0,0,0,32/)
   case ( 5 )
      kpds(1:25) = (/7,82,255,128,33,105,10,1,6,10,18,0,1,3,0,10,&
                     0,1,2,0,21,1,0,0,32/)
   case ( 6 )
      kpds(1:25) = (/7,82,255,128,34,105,10,1,6,10,18,0,1,3,0,10,&
                     0,1,2,0,21,1,0,0,32/)
   case ( 7 )
      kpds(1:25) = (/7,82,255,128,1,1,0,1,6,10,18,0,1,3,0,10,0,1,&
                     2,0,21,-1,0,0,32/)
   case ( 8 )
      kpds(1:25) = (/7,82,255,128,59,1,0,1,6,10,18,0,1,0,3,3,0,1,&
                     2,0,21,6,0,0,32/)
   case ( 9 )
      kpds(1:25) = (/7,82,255,128,214,1,0,1,6,10,18,0,1,0,3,3,0,&
                     1,2,0,21,6,0,0,32/)
   case ( 10 )
      kpds(1:25) = (/7,82,255,128,84,1,0,1,6,10,18,0,1,0,3,3,0,1,&
               2,0,21,1,0,0,32/)
   case default 
      print*,'ERR: get_kpds -- forcing index ',iv,'is out of range (10)',&
             ' (',iam,')'
      stop 344
   end select
end subroutine get_kpds
!EOC

!BOP
! !ROUTINE: get_kgds
! 
! !DESCRIPTION:
!
! Sets the local kgds array for use in interp\_gdas.
!
! !INTERFACE:
subroutine get_kgds(kgds)
   implicit none
! !OUTPUT PARAMETERS:
   integer, intent(out) :: kgds(200)
!EOP

!BOC
   kgds = 0
   kgds(1:20) = (/4,512,256,89463,0,128,-89463,-703,703,128,&
                  0,0,0,0,0,0,0,0,0,255/)

end subroutine get_kgds
!EOC

!BOP
! !ROUTINE: get_gridDesc
! 
! !DESCRIPTION:
!
! Sets the local grid description array for use in interp\_gdas.
!
! !INTERFACE:
subroutine get_gridDesc(gridDesc)
   implicit none
! !OUTPUT PARAMETERS:
   integer, intent(out) :: gridDesc(50)
!EOP

!BOC
   gridDesc = 0
   gridDesc(1:20) = (/4,512,256,89463,0,128,-89463,-703,703,128,&
                      0,0,0,0,0,0,0,0,0,255/)
   call lis_log_msg('DBG: get_gridDesc -- FIX THIS - DOES NOT HANDLE GRID CHANGE')

end subroutine get_gridDesc

!BOP
! !ROUTINE: set_lb
! 
! !DESCRIPTION:
!
! Sets the local lb mask array for use in interp\_gdas.
!
! !INTERFACE:
subroutine set_lb(f,lb,nc,nr)

   implicit none

! !INPUT PARAMETERS:
   integer, intent(in) :: nc,nr
   real, dimension(nc*nr), intent(in) :: f
! !OUTPUT PARAMETERS:
   logical*1, intent(out) :: lb(nc*nr)
!EOP

!BOC
   integer :: i

   lb = .false.

   do i = 1, nc*nr
      if ( f(i) > 0.0 ) then
         lb(i) = .true.
      endif
   enddo

end subroutine set_lb
!EOC

#endif

  end module gdasopendap_module
