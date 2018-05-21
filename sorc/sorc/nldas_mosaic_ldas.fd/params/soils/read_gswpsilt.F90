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
#include "misc.h"
!BOP
!
! !ROUTINE: read_gswpsilt
!
! !DESCRIPTION:
!  This subroutine retrieves GSWP-2 silt data

! !REVISION HISTORY:
!  15 Jun 2005: James Geiger; Initial Specification
!
! !INTERFACE:
subroutine read_gswpsilt(array)
! !USES:
#if ( defined USE_NETCDF )
  use netcdf
#endif
  use lisdrv_module, only : lis, gindex
  use lis_indices_module
!EOP      

  implicit none

  real, intent(inout)             :: array(lis_nc_data, lis_nc_data)
  real, allocatable, dimension(:) :: silt
  integer                         :: ncid, status
  integer                         :: siltid
  integer                         :: cindex,rindex,c,r

#if ( defined USE_NETCDF )
  call lis_log_msg('MSG: read_gswpsilt -- Reading GSWP silt file')

  allocate(silt(lis%d%glbnch))
  
  call lis_log_msg('MSG: read_gswpsilt -- Reading silt file: '// &
                   trim(lis%p%sifile))

  status = nf90_open(path=lis%p%sifile, mode=nf90_nowrite, ncid=ncid)
  status = nf90_inq_varid(ncid, "Silt", siltid)
  status = nf90_get_var(ncid, siltid, silt)
  status = nf90_close(ncid)

!GSWP vector to 2D:
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if ( gindex(cindex,rindex) /= -1 ) then 
           array(cindex,rindex) = silt(gindex(cindex,rindex))
        endif
     enddo
  enddo
!EOC
#endif

end subroutine read_gswpsilt

