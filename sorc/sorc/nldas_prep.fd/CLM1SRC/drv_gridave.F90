#include <misc.h>

function drv_gridave (nch, mask, fgr, var)
!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
! The following function simply averages tile spatial arrays.
!
! REVISION HISTORY:
!  6 May 1999: Paul Houser; initial code
!  7 Dec 2000: Mariana Vertenstein
!=========================================================================

  use precision
  implicit none

!=== Arguments =============================================================

  integer , intent(in) :: nch        ! number of tiles
  integer , intent(in) :: mask(nch)  ! mask (1=not used, 0=used)
  real(r8), intent(in) :: fgr(nch)   ! Fraction of veg class in grid
  real(r8), intent(in) :: var(nch)   ! CLM Variable to average

!=== Local Variables =====================================================
  
  integer  :: t                ! tile index
  real(r8) :: drv_gridave      ! Spatial average function  

!=== End Variable List ===================================================

  drv_gridave = 0.0
  do t = 1,nch
     if (mask(t) == 1) then
        drv_gridave = drv_gridave + var(t)*fgr(t)
     endif
  enddo
  drv_gridave=drv_gridave/float(nch)

end function drv_gridave













