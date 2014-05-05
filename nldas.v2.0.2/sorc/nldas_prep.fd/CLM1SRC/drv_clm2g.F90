#include <misc.h>

subroutine drv_clm2g (drv,grid,tile,clm1)

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
!  Transfer variables from CLM tile space to grid space. 
!  A few examples are given here - the user may need to add more
!  Note that tile to grid space transfers for written output can and
!  are simplified and independent from this subroutine.
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!=========================================================================
! $Id: drv_clm2g.F90,v 1.1.1.1 2003/02/06 16:10:46 jgottsch Exp $
!=========================================================================

  use precision
  use drv_module          ! 1-D Land Model Driver variables
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! 1-D clm1 variables
  use drv_gridmodule      ! Grid-space variables
  implicit none

!=== Arguments ===========================================================

  type (drvdec)      :: drv              
  type (clm_tiledec) :: tile(drv%nch)
  type (clm11d)       :: clm1 (drv%nch)
  type (clm_griddec) :: grid(drv%nc,drv%nr)   

!=== End Variable Definition =============================================

  call drv_t2gr(clm1%qflx_evap_tot ,grid%qflx_evap_tot ,drv%nc,drv%nr,drv%nch,tile%fgrd,tile%col,tile%row)
  call drv_t2gr(clm1%eflx_lh_tot   ,grid%eflx_lh_tot   ,drv%nc,drv%nr,drv%nch,tile%fgrd,tile%col,tile%row)
  call drv_t2gr(clm1%eflx_lwrad_out,grid%eflx_lwrad_out,drv%nc,drv%nr,drv%nch,tile%fgrd,tile%col,tile%row)
  call drv_t2gr(clm1%t_ref2m       ,grid%t_ref2m       ,drv%nc,drv%nr,drv%nch,tile%fgrd,tile%col,tile%row)
  call drv_t2gr(clm1%t_rad         ,grid%t_rad         ,drv%nc,drv%nr,drv%nch,tile%fgrd,tile%col,tile%row)

  return
end subroutine drv_clm2g





