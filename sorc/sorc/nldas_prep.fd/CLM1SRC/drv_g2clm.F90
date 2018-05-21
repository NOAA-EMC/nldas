#include <misc.h>

subroutine drv_g2clm(u,drv,grid,tile,clm1)

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
!  Transfer variables into CLM tile space from grid space. 
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!=========================================================================
! $Id: drv_g2clm.F90,v 1.1.1.1 2003/02/06 16:10:44 jgottsch Exp $
!=========================================================================

  use precision
  use drv_module          ! 1-D Land Model Driver variables
  use drv_tilemodule      ! Tile-space variables
  use drv_gridmodule      ! Grid-space variables
  use clm1type             ! CLM tile variables
  use clm1_varcon, only : istdlak, istslak
  implicit none

!=== Arguments ===========================================================

  type (drvdec)      :: drv              
  type (clm_tiledec) :: tile(drv%nch)
  type (clm_griddec) :: grid(drv%nc,drv%nr)   
  type (clm11d)       :: clm1(drv%nch)     

!=== Local Variables =====================================================

  integer  :: j,c,r,t     ! Loop counters
  real(r8) :: u         ! Tempoary UNDEF Variable  

!=== End Variable Definition =============================================

!=== Transfer variables that are identical for each tile in a grid 
!===  to tile space  - for some aplications, the assumption
!===  of spatially constant information across a grid for these
!===  variables may be incorrect, and should be modified by the user.
  do t=1,drv%nch
  c=tile(t)%col
  r=tile(t)%row

! CLM Forcing parameters (read into 2-D grid module variables)

  if (grid(c,r)%forc_hgt_u      /= u) clm1(t)%forc_hgt_u =grid(c,r)%forc_hgt_u  
  if (grid(c,r)%forc_hgt_t      /= u) clm1(t)%forc_hgt_t =grid(c,r)%forc_hgt_t  
  if (grid(c,r)%forc_hgt_q      /= u) clm1(t)%forc_hgt_q =grid(c,r)%forc_hgt_q  

! CLM Vegetation parameters (read into 2-D grid module variables)

  if (grid(c,r)%dewmx           /= u) clm1(t)%dewmx      =grid(c,r)%dewmx    

! CLM Soil parameters	(read into 2-D grid module variables)

  if (grid(c,r)%smpmax          /= u) clm1(t)%smpmax     =grid(c,r)%smpmax    
  if (grid(c,r)%scalez          /= u) tile(t)%scalez    =grid(c,r)%scalez   
  if (grid(c,r)%hkdepth         /= u) tile(t)%hkdepth   =grid(c,r)%hkdepth    
  if (grid(c,r)%wtfact          /= u) clm1(t)%wtfact     =grid(c,r)%wtfact    

! Roughness lengths (read into 2-D grid module variables)

  if (grid(c,r)%zlnd            /= u) clm1(t)%zlnd       =grid(c,r)%zlnd    
  if (grid(c,r)%zsno            /= u) clm1(t)%zsno       =grid(c,r)%zsno    
  if (grid(c,r)%csoilc          /= u) clm1(t)%csoilc     =grid(c,r)%csoilc    

! Numerical finite-difference parameters (read into 2-D grid module variables)

  if (grid(c,r)%capr            /= u) clm1(t)%capr       =grid(c,r)%capr    
  if (grid(c,r)%cnfac           /= u) clm1(t)%cnfac      =grid(c,r)%cnfac    
  if (grid(c,r)%smpmin          /= u) clm1(t)%smpmin     =grid(c,r)%smpmin    
  if (grid(c,r)%ssi             /= u) clm1(t)%ssi        =grid(c,r)%ssi    
  if (grid(c,r)%wimp            /= u) clm1(t)%wimp       =grid(c,r)%wimp    
  if (grid(c,r)%pondmx          /= u) clm1(t)%pondmx     =grid(c,r)%pondmx    
  enddo
  return
end subroutine drv_g2clm



