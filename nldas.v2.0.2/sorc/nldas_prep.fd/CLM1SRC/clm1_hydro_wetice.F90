#include <misc.h>

subroutine clm1_hydro_wetice (clm1)

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
!  Calculate hydrology for ice and wetland
!
! REVISION HISTORY:
!  7 November 2000: Mariana Vertenstein; Initial code
!
!=========================================================================
! $Id: clm1_hydro_wetice.F90,v 1.1.1.1 2003/02/06 16:10:47 jgottsch Exp $
!=========================================================================

  use precision
  use clm1type
  use clm1_varcon, only : istwet, istice
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1	 !CLM 1-D Module

!=========================================================================

! Wetland and land ice runoff

  clm1%qflx_drain  = 0. 
  clm1%qflx_surf   = 0.
  clm1%qflx_infl   = 0.
  clm1%qflx_qirr   = 0.
  clm1%qflx_qrgwl  = clm1%forc_rain + clm1%forc_snow - clm1%qflx_evap_tot - (clm1%endwb - clm1%begwb)/clm1%dtime

end subroutine clm1_hydro_wetice
