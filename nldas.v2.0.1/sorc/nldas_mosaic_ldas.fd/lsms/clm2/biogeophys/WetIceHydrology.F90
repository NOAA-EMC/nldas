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

subroutine WetIceHydrology (clm)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Calculate hydrology for ice and wetland
!
! Method:
!
! Author:
! 7 November 2000: Mariana Vertenstein; Initial code
!
!-----------------------------------------------------------------------
! $Id: WetIceHydrology.F90,v 1.6 2004/11/24 22:56:50 jim Exp $
!-----------------------------------------------------------------------

  use precision
  use clmtype
  use clm_varcon, only : istwet, istice
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

!----Local Variables----------------------------------------------------
! None
!----End Variable List--------------------------------------------------

!
! Wetland and land ice runoff
!

  clm%qflx_drain  = 0. 
  clm%qflx_surf   = 0.
  clm%qflx_infl   = 0.
  clm%qflx_qrgwl  = clm%forc_rain + clm%forc_snow - clm%qflx_evap_tot - &
                    (clm%endwb - clm%begwb)/clm%dtime

end subroutine WetIceHydrology
