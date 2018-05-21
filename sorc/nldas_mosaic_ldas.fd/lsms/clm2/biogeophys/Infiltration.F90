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

subroutine Infiltration (clm) 

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
! Calculate infiltration into surface soil layer (minus the evaporation)
!
! Method:
! The original code on soil moisture and runoff were provided by 
! R. E. Dickinson in July 1996.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: Infiltration.F90,v 1.6 2004/11/24 22:56:29 jim Exp $
!-----------------------------------------------------------------------

  use precision
  use clmtype
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

!----Local Variables----------------------------------------------------
!None
!----End Variable List--------------------------------------------------

!
! Infiltration into surface soil layer (minus the evaporation)
!

  if (clm%snl+1 >= 1) then
     clm%qflx_infl = clm%qflx_top_soil - clm%qflx_surf - clm%qflx_evap_grnd
  else
     clm%qflx_infl = clm%qflx_top_soil - clm%qflx_surf
  endif
  
end subroutine Infiltration
