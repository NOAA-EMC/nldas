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

function StabilityFunc(k, zeta)

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
! Stability function for rib < 0.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: StabilityFunc.F90,v 1.6 2004/11/24 22:56:43 jim Exp $
!-----------------------------------------------------------------------

  use precision
  use shr_const_mod, only: SHR_CONST_PI
  implicit none

!----Local Variables----------------------------------------------------

  integer k         !
  real(r8) zeta     ! dimensionless height used in Monin-Obukhov theory
  real(r8) StabilityFunc  ! stability function for unstable case
  real(r8) chik     ! 

!=== End Variable List ===================================================

  chik = (1.-16.*zeta)**0.25
  if (k == 1) then
     StabilityFunc = 2.*log((1.+chik)*0.5) &
          + log((1.+chik*chik)*0.5)-2.*atan(chik)+SHR_CONST_PI*0.5
  else
     StabilityFunc = 2.*log((1.+chik*chik)*0.5)
  endif

end function StabilityFunc
