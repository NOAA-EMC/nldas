#include <misc.h>

subroutine clm1_condch (ra,            rb,          rd,          &
                       wta,           wtl,         wtg,         &    
                       wta0,          wtl0,        wtg0,        &  
                       wtal,          wtga,        wtgl,      clm1)

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
!  Provides dimensional and non-dimensional sensible heat
!  conductances for canopy and soil flux calculations.
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!=========================================================================
! $Id: clm1_condch.F90,v 1.1.1.1 2003/02/06 16:10:43 jgottsch Exp $
!=========================================================================

! Declare Modules and data structures

  use precision
  use clm1type
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout)  :: clm1	 !CLM 1-D Module

  real(r8), intent(in) ::  &
       ra,              & ! aerodynamical resistance [s/m]
       rb,              & ! leaf boundary layer resistance [s/m]
       rd                 ! thermal resistance between ground and bottom of canopy

  real(r8), intent(out) :: &
       wta,             & ! heat conduactance for air [m/s]
       wtg,             & ! heat conduactance for ground [m/s]
       wtl,             & ! heat conduactance for leaf [m/s]
       wta0,            & ! normalized heat conductance for air [-]
       wtl0,            & ! normalized heat conductance for air [-]
       wtg0,            & ! normalized heat conductance for ground [-]
       wtal,            & ! normalized heat conductance for air and leaf [-]
       wtgl,            & ! normalized heat conductance for leaf and ground [-]
       wtga               ! normalized heat conductance for air and ground  [-]

!=== Local Variables =====================================================

  real(r8)  wtshi            ! heat resistance for air, ground and leaf [s/m]

!=== End Variable List ===================================================

  wta   = clm1%frac_veg_nosno/ra                     ! air
  wtl   = clm1%frac_veg_nosno*(clm1%elai+clm1%esai)/rb ! leaf
  wtg   = clm1%frac_veg_nosno/rd                     ! ground
  wtshi = 1./(wta+wtl+wtg)

  wtl0  = wtl*wtshi         ! leaf
  wtg0  = wtg*wtshi         ! ground
  wta0  = wta*wtshi         ! air

  wtgl  = wtl0+wtg0         ! ground + leaf
  wtga  = wta0+wtg0         ! ground + air
  wtal  = wta0+wtl0         ! air + leaf

end subroutine clm1_condch 
