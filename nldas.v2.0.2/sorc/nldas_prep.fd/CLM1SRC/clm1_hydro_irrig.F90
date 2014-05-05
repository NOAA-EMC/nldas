#include <misc.h>

subroutine clm1_hydro_irrig (clm1)

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
! Irrigate crops to depth of 1 m and modify runoff accordingly
!
! REVISION HISTORY:
!  7 November 2000: Mariana Vertenstein; Initial code
!
!=========================================================================
! $Id: clm1_hydro_irrig.F90,v 1.1.1.1 2003/02/06 16:10:45 jgottsch Exp $
!=========================================================================

  use precision
  use clm1type
  use clm1_varpar, only : nlevsoi
  use clm1_varcon, only : denh2o
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1	 !CLM 1-D Module

!=== Local Variables =====================================================

  integer i        !index
  real(r8) wtold   !temporary var
  real(r8) wtnew   !temporary var

!=========================================================================

! To saturate the soil, make the liquid water volume equal the
! effective porosity. Irrigate soil to a depth of 30 cm.

  clm1%qflx_qirr = 0.
  if (clm1%irrig .and. (clm1%elai+clm1%esai)>0) then
     wtold = 0.
     wtnew = 0.
     do i = 1, nlevsoi
        if (clm1%zi(i) <= 0.3) then
           wtold = wtold + clm1%h2osoi_liq(i)
           clm1%h2osoi_liq(i) = clm1%eff_porosity(i)*clm1%dz(i)*denh2o
           wtnew = wtnew + clm1%h2osoi_liq(i)
        end if
     end do
     clm1%qflx_qirr = (wtnew - wtold) / clm1%dtime
  end if

end subroutine clm1_hydro_irrig





