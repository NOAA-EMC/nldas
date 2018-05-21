#include <misc.h>

subroutine clm1_snowage (clm1)

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
!  Updates snow cover and snow age, based on BATS code.
!
! REVISION HISTORY:
!  Original Code:  Robert Dickinson
!  15 September 1999: Yongjiu Dai; Integration of code into CLM
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!=========================================================================
! $Id: clm1_snowage.F90,v 1.1.1.1 2003/02/06 16:10:45 jgottsch Exp $
!=========================================================================

! Declare Modules and data structures

  use precision
  use clm1type
  use clm1_varcon, only : tfrz
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1  !CLM 1-D Module

!=== Local Variables =====================================================

  real(r8)                 & !
       age1,               & ! snow aging factor due to crystal growth [-]
       age2,               & ! snow aging factor due to surface growth [-]
       age3,               & ! snow aging factor due to accum of other particles [-]
       arg,                & ! temporary variable used in snow age calculation [-]
       arg2,               & ! temporary variable used in snow age calculation [-]
       dela,               & ! temporary variable used in snow age calculation [-]
       dels,               & ! temporary variable used in snow age calculation [-]
       sge                   ! temporary variable used in snow age calculation [-]

!=== End Variable List ===================================================

  if (clm1%h2osno <= 0.) then

     clm1%snowage = 0.

  else if (clm1%h2osno > 800.) then   ! Over Antarctica

     clm1%snowage = 0.

  else                               ! Away from Antarctica 

     age3  = 0.3
     arg   = 5.e3*(1./tfrz-1./clm1%t_grnd)
     arg2  = min(0.,10.*arg)
     age2  = exp(arg2)
     age1  = exp(arg)
     dela  = 1.e-6*clm1%dtime*(age1+age2+age3)
     dels  = 0.1*max(0.0, clm1%h2osno-clm1%h2osno_old)
     sge   = (clm1%snowage+dela)*(1.0-dels)
     clm1%snowage   = max(0.0,sge)

  endif

end subroutine clm1_snowage
