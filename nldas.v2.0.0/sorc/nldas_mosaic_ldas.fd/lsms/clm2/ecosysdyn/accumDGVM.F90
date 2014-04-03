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

module accumDGVM

#if (defined DGVM)

  use precision
  implicit none

! public methods

  public :: accumDGVMini

!=======================================================================
CONTAINS
!=======================================================================

  subroutine accumDGVMini
    
    use precision
    use clm_varder
    use accumulMod, only : accext
    use clm_varmap, only : begpatch,endpatch
    use time_manager, only: get_nstep
    implicit none
    
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize time variant DGVM variables
! 
! Method: 
! This routine must be called after the restart file is read
! because the time manager is initialized in the restart file
! and that is needed to obtain the time step
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: accumDGVM.F90,v 1.6 2004/11/24 22:57:03 jim Exp $
!-----------------------------------------------------------------------

! ------------------------ local variables ---------------------------
    integer nstep ! time step
!-----------------------------------------------------------------------

! Initialize variables to be time accumulated for various purposes 

    nstep = get_nstep()
    call accext ('T10     ', clm(begpatch:endpatch)%t10    , nstep)
    call accext ('TDA     ', clm(begpatch:endpatch)%t_mo   , nstep)
    call accext ('AGDD0   ', clm(begpatch:endpatch)%agdd0  , nstep)
    call accext ('AGDD5   ', clm(begpatch:endpatch)%agdd5  , nstep)
    call accext ('FNPSN10 ', clm(begpatch:endpatch)%fnpsn10, nstep)
    call accext ('PREC365 ', clm(begpatch:endpatch)%prec365, nstep)
    call accext ('AGDDTW  ', clm(begpatch:endpatch)%agddtw , nstep)
    call accext ('AGDD    ', clm(begpatch:endpatch)%agdd   , nstep)
    
    return
  end subroutine accumDGVMini

#endif

end module accumDGVM
