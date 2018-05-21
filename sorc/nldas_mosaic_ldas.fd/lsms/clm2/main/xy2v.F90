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

subroutine xy2v (nx, ny, fldxy, ki, kf, fldv)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! convert a grid-average field to subgrid patch vector
! 
! Method: 
! This code converts a grid-average field [fldxy] dimensioned
! [lsmlon] x [lsmlat] to a subgrid patch vector [fldv] for 
! [numpatch] subgrid patches. 
!
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: xy2v.F90,v 1.6 2004/11/24 22:57:24 jim Exp $
!-----------------------------------------------------------------------

  use precision
  use clm_varmap, only : patchvec
  implicit none

! ------------------------ arguments----------------------------------
  integer , intent(in)  :: nx, ny          !x-y dimension	
  integer , intent(in)  :: ki, kf          !beginning and end patch indices
  real(r8), intent(in)  :: fldxy(nx,ny)    !gridded input
  real(r8), intent(out) :: fldv(ki:kf)     !subgrid vector output
! --------------------------------------------------------------------

! ------------------------ local variables ----------------------
  integer i,j,k             !indices
! ---------------------------------------------------------------

  do k = ki,kf
     i = patchvec%ixy(k)
     j = patchvec%jxy(k)
     fldv(k) = fldxy(i,j)
  end do

  return
end subroutine xy2v





