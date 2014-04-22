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

subroutine v2xy (fldv, fldxyini, fldxy)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Perfrom grid-average from subgrid patch vector
! 
! Method: 
! Subgrid patch to grid average mapping: average a subgrid input vector 
! [fldv] of length [numpatch] to a 2-d [lsmlon] x [lsmlat] output array [fldxy]
! setting non-land points to [nonland]. Averaging is only done for points 
! that are not equal to "spval".
!
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: v2xy.F90,v 1.6 2004/11/24 22:57:23 jim Exp $
!-----------------------------------------------------------------------

  use precision
  use clm_varder 
!  use clm_varsur, only : landmask
!  use clm_varpar, only : lsmlon, lsmlat
  use clm_varmap, only : numpatch, patchvec
  use clm_varcon, only : spval
  implicit none

! ------------------------ arguments----------------------------------
  real(r8), intent(in)  :: fldv(numpatch)       !subgrid vector input
  real(r8), intent(in)  :: fldxyini             !initial value of fldxy
! quick fix
  real(r8), intent(out) :: fldxy(2,2) !gridded output
!  real(r8), intent(out) :: fldxy(lsmlon,lsmlat) !gridded output
! --------------------------------------------------------------------

! ------------------------ local variables ----------------------
!  integer :: i                    !longitude index     
!  integer :: j                    !latitude index      
!  integer :: k                    !subgrid patch index 
!  real(r8):: sumwt(lsmlon,lsmlat) !sum of wt
! ---------------------------------------------------------------

! Loop over subgrid patches to create grid average. 
  print*, 'not supposed to be in v2xy..'
#if 0 
  fldxy(:,:) = fldxyini
  sumwt(:,:) = 0.
  do k = 1, numpatch
     if (fldv(k) /= spval) then
        i = patchvec%ixy(k)      !longitude index for land point
        j = patchvec%jxy(k)      !latitude index for land point
        if (sumwt(i,j)==0.) fldxy(i,j) = 0.
        fldxy(i,j) = fldxy(i,j) + patchvec%wtxy(k)*fldv(k)
        sumwt(i,j) = sumwt(i,j) + patchvec%wtxy(k)
     endif
  end do
  where (landmask(:,:) == 1 .and. sumwt(:,:) /= 0.)
     fldxy(:,:) = fldxy(:,:)/sumwt(:,:)
  endwhere
#endif
  return
end subroutine v2xy
