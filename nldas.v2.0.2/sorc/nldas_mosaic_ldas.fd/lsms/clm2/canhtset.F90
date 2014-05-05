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
!BOP
!
! !ROUTINE:canhtset.F90
!
! !DESCRIPTION:
!  This subroutine reads in canopy height information into CLM2 for now
!
! !REVISION HISTORY:
!  15 Nov 2002: Jon Gottschalck; Initial code
!
! !INTERFACE:
subroutine canhtset ()
! !USES:
  use clm_varder          ! CLM2 tile variables
  use clm_varpar     , only : maxpatch
  use clm_varmap     , only : numpatch
  use clm_varctl     , only : clmdrv
  use lis_module
!EOP
  implicit none

!=== Arguments ===========================================================

!=== Local variables

  integer :: i,t
  real    :: pthtop(maxpatch),pthbot(maxpatch)  ! Canopy top and bottom height for the 13 UMD vegetation types

  print*, 'DBG: canhtset -- numpatch,maxpatch',numpatch,maxpatch,' (',iam,')'
!BOC
!---------------------------------------------------------------------
! Open canopy heights file and read into temporary variables
!---------------------------------------------------------------------
  open(57,file=clmdrv%clm2_chtfile,status='old')
  read(57,*) (pthtop(i),i=1,maxpatch)
  read(57,*) (pthbot(i),i=1,maxpatch)
  close(57)
!---------------------------------------------------------------------
! Assign canopy top and height values for each 
! vegetation type to CLM2 tiles
!---------------------------------------------------------------------
  do t=1,numpatch
   clm(t)%htop = pthtop(clm(t)%itypveg)
   clm(t)%hbot = pthbot(clm(t)%itypveg)
  enddo
  return
!EOC
end subroutine canhtset


