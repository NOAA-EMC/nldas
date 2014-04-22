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
! !ROUTINE: vic_tilemapping.F90
!
! !DESCRIPTION:
! Returns the column and row indices of a tile given its grid number index
!
! !REVISION HISTORY:
! 
! 27 Feb 2004; James Geiger : Initial Specification
! 
! !INTERFACE:
subroutine vic_tilemapping(k,c,r)
!USES:
   use lisdrv_module, only: tile

   implicit none

! !INPUT PARAMETERS:
   integer, intent(in) :: k     ! grid number index

! !OUTPUT PARAMETERS:
   integer, intent(out) :: c,r  ! column number, row number

!EOP

!BOC

   c = tile(k)%col
   r = tile(k)%row

!EOC

end subroutine vic_tilemapping
