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
!BOP
!
! !ROUTINE: noah_setvegparms.F90
!
! !DESCRIPTION:
!  This subroutine retrieves NOAH parameters - Significant F90 revisions
!   below this subroutine will be required in the future.  
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added NOAH LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!
! !INTERFACE:
subroutine noah_setvegparms
! !USES:
  use lisdrv_module, only : tile,lis
  use noah_varder      ! NOAH tile variables
!EOP      
  implicit none

  integer :: n,i,j
  integer :: ios1
  real :: value(lis%p%nt,noahdrv%noah_nvegp)
!BOC  
  
  do n=1,lis%d%nch
     call mapvegc(tile(n)%vegt)
     noah(n)%vegt = tile(n)%vegt
  enddo   
!-----------------------------------------------------------------------
! Get Vegetation Parameters for NOAH Model in Tile Space
! Read in the NOAH Static Vegetation Parameter Files
!-----------------------------------------------------------------------
  open(unit=11,file=noahdrv%noah_vfile,status='old')
  
  do j=1,noahdrv%noah_nvegp
     read(11,*)(value(i,j),i=1,lis%p%nt)
  enddo
  close(11)
!-----------------------------------------------------------------------
! Assign STATIC vegetation parameters to each tile based on the
! type of vegetation present in that tile.
! These parameters will be stored in one long array--structured
! as follows: Tile 1, all the parameters (1 through numparam)
! then Tile 2, all the parameters. 
! Then Tile 3, all the parameters etc.
!-----------------------------------------------------------------------

  do i=1,lis%d%nch
     do j=1,noahdrv%noah_nvegp
        noah(i)%vegp(j)=value(tile(i)%vegt,j)
     enddo 
  enddo 
end subroutine noah_setvegparms
