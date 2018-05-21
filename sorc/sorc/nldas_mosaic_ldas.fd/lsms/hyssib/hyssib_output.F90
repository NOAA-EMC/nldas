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
! !ROUTINE: hyssib_output.F90
!
! !DESCRIPTION: This subroutines sets up methods to write HY-SSiB output 
!
! !REVISION HISTORY:
!    Dec 2003: Luis-Gustavo Goncalves, Initial version
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
subroutine hyssib_output
! !USES:
  use lisdrv_module, only : lis, tile, glbgindex
  use hyssib_varder, only : hyssibdrv
  use spmdMod, only : masterproc
!EOP
  integer :: i 
  real :: var(lis%d%glbnch)
  
  if (mod(lis%t%gmt,hyssibdrv%writeint).eq.0) then
     call hyssib_gather()
     if (masterproc) then 
        call hyssib_almaout()
     endif
  endif
  
end subroutine hyssib_output

