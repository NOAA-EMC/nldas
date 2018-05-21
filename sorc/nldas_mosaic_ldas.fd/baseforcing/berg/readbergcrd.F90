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
! !ROUTINE: readbergcrd.F90
!
! !DESCRIPTION:
!  Routine to read BERG specific parameters from the card file. 
!
! !REVISION HISTORY:
! 26Jan2004; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readbergcrd(bergdrv)
! !USES:
  use bergdrv_module
  use lisdrv_module, only : lis
!EOP
  implicit none
  integer :: lsm
  type(bergdrvdec) :: bergdrv
  namelist /berg/bergdrv
!BOC
  open(11,file='lis.crd',form='formatted',status='old')
  read(unit=11,NML=berg)
  print*,'Using BERG forcing'
  bergdrv%fmodeltime1 = 3000.0
  bergdrv%fmodeltime2 = 0.0
  bergdrv%nmif = lis%f%nmif
  close(11)
!EOC
end subroutine readbergcrd
