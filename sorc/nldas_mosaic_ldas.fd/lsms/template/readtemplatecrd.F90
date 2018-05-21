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
! !ROUTINE: readtemplatecrd.F90
!
! !DESCRIPTION:
!  Routine to read Template specific parameters from the card file. 
!
! !REVISION HISTORY:
! 21 Jul 2004; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readtemplatecrd(templatedrv)
! !USES:
  use templatedrv_module
!EOP
  implicit none
  integer :: lsm
  type(templatedrvdec) :: templatedrv
  namelist /template/templatedrv
!BOC
  open(11,file='lis.crd',form='formatted',status='old')
  read(unit=11,NML=template)
  print*,'Running TEMPLATE LSM:'

  close(11)
!EOC
end subroutine readtemplatecrd
