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
! !ROUTINE: readperscrd.F90
!
! !DESCRIPTION:
!  Routine to read PERS specific parameters from the card file. 
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readperscrd(persdrv)
! !USES:
  use persdrv_module
!EOP
  implicit none
  integer :: lsm
  type(persdrvdec) :: persdrv
  namelist /pers/persdrv
!BOC
  open(11,file='lis.crd',form='formatted',status='old')
  read(unit=11,NML=pers)
  print*,'Using PERS forcing'
  print*, 'PERS forcing directory :',persdrv%PERSDIR
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
  persdrv%perstime = 0.0
  close(11)
!EOC
end subroutine readperscrd
