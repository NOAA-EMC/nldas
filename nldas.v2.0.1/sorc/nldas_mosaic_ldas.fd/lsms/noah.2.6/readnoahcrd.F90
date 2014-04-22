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
! !ROUTINE: readnoahcrd.F90
!
! !DESCRIPTION:
!  Routine to read Noah specific parameters from the card file. 
!
! !REVISION HISTORY:
! 14 Oct 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readnoahcrd(noahdrv)
! !USES:
  use noahdrv_module
!EOP
  implicit none
  integer :: lsm
  type(noahdrvdec) :: noahdrv
  namelist /noah/noahdrv
!BOC
  open(11,file='lis.crd',form='formatted',status='old')
  read(unit=11,NML=noah)
  print*,'Running NOAH LSM:'
  print*,'NOAH Active Restart File: ', noahdrv%NOAH_RFILE
  noahdrv%noah_gfractime = 0.0
  noahdrv%noah_albtime = 0
  noahdrv%noah_albdchk = 0
  noahdrv%noah_gfracdchk = 0
  noahdrv%NOAHOPEN=0
  noahdrv%numout = 0
  noahdrv%NOAH_ZST     = 9

  close(11)
!EOC
end subroutine readnoahcrd
