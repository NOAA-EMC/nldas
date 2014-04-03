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
! !ROUTINE : readclm2crd.F90
!
! !DESCRIPTION:
!  Routine to read CLM specific parameters from the card file. 
!
! !REVISION HISTORY:
! 14 Oct 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readclm2crd(clmdrv)
! !USES:
  use clm2drv_module
!EOP
  implicit none
  type(clmdrvdec) :: clmdrv
  namelist /clm2/clmdrv
!BOC
  open(11,file='lis.crd',form='formatted',status='old')
  read(unit=11,NML=clm2)
  print*,'Running CLM2 LSM:'
  print*,'CLM2 Active Restart File: ', clmdrv%CLM2_RFILE
  clmdrv%clm2open=0
  close(11)
!EOC
end subroutine readclm2crd
