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
! !ROUTINE: readnldascrd.F90
!
! !DESCRIPTION:
!  Routine to read NLDAS2 specific parameters from the card file. 
!
! !REVISION HISTORY:
! 02Feb2004; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readnldas2crd(nldas2drv,gridDesci)
! !USES:
  use nldas2drv_module
!EOP
  implicit none
  integer :: lsm
  type(nldas2drvdec) :: nldas2drv
  namelist /nldas2/nldas2drv
  real, intent(inout) :: gridDesci(50)
!BOC
  open(11,file='lis.crd',form='formatted',status='old')
  read(unit=11,NML=nldas2)
  print*,'Using NLDAS2 forcing'
  print*, 'NLDAS2 forcing directory :',nldas2drv%NLDAS2DIR
  nldas2drv%nldas2time1 = 3000.0
  nldas2drv%nldas2time2 = 0.0

  gridDesci(1) = 0
  gridDesci(2) = nldas2drv%ncold
  gridDesci(3) = nldas2drv%nrold
  gridDesci(4) = 25.063
  gridDesci(5) = -124.938
  gridDesci(6) = 128
  gridDesci(7) = 52.938
  gridDesci(8) = -67.063
  gridDesci(9) = 0.125
  gridDesci(10) = 0.125
  gridDesci(20) = 255
  close(11)
!EOC
end subroutine readnldas2crd
