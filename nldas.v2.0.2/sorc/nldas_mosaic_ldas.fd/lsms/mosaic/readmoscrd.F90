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
! !ROUTINE : readmoscrd.F90
!
! !DESCRIPTION:
!  Routine to read Mosaic specific parameters from the card file. 
!
! !REVISION HISTORY:
! 15 Oct 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
   subroutine readmoscrd(mosdrv)
! !USES:
      use mosdrv_module
!EOP
      implicit none
      integer :: lsm
      type(mosdrvdec) :: mosdrv
      namelist /mos/mosdrv
      open(11,file='lis.crd',form='formatted',status='old')
      read(unit=11,NML=mos)
      print*,'Running MOS LSM:'
      print*,'MOS Active Restart File: ', mosdrv%MOS_RFILE
      mosdrv%MOSopen = 0
      close(11)
    end subroutine readmoscrd
