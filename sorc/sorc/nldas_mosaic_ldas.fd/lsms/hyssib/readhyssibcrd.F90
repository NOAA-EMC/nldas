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
! !ROUTINE : readhyssibcrd.F90
!
! !DESCRIPTION:
!  Routine to read HY-SSiB specific parameters from the card file. 
!
! !REVISION HISTORY:
! 14 Oct 2003: Sujay Kumar, Initial Code
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      subroutine readhyssibcrd(hyssibdrv)
! !USES:
      use hyssibdrv_module
!EOP
      implicit none

      integer :: lsm
      type(hyssibdrvdec) :: hyssibdrv

      namelist /hyssib/hyssibdrv

      open(11,file='lis.crd',form='formatted',status='old')
      read(unit=11,NML=hyssib)
      print *,'Running HY-SSIB LSM:'
      print *,'HY-SSIB Active Restart File: ', hyssibdrv%HYSSIB_RFILE
      hyssibdrv%hyssib_gfractime = 0.0
      hyssibdrv%hyssib_albtime = 0
      hyssibdrv%hyssib_albdchk = 0
      hyssibdrv%hyssib_gfracdchk = 0
      hyssibdrv%HYSSIBOPEN=0

      close(11)

      end subroutine readhyssibcrd

