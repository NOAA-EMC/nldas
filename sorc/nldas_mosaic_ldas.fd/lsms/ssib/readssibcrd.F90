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
! !ROUTINE: readssibcrd.F90
!
! !DESCRIPTION:
!  Routine to read SSiB specific parameters from the card file. 
!
! !REVISION HISTORY:
! 14 Oct 2003: Sujay Kumar, Initial Code
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
!  5 May 2004: David Mocko, made compatible with SiB-lings
!
! !INTERFACE:
      subroutine readssibcrd(ssibdrv)
! !USES:
      use ssibdrv_module
!EOP
      implicit none

      integer :: lsm
      type(ssibdrvdec) :: ssibdrv

      namelist /ssib/ssibdrv
!BOC
      open(11,file='lis.crd',form='formatted',status='old')
      read(unit=11,NML=ssib)
      print *,'Running SSIB LSM:'
      print *,'SSIB Active Restart File: ', ssibdrv%SSIB_RFILE
      ssibdrv%ssib_gfractime = 0.0
      ssibdrv%ssib_albtime = 0
      ssibdrv%ssib_albdchk = 0
      ssibdrv%ssib_gfracdchk = 0
      ssibdrv%SSIBOPEN=0

      close(11)
      return
!EOC
      end subroutine readssibcrd

