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
! !ROUTINE: ssib_dynsetup.F90
!
! !DESCRIPTION:
!  Updates the time dependent SSiB variables
!
! !REVISION HISTORY:
! 15 Apr 2002: Sujay Kumar, Initial Specification
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
!  6 May 2004: David Mocko, made compatible with SiB-lings
! 
! !INTERFACE:
      subroutine ssib_dynsetup
! !USES:
      use lisdrv_module, only: lis,tile 
      use ssib_varder
      use spmdMod, only : masterproc, npes
      use ssibpardef_module
!EOP
      implicit none

      integer :: t, n,ier
!BOC
#if ( ! defined OPENDAP )
      if (npes.gt.1) then
         call ssib_gather
      endif
      if (masterproc) then
#endif
         call ssib_gfrac
!         call ssib_alb(lis%d, lis%t, tile)  
#if ( ! defined OPENDAP )
      endif
#if (defined SPMD)
      call MPI_BCAST(ssibdrv,1,MPI_SSIBDRV_STRUCT,0, &
          MPI_COMM_WORLD,ier)
#endif
      if ((npes.gt.1).and.((ssibdrv%ssib_gflag.eq.1).or. &
         (ssibdrv%ssib_aflag.eq.1))) then 
         call ssib_scatter
      endif
#endif
      return
!EOC
      end subroutine ssib_dynsetup

