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
! !ROUTINE: hyssib_dynsetup.F90
!
! !DESCRIPTION:
!  Updates the time dependent HY-SSiB variables
!
! !REVISION HISTORY:
! 15 Apr 2002: Sujay Kumar, Initial Specification
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
! 
! !INTERFACE:
      subroutine hyssib_dynsetup()
! !USES:
      use lisdrv_module, only: lis,tile 
      use hyssib_varder
      use spmdMod, only : masterproc, npes
      use hyssibpardef_module
!EOP
      implicit none

      integer :: t, n,ier
      if (npes.gt.1) then
         call hyssib_gather
      endif
      if (masterproc) then 
         call hyssib_gfrac(lis%d, lis%t, tile)
!         call hyssib_alb(ldas%d, ldas%t, tile)  
      endif
#if (defined SPMD)
      call MPI_BCAST(hyssibdrv%hyssib_gflag,1,MPI_INTEGER,0, &
          MPI_COMM_WORLD,ier)
      call MPI_BCAST(hyssibdrv%hyssib_aflag,1,MPI_INTEGER,0, &
          MPI_COMM_WORLD,ier)
#endif
      if ((npes.gt.1).and.((hyssibdrv%hyssib_gflag.eq.1).or. &
         (hyssibdrv%hyssib_aflag.eq.1))) then 
         call hyssib_scatter
      endif

      END SUBROUTINE hySSIB_dynSETUP

