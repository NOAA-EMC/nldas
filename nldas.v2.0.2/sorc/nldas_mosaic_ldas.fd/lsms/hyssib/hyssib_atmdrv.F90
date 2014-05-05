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
! !ROUTINE: hyssib_atmdrv.f
!
! !DESCRIPTION:
!  Transfer forcing from grid to tile space
!
! !REVISION HISTORY:
! 15 Oct 1999: Paul Houser, Initial Code
! 28 Jan 2002: Jon Gottschalck, Added option for different number of
!                               forcing variables
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      SUBROUTINE hyssib_f2t(t, forcing)
! !USES:
      use lisdrv_module, only : lis   ! LIS non-model-specific 1-D variables
      use spmdMod
      use tile_spmdMod
      use hyssib_varder
!EOP     
      implicit none

      real :: forcing(10)
      INTEGER :: F,C,R,T,I,J     ! Loop counters
      INTEGER :: NFORCE          ! Number of forcing variables
      integer :: rc, ier              ! for time manager
      DO F=1,lis%f%NFORCE
         hyssib(t)%forcing(f)=forcing(f)
      enddo
!in parallel mode, scatter/bcast the forcing array? 
!need to put them in the right tiles? 
!      endif

      return
      end

