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
! !ROUTINE: noah_atmdrv.F90
!
! !DESCRIPTION:
!  Transfer forcing from grid to tile space.
!
! !REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Code
!  28 Jan 2002: Jon Gottschalck; Added option for different number of forcing variables
!
! !INTERFACE:
    subroutine mos_f2t(t, forcing)
! !USES:      
      use lisdrv_module , only : lis,tile     ! LDAS non-model-specific 1-D variables
      use spmdMod
      use tile_spmdMod
      use mos_varder
!EOP     
      implicit none
      real :: forcing(10)
      INTEGER :: F,C,R,T,I,J     ! Loop counters
      INTEGER :: NFORCE          ! Number of forcing variables
      integer :: rc, ier              ! for time manager
!BOC      
      do f=1,lis%f%nforce
        mos(t)%forcing(f)=forcing(f)
!        if (forcing(f).eq.564.072021) print *,'564 mos t=',t
!        if (forcing(f).eq.604.538574) print *,'604 mos t=',t
!        if (forcing(f).eq.53.7600021) stop
      enddo
      return
  end subroutine mos_f2t
