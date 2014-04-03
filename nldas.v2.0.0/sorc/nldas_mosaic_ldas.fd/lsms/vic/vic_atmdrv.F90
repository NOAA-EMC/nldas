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
! !ROUTINE: vic_atmdrv.F90
!
! !DESCRIPTION:
!  Transfer forcing from grid to vic tile space.
!
! !REVISION HISTORY:
!   14 Apr 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:
subroutine vic_atmdrv(t, forcing)
! !USES:      
  use lisdrv_module , only : lis     ! LDAS non-model-specific 1-D variables  
!EOP
  IMPLICIT NONE
  real :: forcing(10)
!=== Local Variables =====================================================
  INTEGER :: F,C,R,T,I,J     ! Loop counters
  INTEGER :: NFORCE          ! Number of forcing variables
  integer :: rc, ier              ! for time manager
!=== End Variable Definition =============================================
!BOC
  call vic_f2t(t, forcing,lis%t%ts)
  return
!EOC
end subroutine vic_atmdrv
