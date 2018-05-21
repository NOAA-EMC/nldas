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
! !IROUTINE: bufferpack2d --- Pack a ghost region into a buffer
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: bufferunpack3d --- Pack a ghost region into a buffer
!
! !INTERFACE:  
      subroutine bufferunpack3d( buffer, x1, x2, y1, y2,                &
                                 z1, z2, xfrom, xto, yfrom,             &
                                 yto, zfrom, zto, Out)


! !USES:
      use precision
      implicit none

! !INPUT PARAMETERS:
      integer , intent( in )  :: x1,x2                     ! X limits
      integer , intent( in )  :: y1,y2                     ! Y limits
      integer , intent( in )  :: z1,z2                     ! Z limits
      real(r8), intent( in )  :: buffer(x1:x2,y1:y2,z1:z2) ! Packed buffer
      integer , intent( in )  :: xfrom,xto                 ! X dims
      integer , intent( in )  :: yfrom,yto                 ! Y dim
      integer , intent( in )  :: zfrom,zto                 ! Z dim

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent( inout )  :: out(xfrom:xto,yfrom:yto,zfrom:zto) ! Out

! !DESCRIPTION:
!
!     This routine puts a 3-D ghost region at the end of a buffer, first in
!     X then in Y then in Z.
!
! !LOCAL VARIABLES:
      integer  I, J, K
!
! !REVISION HISTORY:
!   99.09.30   Sawyer     Creation
!   99.10.18   Sawyer     FVCCM3 format (no capitalization)
!   00.07.08   Sawyer     Use precision module
!   01.02.02   Sawyer     Removed SGI directives. OpenMP only; free format
!
!EOP
!-----------------------------------------------------------------------
!BOC

!$omp  parallel do&
!$omp& default(shared)&
!$omp& private(i,j,k)
      do K = z1, z2
        do J = y1, y2
          do I = x1, x2
            out( I,J,K ) = buffer( I, J, K )
          enddo
        enddo
      enddo
      return
!EOC
      end subroutine bufferunpack3d
!-----------------------------------------------------------------------

