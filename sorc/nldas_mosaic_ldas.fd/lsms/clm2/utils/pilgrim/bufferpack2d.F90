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
!
! !INTERFACE:  
      subroutine bufferpack2d ( In, xfrom, xto, yfrom, yto,             &
                                x1, x2, y1, y2, buffer )

! !USES:
      use precision
      implicit none

! !INPUT PARAMETERS:
      integer , intent( in )     :: xfrom,xto                 ! X dims
      integer , intent( in )     :: yfrom,yto                 ! Y dim
      real(r8), intent( in )     :: in(xfrom:xto,yfrom:yto)   ! In array
      integer , intent( in )     :: x1,x2                     ! X limits
      integer , intent( in )     :: y1,y2                     ! Y limits

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent( inout )  :: buffer(x1:x2,y1:y2)       ! Packed buffer

! !DESCRIPTION:
!
!     This routine puts a 2-D ghost region at the end of a buffer, first in
!     X then in Y.
!
! !LOCAL VARIABLES:
      integer  I, J
!
! !REVISION HISTORY:
!   99.09.30   Sawyer     Creation
!   99.10.18   Sawyer     FVCCM3 format (no capitalization)
!   00.07.08   Sawyer     Use precision module
!   01.02.02   Sawyer     Now free format
!
!EOP
!-----------------------------------------------------------------------
!BOC

!
! Fill the buffer, first in X then in Y
!
      do J = y1, y2
        do I = x1, x2
          buffer( I,J ) = in( I, J )
        enddo
      enddo
      return
!EOC
      end subroutine bufferpack2d
!-----------------------------------------------------------------------
