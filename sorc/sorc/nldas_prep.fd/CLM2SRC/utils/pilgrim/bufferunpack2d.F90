!-----------------------------------------------------------------------
!BOP
! !IROUTINE: bufferunpack2d --- unpack a ghost region into an array
!
! !INTERFACE:  
      subroutine bufferunpack2d ( buffer, x1, x2, y1, y2,               &
                                  xfrom, xto, yfrom, yto, Out )

! !USES:
      use precision
      implicit none 

! !INPUT PARAMETERS:
      integer , intent( in )  :: x1,x2                        ! X limits
      integer , intent( in )  :: y1,y2                        ! Y limits
      real(r8), intent( in )  :: buffer(x1:x2,y1:y2)          ! Packed buffer
      integer , intent( in )  :: xfrom,xto                    ! X dims
      integer , intent( in )  :: yfrom,yto                    ! Y dim

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent( inout )  :: out(xfrom:xto,yfrom:yto)  ! Out array

! !DESCRIPTION:
!
!     This routine empties the buffer into the appropriate position
!     in the array.
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
          out( I,J ) = buffer( I,J )
        enddo
      enddo
      return
!EOC
      end subroutine bufferunpack2d
!-----------------------------------------------------------------------
