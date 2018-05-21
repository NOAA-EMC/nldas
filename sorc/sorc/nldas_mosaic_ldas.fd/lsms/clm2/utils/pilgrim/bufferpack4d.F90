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
! !IROUTINE: bufferpack4d --- Pack a ghost region into a buffer
!
! !INTERFACE:  
 subroutine bufferpack4d( In, xfrom, xto, yfrom, yto,              &
                          zfrom, zto, Nfrom, Nto, x1,              &
                          x2, y1, y2, z1, z2,                      &
                          n1, n2, nc, buffer )

! !USES:
   use precision
   implicit none

! !INPUT PARAMETERS:
   integer , intent( in )  :: xfrom,xto                 ! X dims
   integer , intent( in )  :: yfrom,yto                 ! Y dims
   integer , intent( in )  :: zfrom,zto                 ! Z dims
   integer , intent( in )  :: nfrom,nto                 ! N dims
   integer , intent( in )  :: nc                        ! total # of tracers
! SJL
   real(r8),intent( in )::in(xfrom:xto,zfrom:zto,nc,yfrom:yto)

   integer , intent( in )  :: x1,x2                     ! X limits
   integer , intent( in )  :: y1,y2                     ! Y limits
   integer , intent( in )  :: z1,z2                     ! Z limits
   integer , intent( in )  :: n1,n2                     ! N limits

! !INPUT/OUTPUT PARAMETERS:
   real(r8), intent( inout ):: buffer(x1:x2,y1:y2,z1:z2,n1:n2)

! !DESCRIPTION:
!
!     This routine puts a 4-D ghost region at the end of a buffer, 
!     packed in the order of the indices.
!
! !LOCAL VARIABLES:
      integer  I, J, K, L
!
! !REVISION HISTORY:
!   99.09.30   Sawyer     Creation
!   99.10.18   Sawyer     FVCCM3 format (no capitalization)
!   00.07.08   Sawyer     Use precision module
!   01.02.02   Sawyer     Removed SGI directives. OpenMP only; free format
!   01.04.27   S-J Lin    Re-indexing input array to be compatible with CCM
!
!EOP
!-----------------------------------------------------------------------
!BOC

!
! Fill the buffer, first in X then in Y
!
   do L = n1, n2

!$omp  parallel do private(i, j, k)

      do k = z1, z2
         do j = y1, y2
            do i = x1, x2
               buffer(i,j,k,L) = in(i,k,L,j)
            enddo
         enddo
      enddo
   enddo
!EOC
 end subroutine bufferpack4d
!-----------------------------------------------------------------------
