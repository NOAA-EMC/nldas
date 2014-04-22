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
subroutine paramb (nch,tex,b)

!	this subroutine assigns the dimensionless b parameter based on
!	 values modified from the soil texture class and the table of 
!	 cosby et al. [1984].
!	original code by matt rodell (11 march 2001).

  implicit none

  integer :: nch
  integer :: tex(nch)
  real	:: b(nch)

  real, parameter :: bc = 7.62		! cosby's b values for 
  real, parameter :: bsic = 6.12		! each texture class, 
  real, parameter :: bsc = 9.19		! *minus one standard
  real, parameter :: bcl = 4.43		! deviation*
  real, parameter :: bsicl = 4.39
  real, parameter :: bscl = 3.38
  real, parameter :: bl = 3.59
  real, parameter :: bsil = 3.61
  real, parameter :: bsl = 3.34
  real, parameter :: bls = 2.31
  real, parameter :: bs = 1.41
  
  integer :: i

!	assign b values.
  do i=1,nch
     select case (tex(i))
     case (1)
        b(i) = bc
     case (2)
        b(i) = bsic
     case (3) 
        b(i) = bsc
     case (4)
        b(i) = bcl
     case (5)
        b(i) = bsicl
     case (6)
        b(i) = bscl
     case (7)
        b(i) = bl
     case (8)
        b(i) = bsil
     case (9) 
        b(i) = bsl
     case (10) 
        b(i) = bls
     case (11)
        b(i) = bs
!	    non-land points.
     case default
        b(i) = -0.99
     end select
  end do 		!i
  return
end subroutine paramb

