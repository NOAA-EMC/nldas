	SUBROUTINE PARAMB (NC,NR,TEX,B)

!	This subroutine assigns the dimensionless b parameter based on
!	 values modified from the soil texture class and the table of 
!	 Cosby et al. [1984].
!	Original code by Matt Rodell (11 March 2001).

	implicit none

	integer :: NC,NR
	integer :: TEX(NC,NR)
	real	:: B(NC,NR)

	real, parameter :: bc = 7.62		! Cosby's b values for 
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

	integer :: i,j 

!	Assign b values.
	do j=1,NR
	  do i=1,NC
	    select case (TEX(i,j))
	    case (1)
	      B(i,j) = bc
	    case (2)
	      B(i,j) = bsic
	    case (3) 
	      B(i,j) = bsc
	    case (4)
	      B(i,j) = bcl
	    case (5)
	      B(i,j) = bsicl
	    case (6)
	      B(i,j) = bscl
	    case (7)
	      B(i,j) = bl
	    case (8)
	      B(i,j) = bsil
	    case (9) 
	      B(i,j) = bsl
	    case (10) 
	      B(i,j) = bls
	    case (11)
	      B(i,j) = bs
!	    Non-land points.
	    case default
	      B(i,j) = -0.99
	    end select
	  end do 		!i
	end do			!j

	return
	end

