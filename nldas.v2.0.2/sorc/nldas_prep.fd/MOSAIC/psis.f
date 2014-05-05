	SUBROUTINE PSIS (NC,NR,TEX,PSI)

!	This subroutine assigns the saturated soil potential based
!	 on the soil texture class and the table of Cosby et al. 
!	 [1984] [m].
!	Original code by Matt Rodell (11 March 2001)

	implicit none

	integer :: NC,NR
	integer :: TEX(NC,NR)
	real :: PSI(NC,NR)

	real, parameter :: pc = -0.468		! Cosby's psi values for 
	real, parameter :: psic = -0.324	! each texture class [m].
	real, parameter :: psc = -0.098
	real, parameter :: pcl = -0.263
	real, parameter :: psicl = -0.617
	real, parameter :: pscl = -0.135
	real, parameter :: pl = -0.355
	real, parameter :: psil = -0.759
	real, parameter :: psl = -0.141
	real, parameter :: pls = -0.036
	real, parameter :: ps = -0.069

	integer :: i,j 

!	Assign saturated soil potential values.
	do j=1,NR
	  do i=1,NC
	    select case (TEX(i,j))
	    case (1)
	      PSI(i,j) = pc
	    case (2)
	      PSI(i,j) = psic
	    case (3) 
	      PSI(i,j) = psc
	    case (4)
	      PSI(i,j) = pcl
	    case (5)
	      PSI(i,j) = psicl
	    case (6)
	      PSI(i,j) = pscl
	    case (7)
	      PSI(i,j) = pl
	    case (8)
	      PSI(i,j) = psil
	    case (9) 
	      PSI(i,j) = psl
	    case (10) 
	      PSI(i,j) = pls
	    case (11)
	      PSI(i,j) = ps
!	    Non-land points.
	    case default
	      PSI(i,j) = -0.99
	    end select
	  end do 		!i
	end do			!j

	return
	end

