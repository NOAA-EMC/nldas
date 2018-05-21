	SUBROUTINE HYDCON_TEX (NC,NR,TEX,KSAT)

!	This subroutine assigns the saturated hydraulic conductivity
!	 [m/s] using the values reported by Cosby et al. [1984], based
!	 on the texture class.
!	Original code by Matt Rodell (27 June 2001).

	implicit none

	integer :: NC,NR
	integer :: TEX(NC,NR)
	real	:: KSAT(NC,NR)

	real, parameter :: kc = 0.000000974	! Cosby's mean ksat  
	real, parameter :: ksic = 0.00000134	! values [m/s] for 
	real, parameter :: ksc = 0.00000722	! each texture class
	real, parameter :: kcl = 0.00000245		
	real, parameter :: ksicl = 0.00000203
	real, parameter :: kscl = 0.00000445
	real, parameter :: kl = 0.00000338
	real, parameter :: ksil = 0.00000281
	real, parameter :: ksl = 0.00000523
	real, parameter :: kls = 0.0000141
	real, parameter :: ks = 0.0000466

	integer :: i,j 

!	Assign KSAT values.
	do j=1,NR
	  do i=1,NC
	    select case (TEX(i,j))
	    case (1)
	      KSAT(i,j) = kc
	    case (2)
	      KSAT(i,j) = ksic
	    case (3) 
	      KSAT(i,j) = ksc
	    case (4)
	      KSAT(i,j) = kcl
	    case (5)
	      KSAT(i,j) = ksicl
	    case (6)
	      KSAT(i,j) = kscl
	    case (7)
	      KSAT(i,j) = kl
	    case (8)
	      KSAT(i,j) = ksil
	    case (9) 
	      KSAT(i,j) = ksl
	    case (10) 
	      KSAT(i,j) = kls
	    case (11)
	      KSAT(i,j) = ks
!	    Non-land points.
	    case default
	      KSAT(i,j) = -0.99
	    end select
	  end do 		!i
	end do			!j

	return
	end

