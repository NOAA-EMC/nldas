	SUBROUTINE TEXTURE (NC,NR,SAND,SILT,CLAY,TEX)
!	Finds the USDA texture class based on percentages of sand, silt,
!	 and clay.
!	Original code by Matt Rodell (7 March 2001)

	implicit none

	integer :: NC,NR
	real :: SAND(NC,NR),SILT(NC,NR),CLAY(NC,NR)
	integer :: TEX(NC,NR)

	integer :: i,j
	real :: frac,sa,si,cl

	TEX = 0

	do j=1,NR
	  do i=1,NC

!	    Test for ocean points.
	    if (CLAY(i,j) .lt. 0.00) then
	      TEX(i,j) = -99
	    else

!	      Adjust values at points whose fractions don't sum to 1.00.
		frac = CLAY(i,j) + SAND(i,j) + SILT(i,j)
		cl = CLAY(i,j) / frac
		sa = SAND(i,j) / frac
		si = SILT(i,j) / frac

!	      Identify texture class.
	      if ((cl .ge. 0.40) .and. (si .lt. 0.40) .and. 
     &	       (sa .lt. 0.44)) then
	      	TEX(i,j) = 1	!CLAY
	      end if

	      if ((cl .ge. 0.40) .and. (si .ge. 0.40)) then
	      	TEX(i,j) = 2	!SILTY CLAY
	      end if

	      if ((cl .ge. 0.36) .and. (sa .ge. 0.44)) then
	      	TEX(i,j) = 3	!SANDY CLAY
	      end if

	      if ((cl .ge. 0.28) .and. (cl .lt. 0.40) .and. 
     &	       (sa .ge. 0.20) .and. (sa .lt. 0.44)) then
	      	TEX(i,j) = 4	!CLAY LOAM
	      end if

	      if ((cl .ge. 0.28) .and. (cl .lt. 0.40) .and. 
     &	       (sa .lt. 0.20)) then
	      	TEX(i,j) = 5	!SILTY CLAY LOAM
	      end if

	      if ((cl .ge. 0.20) .and. (cl .lt. 0.36) .and. 
     &	       (si .lt. 0.28) .and. (sa .ge. 0.44)) then
	      	TEX(i,j) = 6	!SANDY CLAY LOAM
	      end if

	      if ((cl .ge. 0.08) .and. (cl .lt. 0.28) .and. 
     &	       (si .ge. 0.28) .and. (si .lt. 0.50) 
     &	       .and. (sa .lt. 0.52)) then
	      	TEX(i,j) = 7	!LOAM
	      end if

	      if ((cl .lt. 0.28) .and. (si .ge. 0.50)) then
	      	TEX(i,j) = 8	!SILTY LOAM (and SILT)
	      end if

	      if ((sa - cl) .gt. 0.70) then
	      	if ((sa - (0.3 * cl)) .gt. 0.87) then
		  TEX(i,j) = 11	!SAND
	      	else
		  TEX(i,j) = 10	!LOAMY SAND
	      	end if
	      end if

	      if (TEX(i,j) .eq. 0) then
	      	TEX(i,j) = 9	!SANDY LOAM
	      end if

	    end if	!CLAY

	  end do 	!i
	end do		!j

	return
	end
