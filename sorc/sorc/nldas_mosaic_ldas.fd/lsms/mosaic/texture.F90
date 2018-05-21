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
subroutine texture (nc, nr, nch, sand,silt,clay,tex)
!	Finds the USDA texture class based on percentages of sand, silt,
!	 and clay.
!	Original code by Matt Rodell (7 March 2001)
  use lisdrv_module, only : tile
  implicit none
  
  integer :: nc,nr,nch
  real :: sand(nc,nr),silt(nc,nr),clay(nc,nr)
  integer :: tex(nch)
  
  integer :: i,j
  real :: frac,sa,si,cl
  
  TEX = 0

  do i=1,nch
!	    Test for ocean points.
     if (CLAY(tile(i)%col,tile(i)%row) .lt. 0.00) then
        tex(i) = -99
     else
        
!	      Adjust values at points whose fractions don't sum to 1.00.
        frac = CLAY(tile(i)%col,tile(i)%row) + & 
             SAND(tile(i)%col,tile(i)%row) + & 
             SILT(tile(i)%col,tile(i)%row)
        cl = CLAY(tile(i)%col,tile(i)%row) / frac
        sa = SAND(tile(i)%col,tile(i)%row) / frac
        si = SILT(tile(i)%col,tile(i)%row) / frac
        
!	      Identify texture class.
        if ((cl .ge. 0.40) .and. (si .lt. 0.40) .and.  &
             (sa .lt. 0.44)) then
           tex(i) = 1	!CLAY
        end if
        
        if ((cl .ge. 0.40) .and. (si .ge. 0.40)) then
           tex(i) = 2	!SILTY CLAY
        end if
        
        if ((cl .ge. 0.36) .and. (sa .ge. 0.44)) then
           tex(i) = 3	!SANDY CLAY
        end if
        
        if ((cl .ge. 0.28) .and. (cl .lt. 0.40) .and.  &
             (sa .ge. 0.20) .and. (sa .lt. 0.44)) then
           tex(i) = 4	!CLAY LOAM
        end if
        
        if ((cl .ge. 0.28) .and. (cl .lt. 0.40) .and.  &
             (sa .lt. 0.20)) then
           tex(i) = 5	!SILTY CLAY LOAM
        end if
        
        if ((cl .ge. 0.20) .and. (cl .lt. 0.36) .and. &
             (si .lt. 0.28) .and. (sa .ge. 0.44)) then
           tex(i) = 6	!SANDY CLAY LOAM
        end if
        
        if ((cl .ge. 0.08) .and. (cl .lt. 0.28) .and.  &
             (si .ge. 0.28) .and. (si .lt. 0.50) &
             .and. (sa .lt. 0.52)) then
           tex(i) = 7	!LOAM
        end if
        
        if ((cl .lt. 0.28) .and. (si .ge. 0.50)) then
           tex(i) = 8	!SILTY LOAM (and SILT)
        end if
        
        if ((sa - cl) .gt. 0.70) then
           if ((sa - (0.3 * cl)) .gt. 0.87) then
              tex(i) = 11	!SAND
           else
              tex(i) = 10	!LOAMY SAND
           end if
        end if
        
        if (tex(i) .eq. 0) then
           tex(i) = 9	!SANDY LOAM
        end if
        
     end if	!CLAY
  end do 	!i

  return
end subroutine texture
