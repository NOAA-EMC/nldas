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
subroutine psis (nch,tex,psi)

!	this subroutine assigns the saturated soil potential based
!	 on the soil texture class and the table of cosby et al. 
!	 [1984] [m].
!	original code by matt rodell (11 march 2001)

  implicit none

  integer :: nch
  integer :: tex(nch)
  real :: psi(nch)

  real, parameter :: pc = -0.468		! cosby's psi values for 
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
  
  integer :: i

!	assign saturated soil potential values.
  do i=1,nch
     select case (tex(i))
     case (1)
        psi(i) = pc
     case (2)
        psi(i) = psic
     case (3) 
        psi(i) = psc
     case (4)
        psi(i) = pcl
     case (5)
        psi(i) = psicl
     case (6)
        psi(i) = pscl
     case (7)
        psi(i) = pl
     case (8)
        psi(i) = psil
     case (9) 
        psi(i) = psl
     case (10) 
        psi(i) = pls
     case (11)
        psi(i) = ps
!	    non-land points.
     case default
        psi(i) = -0.99
     end select
  end do 		!i
  
  return
end subroutine psis

