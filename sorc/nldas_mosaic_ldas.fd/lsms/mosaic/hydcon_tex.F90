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
subroutine hydcon_tex (nch,tex,ksat)

!	this subroutine assigns the saturated hydraulic conductivity
!	 [m/s] using the values reported by cosby et al. [1984], based
!	 on the texture class.
!	original code by matt rodell (27 june 2001).

  implicit none

  integer :: nch
  integer :: tex(nch)
  real	:: ksat(nch)

  real, parameter :: kc = 0.000000974	! cosby's mean ksat  
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
        
!	assign ksat values.
  do i=1,nch
     select case (tex(i))
     case (1)
        ksat(i) = kc
     case (2)
        ksat(i) = ksic
     case (3) 
        ksat(i) = ksc
     case (4)
        ksat(i) = kcl
     case (5)
        ksat(i) = ksicl
     case (6)
        ksat(i) = kscl
     case (7)
        ksat(i) = kl
     case (8)
        ksat(i) = ksil
     case (9) 
        ksat(i) = ksl
     case (10) 
        ksat(i) = kls
     case (11)
        ksat(i) = ks
        !	    non-land points.
     case default
        ksat(i) = -0.99
     end select
  end do 		!i

  return
end subroutine hydcon_tex

