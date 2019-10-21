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
! !MODULE: persdomain_module.F90
! 
! !DESCRIPTION: 
!  Contains routines and variables that define the native domain
!  for PERS precipitation product.
! 
! !INTERFACE:
module persdomain_module
! !USES:
  use persdrv_module
! !ARGUMENTS:
  type(persdrvdec) :: persdrv
  integer :: mi
  real, allocatable :: rlat(:)
  real, allocatable :: rlon(:)
  integer, allocatable :: n11(:)
  integer, allocatable :: n12(:)
  integer, allocatable :: n21(:)
  integer, allocatable :: n22(:)
  real, allocatable ::  w11(:),w12(:)
  real, allocatable ::  w21(:),w22(:)
!EOP  

contains
  
!BOP
!
! !ROUTINE: defnatpers.F90
! 
! !DESCRIPTION: 
!  Defines the kgds array describing the native forcing resolution 
!  for PERS data. 
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine defnatpers()
! !USES: 
    use lisdrv_module, only: lis
    implicit none
! !ARGUMENTS:
    integer :: kgdsi(200)
!EOP
!BOC
    call readperscrd(persdrv)
    kgdsi = 0
    call allocate_pers_ip(lis%d%lnc*lis%d%lnr)
!EOC
  end subroutine defnatpers
!BOP
! !ROUTINE: allocate_pers_ip
!
! !DESCRIPTION: 
! 
! Allocates memory for PERS interpolation variables
! 
! !INTERFACE:
  subroutine allocate_pers_ip(N)
!EOP
    integer :: N
!BOC
    allocate(rlat(n))
    allocate(rlon(n))
    allocate(n11(n))
    allocate(n12(n))
    allocate(n21(n))
    allocate(n22(n))
    allocate(w11(n))
    allocate(w12(n))
    allocate(w21(n))
    allocate(w22(n))
    mo = n
    nn = n
    w11 = 0.0
    w12 = 0.0
    w21 = 0.0
    w22 = 0.0
!EOC
  end subroutine allocate_pers_ip
end module persdomain_module
