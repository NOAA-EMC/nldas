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
! !ROUTINE: pack_string_f2c
! 
! !DESCRIPTION: 
! This routine packs an array of strings into a single single that may
! be passed into C routines.
! 
! !REVISION HISTORY: 
! 02 Mar 2004; James Geiger :Initial Specification
! 
! !INTERFACE:
subroutine pack_string_f2c(n, string_array, string)

   implicit none

! !INPUT PARAMETERS:
   integer, intent(in)                        :: n
   character(len=*), dimension(n), intent(in) :: string_array

! !OUTPUT PARAMETERS:
   character(len=*), intent(out)              :: string
!EOP

   integer :: i, s_index, e_index, length

   length = len(string)

!BOC

   s_index = 1
   e_index = len(trim(adjustl(string_array(1)))) + 1
   do i = 1, n-1
      string(s_index:e_index) = trim(adjustl(string_array(i))) // ' '
      s_index = e_index + 1
      e_index = s_index + len(trim(adjustl(string_array(i+1))))
   enddo
   string(s_index:e_index) = trim(adjustl(string_array(n))) // ' '

   string(e_index+1:e_index+1) = char(0)

!EOC
end subroutine pack_string_f2c
