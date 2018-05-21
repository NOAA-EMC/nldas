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
!
! !ROUTINE: noah_alb
!
! !DESCRIPTION:
! This routine determines which type of albedo data is required and calls
! the appropriate reading routine.
!
! !REVISION HISTORY:
! 14 Jun 2005: James Geiger; Initial Specification
!
!
! !INTERFACE:
subroutine noah_alb
! !USES:
   use noah_varder
!EOP
   implicit none
!BOC
   if ( noahdrv%albedo_type == 1 ) then
      call noah_lis_alb
   elseif ( noahdrv%albedo_type == 2 ) then
      call noah_gswp_alb
   else
      call lis_log_msg("ERR: noah_alb -- Don't know how to read albedo.  "// &
                       "Please check the definition of noahdrv%albedo_type.")
      call endrun
   endif

   return
!EOC
end subroutine noah_alb

