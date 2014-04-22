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
! !ROUTINE: opendap_init_c_struct
! 
! !DESCRIPTION: 
! This routine assigns various opendap variables from Fortran to their C
! counterpart.
! 
! !REVISION HISTORY: 
! 14 Oct 2003; James Geiger :Initial Specification
! 
! !INTERFACE:
subroutine opendap_init_c_struct
! !USES:
#if ( defined OPENDAP )
  use spmdMod, only : iam
  use opendap_module, only : parm_nc, parm_nr,     &
       parm_slat, parm_nlat, &
       parm_wlon, parm_elon, &
       tnroffset
!EOP
  implicit none
!BOC
  call setup_c_struct(iam,                  &
       parm_nc, parm_nr,     &
       parm_slat, parm_nlat, &
       parm_wlon, parm_elon, &
       tnroffset)
!EOC
#endif
end subroutine opendap_init_c_struct
