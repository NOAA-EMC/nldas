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
! !ROUTINE: vic_main.F90
!
! !DESCRIPTION:
!  Initializes vic model state and calls the VIC physics routines
!
! !REVISION HISTORY:
!   14 Apr 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:
subroutine vic_main()
! !USES:
  use lisdrv_module, only: lis,tile 
  use vic_varder, only : vicdrv
  use spmdMod, only : masterproc, npes
!EOP
  IMPLICIT NONE
  integer :: t, n
  integer :: outputflag 
!BOC
  outputflag = 0
! TODO: I don't think this should be here. 
!       It should be called before the first time step and 
!       before the read restart
!  call initialize_model_state(vicdrv%vic_nlayer,   &
!       vicdrv%vic_snowband,vicdrv%vic_nnode,       &
!       vicdrv%vic_quick_flux,vicdrv%vic_grnd_flux, &
!       vicdrv%vic_frozen_soil); 
  if ( lis%t%tscount == 0 .or. lis%t%tscount == 1 ) then 
     outputflag = 1
  endif
  call vic_run(vicdrv%vic_nlayer,      &
               vicdrv%vic_nnode,       &
               vicdrv%vic_snowband,    &
               lis%t%da,               &
               lis%t%doy,              &
               lis%t%hr,               &
               lis%t%mo,               &
               lis%t%yr,               &
               vicdrv%vic_full_energy, &
               vicdrv%vic_frozen_soil, &
               vicdrv%vic_grnd_flux,   &
               vicdrv%vic_quick_flux,  &
               outputflag); 
!EOC
end subroutine vic_main
