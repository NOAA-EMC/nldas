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
! !MODULE: vicdrv_module.F90
!
! !DESCRIPTION:
!  Module for runtime specific VIC variables
!
! !REVISION HISTORY:
!
! 14 Oct 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
module vicdrv_module
! !ARGUMENTS:
  type vicdrvdec
     integer :: vic_nlayer     !Number of soil layers in VIC
     integer :: vic_nnode      !Number of soil thermal nodes in the model
     integer :: vic_snowband   !Number of snow bands
     integer :: vic_rootzones 
     integer :: vic_full_energy !FLAG 
     integer :: vic_grnd_flux   !GRND_FLUX 
     integer :: vic_frozen_soil !FROZEN_SOIL
     integer :: vic_quick_flux  !quick_flux
     integer :: vicopen
     character*40:: vic_sfile   !VIC SOIL FILE
     character*40:: vic_veglibfile !VIC veg library
     character*40:: vic_rfile   !VIC restart file
     ! vic parameter maps
     character*40:: VIC_dsmapfile     !VIC Ds map
     character*40:: vic_dsmaxmapfile  !VIC Dsmax map
     character*40:: vic_wsmapfile     !VIC Ws map
     character*40:: vic_infiltmapfile !VIC infilt map
     character*40:: vic_depth1mapfile !VIC depth1 map
     character*40:: vic_depth2mapfile !VIC depth2 map
     character*40:: vic_depth3mapfile !VIC depth3 map
     real :: writeintvic      !Counts number of output times for VIC	
     real :: vic_initial_surf_temp ! Initial surface temperature ( K )
  end type vicdrvdec
!EOP
end module vicdrv_module
