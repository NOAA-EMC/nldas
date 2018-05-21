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
#include "misc.h"
!BOP
! 
! !ROUTINE: vic_setup.F90
!
! !DESCRIPTION:
! Completes the setup routines for VIC
!
! !REVISION HISTORY:
! 
! 21 Nov 2003; Sujay Kumar : Initial Specification
! 
! !INTERFACE:
subroutine vic_setup()
!USES:
  use lisdrv_module, only : lis, tile
  use vic_varder,    only : vicdrv
  use spmdMod
#if 0
#if ( defined OPENDAP )
  use opendap_module
#endif
#endif
!EOP
  IMPLICIT NONE
  integer :: t, t1, t2, n, ntype, ierr
  integer :: nc, nr
  real    :: initial_surf_temp
!=== End Variable List ===================================================
!BOC
#if 0
#if ( defined OPENDAP )
  nc = lis%d%lnc
  nr = lis%d%lnr
#else
  nc = lis%d%gnc
  nr = lis%d%gnr
#endif
#endif
#if ( ! defined OPENDAP )
  if(masterproc) then
#endif
     t = len(trim(vicdrv%vic_veglibfile))
     call read_veglib(trim(vicdrv%vic_veglibfile),t, &
          ntype,vicdrv%vic_rootzones)

     t = len(trim(vicdrv%vic_sfile))
     t1 = len(trim(lis%p%safile))
     t2 = len(trim(lis%p%clfile))

     call read_soils()
!     call read_soilparam(trim(vicdrv%vic_sfile),lis%d%nch, &
!          t,vicdrv%vic_nlayer,vicdrv%vic_nnode,            &
!          trim(lis%p%safile),t1,trim(lis%p%clfile),t2,     &
!          nc,nr,tile%col,tile%row)
     
     ! override parameter values with parameter maps
!     call read_parammap()
!     call read_parammap(trim(vicdrv%vic_dsmapfile),          &
!                        len(trim(vicdrv%vic_dsmapfile)),     &
!                        trim(vicdrv%vic_dsmaxmapfile),       &
!                        len(trim(vicdrv%vic_dsmaxmapfile)),  &
!                        trim(vicdrv%vic_wsmapfile),          &
!                        len(trim(vicdrv%vic_wsmapfile)),     &
!                        trim(vicdrv%vic_infiltmapfile),      &
!                        len(trim(vicdrv%vic_infiltmapfile)), &
!                        trim(vicdrv%vic_depth1mapfile),      &
!                        len(trim(vicdrv%vic_depth1mapfile)), &
!                        trim(vicdrv%vic_depth2mapfile),      &
!                        len(trim(vicdrv%vic_depth2mapfile)), &
!                        trim(vicdrv%vic_depth3mapfile),      &
!                        len(trim(vicdrv%vic_depth3mapfile)), &
!                        lis%d%nch,nc,nr,tile%col,tile%row,vicdrv%vic_nlayer)
#if ( ! defined OPENDAP )
  endif
#endif
  
#if ( ( ! defined OPENDAP ) && ( defined SPMD ) )
  if ( npes > 1 ) then
     call MPI_BCAST(ntype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
     call vic_bcast(ntype)   ! veg_lib_struct
     call vic_scatter() ! soil_con_struct
  endif
#endif

  call read_vegparam(tile%vegt)
  call calc_root_fractions(vicdrv%vic_nlayer,vicdrv%vic_rootzones)
!  call make_dist_prcp(vicdrv%vic_snowband)
  call vic_totinit()
  
!  call vic_soilparamalloc(vicdrv%vic_nlayer,vicdrv%vic_nnode, &
!       vicdrv%vic_snowband);

  if ( lis%o%startcode .eq. 3 .or. lis%o%startcode .eq. 2 ) then
    print*,'DBG: vic_setup -- nlayer, snowband, nnode, quick_flux, '// &
            'grnd_flux, frozen_soil', &
            vicdrv%vic_nlayer,        &
            vicdrv%vic_snowband,      &
            vicdrv%vic_nnode,         &
            vicdrv%vic_quick_flux,    &
            vicdrv%vic_grnd_flux,     &
            vicdrv%vic_frozen_soil, ' (',iam, ')'

    ! VIC expects temperatures to be in degrees Celcius
    initial_surf_temp = vicdrv%vic_initial_surf_temp - 273.15
    call initialize_model_state(vicdrv%vic_nlayer,     &
                                vicdrv%vic_snowband,   &
                                vicdrv%vic_nnode,      &
                                vicdrv%vic_quick_flux, &
                                vicdrv%vic_grnd_flux,  &
                                vicdrv%vic_frozen_soil,&
                                initial_surf_temp)
  endif

  print*,'MSG: vic_setup -- Done',' (',iam,')'
!EOC
end subroutine vic_setup
