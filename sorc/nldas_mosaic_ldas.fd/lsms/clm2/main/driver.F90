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

subroutine driver ()

!-----------------------------------------------------------------------
!
! Purpose:
! clm model driver
!
! Method:
! Calling sequence:
!
! -> histend   Determines if current time step is the end of history interval
!
! -> calendr   Generate the calendar day (1.00 -> 365.99), month (1 -> 12),
!              and day (1 -> 31) used to calculate the surface albedos and
!              leaf and stem areas for the next time step
!
! -> loop over patch points calling for each patch point:
!    -> Hydrology1          canopy interception and precip on ground
!    -> Biogeophysics1      leaf temperature and surface fluxes
!    -> Biogeophysics_Lake  lake temperature and surface fluxes
!    -> Biogeophysics2      soil/snow and ground temp and update surface fluxes
!    -> Hydrology2          surface and soil hydrology
!    -> Hydrology_Lake      lake hydrology
!    -> Biogeochemistry     surface biogeochemical fluxes (LSM)
!    -> EcosystemDyn:       ecosystem dynamics: phenology, vegetation, soil carbon
!    -> SurfaceAlbedo:      albedos for next time step
!      -> SnowAlbedo:       snow albedos: direct beam
!      -> SnowAlbedo:       snow albedos: diffuse
!      -> SoilAlbedo:       soil/lake albedos
!      -> TwoStream:        absorbed, reflected, transmitted solar fluxes (vis dir)
!      -> TwoStream:        absorbed, reflected, transmitted solar fluxes (vis dif)
!      -> TwoStream:        absorbed, reflected, transmitted solar fluxes (nir dir)
!      -> TwoStream:        absorbed, reflected, transmitted solar fluxes (nir dif)
!    -> BalanceCheck        check for errors in energy and water balances
!    -> histUpdate:         accumulate history fields over history time interval
!
!  -> Rtmriverflux          calls RTM river routing model
!
!  -> histHandler           write history and restart files if appropriate
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id: driver.F90,v 1.17 2004/11/30 20:35:12 jim Exp $
!-----------------------------------------------------------------------
#include "misc.h"
  use precision
  use clm_varder
  use clm_varcon     , only : doalb, eccen, obliqr, lambm0, mvelpp, & 
       denh2o, denice, hvap, hsub, hfus, istwet
  use lisdrv_module, only : lis
!<debug>
!  use clm_varpar    , only : maxpatch
!</debug>
  use clm_varmap    , only : begpatch, endpatch, numpatch, numland, landvec, patchvec
  use clm_varctl    , only : fsurdat, wrtdia, csm_doflxave
  use histHandlerMod, only : histHandler,histend,do_restwrite,lrestwrt,mcdate_i,mcsec_i,mdcur_i,mscur_i
  use histFileMod, only : nhist, nbeghis
!  use restFileMod   , only : restwrt 
  use inicFileMod   , only : inicwrt, do_inicwrite 
  use mvegFileMod   , only : interpmonthlyveg
  use time_manager  , only : get_step_size, get_curr_calday, get_curr_date, get_nstep
#if (defined RTM)
  use RtmMod        , only : Rtmriverflux
#endif
#if (defined SPMD)
  use spmdMod       , only : masterproc, npes, compute_mpigs_patch
  use mpishorthand  , only : mpir8, mpilog, mpicom 
#else
  use spmdMod       , only : masterproc
#endif
#if (defined COUP_CSM)
  use clm_csmMod    , only : csm_dosndrcv, csm_recv, csm_send, csm_flxave, &
                             dorecv, dosend, csmstop_now
#endif
  use shr_sys_mod   , only : shr_sys_flush
  implicit none
  real :: soilm(endpatch,1:nlevsoi)
  real :: soilmtc(endpatch)
  

! ------------------- arguments -----------------------------------
!  logical , intent(in) :: doalb   !true if time for surface albedo calculation
!  real(r8), intent(in) :: eccen   !Earth's orbital eccentricity
!  real(r8), intent(in) :: obliqr  !Earth's obliquity in radians
!  real(r8), intent(in) :: lambm0  !Mean longitude of perihelion at the vernal equinox (radians)
!  real(r8), intent(in) :: mvelpp  !Earth's moving vernal equinox long. of perihelion + pi (radians)
! -----------------------------------------------------------------

! ---------------------- local variables --------------------------
  integer  :: j,k,t,m           !loop/array indices
!  real :: tvegb(144,76)
  integer  :: yrp1                !year (0, ...) for nstep+1
  integer  :: monp1               !month (1, ..., 12) for nstep+1
  integer  :: dayp1               !day of month (1, ..., 31) for nstep+1
  integer  :: secp1               !seconds into current date for nstep+1   
  real(r8) :: caldayp1            !calendar day for nstep+1
  integer  :: dtime               !timestep size [seconds]
  integer  :: nstep               !timestep index
!  real(r8) :: buf1d(numpatch)     !temporary buffer 
!  real(r8) :: tsxyav              !average ts for diagnostic output
#if (defined SPMD)
!  integer :: numrecvv(0:npes-1)   !vector of items to be received  
!  integer :: displsv(0:npes-1)    !displacement vector
!  integer :: numsend              !number of items to be sent
  integer :: ier                  !error code
#endif
  real(r8) :: cgrnd
  real(r8) :: cgrndl
  real(r8) :: cgrnds
  real(r8) :: tg
  real(r8) :: emg
  real(r8) :: htvp
  real(r8) :: dlrad
  real(r8) :: ulrad
  real(r8) :: tssbef(-5:10)
! -----------------------------------------------------------------
!  call t_startf('clm_driver')

! determine time step
  if(masterproc) then
     nstep = get_nstep(lis%t) 
! ----------------------------------------------------------------------
! Coupler receive
! ----------------------------------------------------------------------

#if (defined COUP_CSM)
!
! Determine if information should be sent/received to/from flux coupler
!
  call csm_dosndrcv(doalb)
!
! Get atmospheric state and fluxes from flux coupler
!
  if (dorecv) then
     call csm_recv()
     if (csmstop_now) then
!        call t_stopf('clm_driver')
        RETURN
     endif
  endif
#endif

! ----------------------------------------------------------------------
! Determine if end of history interval
! ----------------------------------------------------------------------

!  call histend ()

! ----------------------------------------------------------------------
! Calendar information for next time step
! o caldayp1 = calendar day (1.00 -> 365.99) for cosine solar zenith angle
!              calday is based on Greenwich time
! o monp1    = month (1 -> 12) for leaf area index and stem area index
! o dayp1    = day (1 -> 31)   for leaf area index and stem area index
! ----------------------------------------------------------------------

  dtime = get_step_size(lis%t)
  caldayp1 = get_curr_calday(lis%t,offset=dtime)  
!  call get_curr_date(yrp1, monp1, dayp1, secp1, offset=dtime)
endif

#if(defined SPMD)
  call MPI_BCAST(nstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(doalb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(caldayp1,1,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(dtime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
#endif


#if (!defined DGVM)
! ----------------------------------------------------------------------
! Determine weights for time interpolation of monthly vegetation data.
! This also determines whether it is time to read new monthly vegetation and
! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
! weights obtained here are used in subroutine ecosystemdyn to obtain time
! interpolated values.
! ----------------------------------------------------------------------
!===LDAS modification: Commented out because LAI etc data read in elsewhere
!!  if (doalb) call interpMonthlyVeg (fsurdat, monp1, dayp1)
#endif

! ----------------------------------------------------------------------
! LOOP 1
! ----------------------------------------------------------------------
!  call t_startf('clm_loop1')
  do k = begpatch, endpatch    ! begin 1st loop over patches
! Initial set of variables
!     print*,iam,clm(12)%t_veg
     clm(k)%nstep = nstep
     clm(k)%h2osno_old = clm(k)%h2osno  ! snow mass at previous time step
     clm(k)%frac_veg_nosno = clm(k)%frac_veg_nosno_alb
! Determine if will cap snow
     if (clm(k)%h2osno > 1000.) then
        clm(k)%do_capsnow = .true.
     else
        clm(k)%do_capsnow = .false.
     endif
!
! Energy for non-lake points
!
     if (.not. clm(k)%lakpoi) then
!
! Initial set of previous time step variables
!       print*,'lp2',iam,k,clm(k)%t_soisno(clm%snl+1) 
        do j = clm(k)%snl+1, 0       ! ice fraction of snow at previous time step
           clm(k)%frac_iceold(j) = clm(k)%h2osoi_ice(j)/(clm(k)%h2osoi_liq(j)+clm(k)%h2osoi_ice(j))
        enddo
!
! Determine beginning water balance (water balance at previous time step)
!
        clm(k)%begwb = clm(k)%h2ocan + clm(k)%h2osno
        do j = 1, nlevsoi
           clm(k)%begwb = clm(k)%begwb + clm(k)%h2osoi_ice(j) + clm(k)%h2osoi_liq(j)

        enddo
! Determine canopy interception and precipitation onto ground surface.
! Determine the fraction of foliage covered by water and the fraction
! of foliage that is dry and transpiring. Initialize snow layer if the
! snow accumulation exceeds 10 mm.
        call Hydrology1(clm(k))
! Determine leaf temperature and surface fluxes based on ground
! temperature from previous time step.
!

        call Biogeophysics1(clm(k),cgrnd,cgrndl,cgrnds,tg,emg,htvp, dlrad,ulrad,tssbef)
!        if(clm(k)%kpatch == 841 .OR. clm(k)%kpatch==1559) then 
!        print*,iam,k,dtime
!        print*,'aft',iam,k,clm(k)%fsa,clm(k)%forc_solad,clm(k)%fabd,clm(k)%fabi
!        print*,iam,k,clm(k)%albgri,clm(k)%albgrd
!        endif
     else if (clm(k)%lakpoi) then
!
! Determine lake temperature and surface fluxes
!
        call Biogeophysics_Lake (clm(k))

     endif
 
 44   format(i6,1x,i4,1x,i4,1x,i2,1x,f6.2,1x,4(f8.4,1x),f12.8,1x,f12.8)
     if (.not. clm(k)%lakpoi) then
!
! Surface biogeochemical fluxes: co2 respiration and plant production
!
#if (defined BGC)
        print*,'**** calling biogeochem...'
        call Biogeochemistry (clm(k))
#endif
!
! Ecosystem dynamics: phenology, vegetation, soil carbon.
! Also updates snow fraction
!
        call EcosystemDyn (clm(k), doalb, .false.)
     else if (clm(k)%lakpoi) then

!        clm(k)%fpsn = 0.
!        clm(k)%frmf = 0. !put these lines here to avoid psn = NaN

     endif
!
! Albedos for next time step
!

     if (doalb) then
!      if(clm(k)%kpatch==1559) then 
!        print*, eccen,mvelpp,lambm0,obliqr
!      endif
        call SurfaceAlbedo (clm(k), caldayp1, eccen, obliqr, lambm0, mvelpp)
     endif
!
! THIS WILL EVENTUALLY MARK THE END OF THE PATCH LOOP AND
! THE BEGINNING OF THE SINGLE COLUMN SOIL LOOP(S)
!
! Determine soil/snow temperatures including ground temperature and
! update surface fluxes for new ground temperature.
!
     if (.not. clm(k)%lakpoi) then
        call Biogeophysics2(clm(k),cgrnd,cgrndl,cgrnds,tg,emg,htvp, dlrad,ulrad,tssbef)
     endif
  end do     
!  call t_stopf('clm_loop1')
! ----------------------------------------------------------------------
! Coupler send
! ----------------------------------------------------------------------

#if (defined COUP_CSM)
!
! Average fluxes over interval if appropriate
! Surface states sent to the flux coupler states are not time averaged
!
     if (csm_doflxave) call csm_flxave()
!
! Send fields to flux coupler
! Send states[n] (except for snow[n-1]), time averaged fluxes for [n,n-1,n-2],
! albedos[n+1], and ocnrof_vec[n]
!
     if (dosend) call csm_send()

#endif

! ----------------------------------------------------------------------
! LOOP 2
! ----------------------------------------------------------------------

!  call t_startf('clm_loop2')
  do k = begpatch, endpatch   ! begin 2nd loop over patches
! Vertical (column) soil and surface hydrology
!
     if (.not. clm(k)%lakpoi) call Hydrology2 (clm(k))
     
! Lake hydrology
!
     if (clm(k)%lakpoi) call Hydrology_Lake (clm(k))
     
! Update Snow Age (needed for surface albedo calculation - but is
! really a column type property
!
     call SnowAge (clm(k))
     
! Fraction of soil covered by snow - really a column property
!
     clm(k)%frac_sno = clm(k)%snowdp/(0.1 + clm(k)%snowdp)
     
!
! Check the energy and water balance
!
     call BalanceCheck (clm(k))
!      print*,'sf..',iam, clm(k)%fabd(1),clm(k)%fabd(2)
!     print*,'dr',iam, k , clm(k)%t_grnd
  end do    
!  call t_stopf('clm_loop2')

! ----------------------------------------------------------------------
! Write global average diagnostics to standard output
! ----------------------------------------------------------------------
  do k=begpatch, endpatch
   
      clm(k)%count=clm(k)%count+1
      clm(k)%totfsa=clm(k)%totfsa+clm(k)%fsa
      clm(k)%toteflx_lwrad_net=clm(k)%toteflx_lwrad_net+ & 
                 clm(k)%eflx_lwrad_net
      clm(k)%toteflx_lh_tot=clm(k)%toteflx_lh_tot+clm(k)%eflx_lh_tot
      clm(k)%toteflx_sh_tot=clm(k)%toteflx_sh_tot+clm(k)%eflx_sh_tot
      clm(k)%toteflx_soil_grnd=clm(k)%toteflx_soil_grnd+ &
       clm(k)%eflx_soil_grnd
      clm(k)%totqflx_snomelt=clm(k)%totqflx_snomelt+clm(k)%qflx_snomelt
      clm(k)%totsnow=clm(k)%totsnow+clm(k)%forc_snow
      clm(k)%totrain=clm(k)%totrain+clm(k)%forc_rain
      clm(k)%totqflx_evap=clm(k)%totqflx_evap+clm(k)%qflx_evap_tot
      clm(k)%totqflx_surf=clm(k)%totqflx_surf + &
                          clm(k)%qflx_surf + clm(k)%qflx_qrgwl
      clm(k)%totqflx_drain=clm(k)%totqflx_drain + clm(k)%qflx_drain

      clm(k)%totqflx_ecanop=clm(k)%totqflx_ecanop+ &          
       (hvap*(clm(k)%qflx_evap_veg-clm(k)%qflx_tran_veg))
      clm(k)%totqflx_tran_veg=clm(k)%totqflx_tran_veg+ &
       clm(k)%qflx_tran_veg
      clm(k)%totqflx_evap_grnd=clm(k)%totqflx_evap_grnd+ &
       clm(k)%qflx_evap_grnd
       clm(k)%totqflx_sub_snow=clm(k)%totqflx_sub_snow+ &
            clm(k)%qflx_sub_snow*(-hsub)
       clm(k)%canopint = clm(k)%h2ocan
  enddo
  soilmtc = 0.0
  if(lis%t%tscount==0 .or. lis%t%tscount ==1 &
       .or. lis%f%rstflag.eq.1 ) then 
     do m=1,nlevsoi 
        do t=begpatch, endpatch
           soilm(t,m)=clm(t)%h2osoi_liq(m)+clm(t)%h2osoi_ice(m)
        enddo
     enddo
     do m=1,nlevsoi
        do t=begpatch, endpatch
           soilmtc(t)=soilmtc(t)+soilm(t,m)
        enddo
     enddo
     do t=begpatch,endpatch
        clm(t)%soilmtc_prev = soilmtc(t)
        clm(t)%h2osno_prev = clm(t)%h2osno
     enddo
  endif
  !!    print*, "AFTER restwrt ", nbeghis(1), lrestwrt

!===LDAS modification: No writing of initialization files yet
!!  if (do_inicwrite()) call inicwrt ()

!  call t_stopf('clm_output')

  return
end subroutine driver

