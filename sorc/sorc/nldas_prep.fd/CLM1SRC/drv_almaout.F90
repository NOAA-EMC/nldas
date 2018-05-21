#include <misc.h>

subroutine drv_almaout (drv,tile,clm1)

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! drv_almaout.F90:
!
! DESCRIPTION:
!  This subroutine writes ALMA standard output.  Variables that are not
!  captured in the CLM are set to undefined values.
!
! REVISION HISTORY:
!  15 December 2000:  Jon Radakovich; Initial Code
!=========================================================================

! Declare Modules and data structures
  use precision
  use drv_module          ! 1-D Land Model Driver variables
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! 1-D CLM variables
  use clm1_varpar, ONLY : nlevsoi, nlevsno
  use clm1_varcon, ONLY : hfus,hsub,denice,hvap,istwet,denh2o,tfrz

  implicit none
  type (drvdec)       :: drv
  type (clm_tiledec)  :: tile(drv%nch)
  type (clm11d)       :: clm1(drv%nch)

!=== Local Variables =====================================================
  integer :: n,t,m,i,count      ! Tile space counter
  integer :: mask(drv%nch)      ! Water mask
  real(r8) &
       qg(drv%nch),             & !Ground heat flux Qg [W/m2]
       qf(drv%nch),             & !Energy of fusion Qf [W/m2]
       qv(drv%nch),             & !Energy of sublimation Qv [W/m2]
       qtau(drv%nch),           & !Momentum flux Qtau [N/m2]
       qa(drv%nch),             & !Advective energy Qa [W/m2]
       delsoilheat(drv%nch),    & !Change in soil heat storage DelSoilHeat[W/m2]
       delcoldcont(drv%nch),    & !Change in snow cold content DelColdCont[W/m2]
       qs(drv%nch),             & !Surface runoff Qs [kg/m2s]
       qrec(drv%nch),           & !Recharge Qrec [kg/m2s]
       delsoilmoist(drv%nch),   & !Change in soil moisture DelSoilMoist [kg/m2]
       delswe(drv%nch),         & !Change in snow water equivalent DelSWE [kg/m2]
       delsurfstor(drv%nch),    & !Change in surface water storage DelSurfStor [kg/m2]
       delinterc(drv%nch),      & !Change in interception storage DelIntercept [kg/m2]
       totaldepth(drv%nch),     & !Total depth of soil layers [m]
       snowt(drv%nch),          & !Snow temperature SnowT [K]
       avgsurft(drv%nch),       & !Average surface temperature AvgSurfT [K]
       albedo(drv%nch),         & !Surface albedo Albedo [-]
       surfstor(drv%nch),       & !Surface water storage SurfStor [kg/m2]
       avgwatsat(drv%nch),      & !Depth averaged volumetric soil water at saturation (porosity) [kg/m2]
       swetwilt(drv%nch),       & !Depth averaged wilting point [kg/m2]
       swetint(drv%nch),        & !Depth averaged h2osoi_liq [kg/m2]
       soilwet(drv%nch),        & !Total soil wetness [-]
       ecanop(drv%nch),         & !Interception evaporation ECanop [kg/m2s] 
       ewater(drv%nch),         & !Open water evaporation EWater [kg/m2s]
       rootmoist(drv%nch),      & !Root zone soil moisture RootMoist [kg/m2s] 
       dis(drv%nch),            & !Simulated river discharge [m3/2]
       icefrac(drv%nch),        & !Ice-covered fraction IceFrac [-]
       icet(drv%nch),           & !Sea-ice thickness IceT [m]
       testdepth,               & !Test variable used for calculation of Fdepth and Tdepth [m]
       fdepth(drv%nch),         & !Frozen soil depth Fdepth [m]
       tdepth(drv%nch),         & !Depth to soil thaw Tdepth [m]
       salbedo(drv%nch)           !Snow albedo [-]

  real(r8) :: drv_gridave                              ! Spatial Averaging Function
  real(r8) :: dz        (drv%nch,-nlevsno+1:nlevsoi)
  real(r8) :: t_soisno  (drv%nch,-nlevsno+1:nlevsoi)
  real(r8) :: h2osoi_liq(drv%nch,-nlevsno+1:nlevsoi)
  real(r8) :: h2osoi_ice(drv%nch,-nlevsno+1:nlevsoi)

!=== End Variable List ===================================================

  n=drv%nch

  do t=1,drv%nch 
#if (defined GRID_AVERAGE_NONSOIL)     
     !all points will be grid averaged, inluding lakes, wetlands and land ice
     mask(t) = 1.
#else
     ! lakes, wetlands and land-ice will not be grid averaged
     if (clm1(t)%lakpoi) then
        mask(t) = 0
     else
        mask(t) = 1
     endif
#endif
  end do

  do t = 1,drv%nch 
     do i = -nlevsno+1,nlevsoi
        dz(t,i)         = clm1(t)%dz(i)
        t_soisno(t,i)   = clm1(t)%t_soisno(i)
        h2osoi_liq(t,i) = clm1(t)%h2osoi_liq(i)
        h2osoi_ice(t,i) = clm1(t)%h2osoi_ice(i)
     enddo
  enddo

! ALMA General Energy Balance Components
  qf=hfus*clm1%qflx_snomelt
  qv=hsub*clm1%qflx_sub_snow

! Calculation for Qg & DelColdCont  
! Even if snow if present Qg is equal to the heat flux between the soil/air interface
! DelColdCont is defined as the change in internal energy of the snow pack
! over a timestep, which is zero when there is no snow.
  do t=1,drv%nch
      if(clm1(t)%snl < 0)then
         qg(t)=clm1(t)%diffusion
         delcoldcont(t)=clm1(t)%eflx_soil_grnd-clm1(t)%diffusion
      else
         qg(t)=clm1(t)%eflx_soil_grnd
         delcoldcont(t)=0.
      endif
  enddo

  qtau=sqrt((clm1%taux*clm1%taux)+(clm1%tauy*clm1%tauy))

! Qa is the heat transferred to snow cover by rain, not represented in clm1
  qa=0.

! DelSoilHeat is always zero because Qg is calculated at the soil/air interface
  delsoilheat=0.

  write(57) drv_gridave (n,mask,tile%fgrd,clm1%fsa)                 !Net shortwave radiation SWnet [W/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,(-1)*clm1%eflx_lwrad_net) !Net longwave radiation LWnet [W/m2] 
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%eflx_lh_tot)         !Latent heat flux Qle [W/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%eflx_sh_tot)         !Sensible heat flux Qh [W/m2]  
  write(57) drv_gridave (n,mask,tile%fgrd,qg)                      !Ground heat flux Qg [W/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,qf)                      !Energy of fusion Qf [W/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,qv)                      !Energy of sublimation Qv [W/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,qtau)                    !Momentum flux Qtau [N/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,qa)                      !Advective energy Qa [W/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,delsoilheat)             !Change in soil heat storage DelSoilHeat[W/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,delcoldcont)             !Change in snow cold content DelColdCont[W/m2]


! ALMA General Water Balance Components
  qs=clm1%qflx_surf+clm1%qflx_qrgwl-clm1%qflx_qirr
  qrec=-9999.0   !Recharge from river to flood plain not a capability in clm1
  do t=1,drv%nch
     delsoilmoist(t)=0.
     totaldepth(t)=0.
     do i=1,nlevsoi
        delsoilmoist(t)=delsoilmoist(t)+ &
                        (h2osoi_ice(t,i)+h2osoi_liq(t,i))-clm1(t)%h2osoi_liq_old(i)
        totaldepth(t)=totaldepth(t)+dz(t,i)
     enddo
     delsoilmoist(t)=delsoilmoist(t)
  enddo
  delswe=clm1%h2osno-clm1%h2osno_old
  delsurfstor=0.
  delinterc=clm1%h2ocan-clm1%h2ocan_old

  write(57) drv_gridave (n,mask,tile%fgrd,clm1%forc_snow)           !Snowfall rate Snowf [kg/m2s]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%forc_rain)           !Rainfall rate Rainf [kg/m2s]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%qflx_evap_tot)       !Total evapotranspiration Evap [kg/m2s]
  write(57) drv_gridave (n,mask,tile%fgrd,qs)                      !Surface runoff Qs [kg/m2s]
  write(57) drv_gridave (n,mask,tile%fgrd,qrec)                    !Recharge Qrec [kg/m2s]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%qflx_drain)          !Subsurface runoff Qsb [kg/m2s]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%qflx_snomelt)        !Rate of snowmelt Qsm [kg/m2s]  
  write(57) drv_gridave (n,mask,tile%fgrd,delsoilmoist)            !Change in soil moisture DelSoilMoist [kg/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,delswe)                  !Change in snow water equivalent DelSWE [kg/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,delsurfstor)             !Change in surface water storage DelSurfStor [kg/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,delinterc)               !Change in interception storage DelIntercept [kg/m2]

! ALMA Surface State Variables

! SnowT is the snow surface temperature, i.e. top layer t_soisno
  do t=1,drv%nch
     snowt(t)=0.
     if (clm1(t)%itypwat/=istwet)then 
        if(clm1(t)%snl < 0)then
           snowt(t)=t_soisno(t,clm1(t)%snl+1)
        endif
     endif
     if(snowt(t)==0.)snowt(t)=-9999.0  !SnowT is undefined when there is no snow
  enddo

! AvgSurfT is the average surface temperature which depends on
! the snow temperature, bare soil temperature and canopy temperature
  do t=1,drv%nch
     if(snowt(t).ne.-9999.0)then
        avgsurft(t)=clm1(t)%frac_sno*snowt(t)+clm1(t)%frac_veg_nosno*clm1(t)%t_veg+  &
                    (1-(clm1(t)%frac_sno+clm1(t)%frac_veg_nosno))*clm1(t)%t_grnd
     else
        avgsurft(t)=clm1(t)%frac_veg_nosno*clm1(t)%t_veg+ &
                    (1-clm1(t)%frac_veg_nosno)*clm1(t)%t_grnd
     endif
  enddo

  do t=1,drv%nch
     albedo(t)=clm1(t)%surfalb
  enddo

!Surface water storage not captured in CLM
  surfstor=0.

  write(57) drv_gridave (n,mask,tile%fgrd,snowt)                   !Snow temperature SnowT [K]  
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%t_veg)               !Vegetation canopy temperature VegT [K]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%t_grnd)              !Temperature of bare soil BaresoilT [K]
  write(57) drv_gridave (n,mask,tile%fgrd,avgsurft)                !Average surface temperature AvgSurfT [K]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%t_rad)               !Surface radiative temperature RadT [K]
  write(57) drv_gridave (n,mask,tile%fgrd,albedo)                  !Surface albedo Albedo [-]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%h2osno)              !Snow water equivalent SWE [kg/m2]
  write(57) drv_gridave (n,mask,tile%fgrd,surfstor)                !Surface water storage SurfStor [kg/m2]

! ALMA Subsurface State Variables

!Average layer soil moisture (liquid+ice) SoilMoist [kg/m2]
  do i=1,nlevsoi
     write(57) drv_gridave (n,mask,tile%fgrd,h2osoi_liq(:,i)+h2osoi_ice(:,i))
  enddo

!Average layer soil temperature SoilTemp [K]
  do i=1,nlevsoi
     write(57) drv_gridave (n,mask,tile%fgrd,t_soisno(:,i))
  enddo

!Average layer liquid moisture LSoilMoist [kg/m2]
  do i=1,nlevsoi
     write(57) drv_gridave (n,mask,tile%fgrd,h2osoi_liq(:,i))
  enddo

!Total soil wetness SoilWet [-]
!SoilWet = (vertically averaged SoilMoist - wilting point)/
!          (vertically averaged layer porosity - wilting point)
!where average SoilMoist is swetint, the wilting point is swetwilt,
!and avgwatsat is average porosity.
!totaldepth represents the totaldepth of all of the layers
  do t=1,drv%nch
     swetwilt(t)=0.
     swetint(t)=0.
     totaldepth(t)=0.
     avgwatsat(t)=0.
     do i=1,nlevsoi
        swetwilt(t)=swetwilt(t) + dz(t,i)*(clm1(t)%watsat(i)*((-1)*clm1(t)%smpmax/clm1(t)%sucsat(i))**(-1/clm1(t)%bsw(i)))
        avgwatsat(t)=avgwatsat(t)+dz(t,i)*clm1(t)%watsat(i)
        totaldepth(t)=totaldepth(t)+clm1(t)%dz(i)
        swetint(t)=swetint(t)+h2osoi_liq(t,i)
     enddo
     swetwilt(t)=swetwilt(t)/totaldepth(t)
     avgwatsat(t)=avgwatsat(t)/totaldepth(t)
     swetint(t)=(swetint(t)/denh2o)/totaldepth(t)     
     soilwet(t)=(swetint(t)-swetwilt(t))/(avgwatsat(t)-swetwilt(t))
  enddo
  write(57) drv_gridave (n,mask,tile%fgrd,soilwet)
    
! ALMA Evaporation Components
!Ecanop is the total evaporation from vegetation - vegetation transpiration
  ecanop=clm1%qflx_evap_veg-clm1%qflx_tran_veg
!Ewater is not represented in the clm1
  ewater=0.

!Rootmoist is the total soil moisture available for evapotranspiration
  do t=1,drv%nch
     rootmoist(t)=0.
     do i=1,nlevsoi
        rootmoist(t)=rootmoist(t)+h2osoi_liq(t,i)*clm1(t)%rootfr(i)
     enddo
  enddo

  write(57) drv_gridave (n,mask,tile%fgrd,ecanop)               !Interception evaporation ECanop [kg/m2s]  
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%qflx_tran_veg)    !Vegetation transpiration TVeg [kg/m2s]  
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%qflx_evap_grnd)   !Bare soil evaporation ESoil [kg/m2s]  
  write(57) drv_gridave (n,mask,tile%fgrd,ewater)               !Open water evaporation EWater [kg/m2s]
  write(57) drv_gridave (n,mask,tile%fgrd,rootmoist)            !Root zone soil moisture RootMoist [kg/m2s]  
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%h2ocan)           !Total canopy water storage CanopInt [kg/m2]  
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%qflx_sub_snow)    !Snow sublimation SubSnow [kg/m2s]  
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%acond)            !Aerodynamic conductance ACond [m/s]


! ALMA Streamflow 
! Dis is not captured in the clm1
  dis=-9999.0
  write(57) drv_gridave (n,mask,tile%fgrd,dis)                  !Simulated river discharge [m3/s]


! ALMA Cold Season Processes
!Icefrac and IceT are not captured in clm1 because there is no representation of sea-ice
  icefrac=-9999.0
  icet=-9999.0

!Fdepth is the frozen soil depth, which is undefined when no frozen soil is calculated
!Fdepth is calculated from the top down
  do t=1,drv%nch
     fdepth(t)=0.
     do i=1,nlevsoi
        testdepth=fdepth(t)
        if(t_soisno(t,i)<=tfrz)then
           fdepth(t)=fdepth(t)+dz(t,i)
        elseif(t_soisno(t,i)>tfrz)then
           EXIT  !If a layer is above freezing then the if statement is exited
        endif
!If the Fdepth does not change then a layer above freezing was encountered and the do loop
!must be exited
        if(testdepth.eq.fdepth(t))EXIT   
     enddo
     if(fdepth(t)==0.)fdepth(t)=-9999.0 
  enddo

!Tdepth is the thawed soil depth, which is undefined if the entire layer is thawed and zero
!when the top layer is frozen
  do t=1,drv%nch
     tdepth(t)=0.
     count=0
     do i=1,nlevsoi
        testdepth=tdepth(t)
        if(t_soisno(t,i)>tfrz)then
           count=count+1
           tdepth(t)=tdepth(t)+dz(t,i)
        elseif(t_soisno(t,i)<=tfrz)then
           EXIT  !If a layer is below freezing then the if statement is exited
        endif    
!If the Tdepth does not change then a layer below freezing was encountered and the do loop
!must be exited   
        if(testdepth.eq.tdepth(t))EXIT
     enddo
     if(count==nlevsoi)then   !Test to see if all layers are thawed
        tdepth(t)=-9999.0
     endif
  enddo

  do t=1,drv%nch
     salbedo(t)=clm1(t)%snoalb
  enddo
        
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%frac_sno)            !Snow covered fraction SnowFrac [-]
  write(57) drv_gridave (n,mask,tile%fgrd,icefrac)                 !Ice-covered fraction IceFrac [-]
  write(57) drv_gridave (n,mask,tile%fgrd,icet)                    !Sea-ice thickness IceT [m]
  write(57) drv_gridave (n,mask,tile%fgrd,fdepth)                  !Frozen soil depth Fdepth [m]
  write(57) drv_gridave (n,mask,tile%fgrd,tdepth)                  !Depth to soil thaw Tdepth [m]
  write(57) drv_gridave (n,mask,tile%fgrd,salbedo)                 !Snow albedo [-]
  write(57) drv_gridave (n,mask,tile%fgrd,clm1%snowdp)              !Depth of snow layer SnowDepth [m]

end subroutine drv_almaout
