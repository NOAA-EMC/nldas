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
! !ROUTINE: clm2_binout.F90
!
! !DESCRIPTION:  
!  LIS CLM2 data writer: Write CLM2 output in binary format
!
! !REVISION HISTORY:
! 02 Dec  2003: Sujay Kumar; Initial Version
! 
! !INTERFACE:
subroutine clm2_binout(ftn)
! !USES:
  use lisdrv_module, only : lis
  use drv_output_mod, only : drv_writevar_bin
  use clm_varcon, only : denh2o, denice, hvap, hsub, hfus, istwet 
  use clm_varpar, only : nlevsoi
  use clm_varmap, only : patchvec
  use clm_varder
!EOP
  implicit none
  real :: snowmelt(lis%d%glbnch)
  real :: snowtemp(lis%d%glbnch)
  real :: totaldepth(lis%d%glbnch)
  real :: snowt(lis%d%glbnch)
  real :: asurft(lis%d%glbnch)
  real :: soilm(lis%d%glbnch,1:nlevsoi)
  real :: soilmtc(lis%d%glbnch)
  real :: delsoilmoist(lis%d%glbnch)
  real :: delswe(lis%d%glbnch)
  real :: soilmr(lis%d%glbnch)
  real :: soilwtc(lis%d%glbnch)
  real :: avgwatsat(lis%d%glbnch)
  real :: swetint(lis%d%glbnch)
  real :: swetintr(lis%d%glbnch)
  
  real :: cantrn(lis%d%glbnch)
  real :: bare(lis%d%glbnch)
  real :: snowevp(lis%d%glbnch)
  real :: potevp(lis%d%glbnch)
  real :: vmean,vstdev,vmin,vmax
  real :: tempvar(lis%d%glbnch)
  integer :: ftn,c,t,m,i
!BOC
  soilmtc=0.0
  delsoilmoist = 0.0
  delswe = 0.0
  soilmr=0.0
  soilwtc=0.0

!----------------------------------------------------------------------
! Net Surface Shortwave (absorbed) Radiation (W/m2)
!----------------------------------------------------------------------
  clm%totfsa=clm%totfsa/float(clm%count)
  call drv_writevar_bin(ftn,clm%totfsa) 
!----------------------------------------------------------------------
! Net Surface Longwave Radiation (W/m2)
!----------------------------------------------------------------------
  clm%toteflx_lwrad_net=-1.0*clm%toteflx_lwrad_net/float(clm%count)
  call drv_writevar_bin(ftn,clm%toteflx_lwrad_net)
!----------------------------------------------------------------------
! Latent Heat Flux (W/m2)       
!----------------------------------------------------------------------
  clm%toteflx_lh_tot=clm%toteflx_lh_tot/float(clm%count)
  call drv_writevar_bin(ftn,clm%toteflx_lh_tot)
!----------------------------------------------------------------------
! Sensible Heat Flux (W/m2)       
!----------------------------------------------------------------------
  clm%toteflx_sh_tot=clm%toteflx_sh_tot/float(clm%count)
  call drv_writevar_bin(ftn,clm%toteflx_sh_tot)
!----------------------------------------------------------------------
! Ground Heat Flux (W/m2)       
!----------------------------------------------------------------------
  clm%toteflx_soil_grnd=clm%toteflx_soil_grnd/float(clm%count) 
  call drv_writevar_bin(ftn,clm%toteflx_soil_grnd)
!----------------------------------------------------------------------
! General Water Balance Components (Time Averaged)
! Snowfall (kg/m2s)
!----------------------------------------------------------------------
  clm%totsnow = clm%totsnow/float(clm%count) 
  call drv_writevar_bin(ftn,clm%totsnow)
!----------------------------------------------------------------------
! Rainfall (kg/m2s)       
!----------------------------------------------------------------------
  clm%totrain = clm%totrain/float(clm%count) 
  call drv_writevar_bin(ftn,clm%totrain)
!----------------------------------------------------------------------
! Total Evaporation (kg/m2s)
!----------------------------------------------------------------------
  clm%totqflx_evap = clm%totqflx_evap/float(clm%count) 
  call drv_writevar_bin(ftn,clm%totqflx_evap)
!----------------------------------------------------------------------
! Surface Runoff (kg/m2s)
!----------------------------------------------------------------------
  clm%totqflx_surf = clm%totqflx_surf/float(clm%count) 
  call drv_writevar_bin(ftn,clm%totqflx_surf)
!----------------------------------------------------------------------
! Subsurface Runoff (kg/m2s)
!----------------------------------------------------------------------
  clm%totqflx_drain = clm%totqflx_drain/float(clm%count) 
  call drv_writevar_bin(ftn,clm%totqflx_drain)
!----------------------------------------------------------------------
! Snowmelt (kg/m2s)
!----------------------------------------------------------------------
  snowmelt=clm%totqflx_snomelt/float(clm%count) 
  call drv_writevar_bin(ftn,snowmelt)
!----------------------------------------------------------------------
! Calculation of total column soil moisture        
! Total soil moisture (liquid+ice) in each layer
!----------------------------------------------------------------------
  do m=1,nlevsoi 
     do t=1,lis%d%glbnch
        soilm(t,m)=clm(t)%h2osoi_liq(m)+clm(t)%h2osoi_ice(m)
     enddo
  enddo
  
  do m=1,nlevsoi
     do t=1,lis%d%glbnch
        soilmtc(t)=soilmtc(t)+soilm(t,m)
     enddo
  enddo
!----------------------------------------------------------------------
! Change in Soil Moisture (kg/m2)
!----------------------------------------------------------------------
  delsoilmoist = (soilmtc-clm%soilmtc_prev)/float(clm%count)
  call drv_writevar_bin(ftn,delsoilmoist)
!----------------------------------------------------------------------
! Change in Snow water equivalent (kg/m2)     
!----------------------------------------------------------------------
  delswe = (clm%h2osno-clm%h2osno_prev)/float(clm%count)
  call drv_writevar_bin(ftn,delswe)
!----------------------------------------------------------------------  
! Surface State Variables  
! Average Surface Temperature Calculation
!----------------------------------------------------------------------
  do t=1,lis%d%glbnch
     snowt(t)=0.
     if (clm(t)%itypwat/=istwet)then 
        if(clm(t)%snl < 0)then
           snowt(t)=clm(t)%t_soisno(clm(t)%snl+1)
        endif
     endif
     if(snowt(t)==0.)snowt(t)=lis%d%udef  
  enddo
!----------------------------------------------------------------------  
! AvgSurfT is the average surface temperature which depends on
! the snow temperature, bare soil temperature and canopy temperature
!----------------------------------------------------------------------
  do t=1,lis%d%glbnch
     if(snowt(t).ne.lis%d%udef)then
        asurft(t)=clm(t)%frac_sno*snowt(t)+ & 
             clm(t)%frac_veg_nosno*clm(t)%t_veg+  & 
             (1-(clm(t)%frac_sno+clm(t)%frac_veg_nosno))* & 
             clm(t)%t_grnd
     else
        asurft(t)=clm(t)%frac_veg_nosno*clm(t)%t_veg+ & 
             (1-clm(t)%frac_veg_nosno)*clm(t)%t_grnd
     endif
  enddo
  
  clm%totqflx_ecanop=clm%totqflx_ecanop/float(clm%count)
  cantrn=(clm%totqflx_tran_veg/float(clm%count))
  bare=(clm%totqflx_evap_grnd/float(clm%count))
  snowevp=(clm%totqflx_sub_snow/float(clm%count))
  potevp=lis%d%udef
!----------------------------------------------------------------------  
! Snow Temperature Calculation
!----------------------------------------------------------------------
  do t=1,lis%d%glbnch
     snowtemp(t)=0.
     if (clm(t)%itypwat/=istwet)then
        if(clm(t)%snl < 0)then
           totaldepth(t)=0.
           do i=clm(t)%snl+1,0    ! Compute total depth of snow layers
              totaldepth(t)=totaldepth(t)+clm(t)%dz(i)
           enddo
           
           do i=clm(t)%snl+1,0    ! Compute snow temperature
              snowtemp(t)=snowtemp(t)+(clm(t)%t_soisno(i)*clm(t)%dz(i))
           enddo
           snowtemp(t)=snowtemp(t)/totaldepth(t)
        endif
        if(snowtemp(t).eq.0)snowtemp(t)=lis%d%udef
     endif
  enddo
!----------------------------------------------------------------------  
! Snow Temperature (K)      
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,snowtemp)
!----------------------------------------------------------------------
! Canopy Temperature(K)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,clm%t_veg)
!----------------------------------------------------------------------
! Bare Soil Surface Temperature(K)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,clm%t_grnd)
!----------------------------------------------------------------------
! Average Surface Temperature(K)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,asurft)
!----------------------------------------------------------------------
! Effective Radiative Surface Temperature (K)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,clm%t_rad)
!----------------------------------------------------------------------
! Surface Albedo
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,clm%surfalb)
!----------------------------------------------------------------------
! Snow Water Equivalent (kg/m2)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,clm%h2osno)
!----------------------------------------------------------------------  
! Subsurface State Variables
! Average Layer Soil Temperature (K)       
!----------------------------------------------------------------------
  do m=1,nlevsoi
     call drv_writevar_bin(ftn,clm%t_soisno(m))
  enddo
!----------------------------------------------------------------------
! Subsurface State Variables
! Average Layer Soil Moisture (kg/m2)
!----------------------------------------------------------------------
  do m=1,nlevsoi
     do c=1,lis%d%glbnch
        tempvar(c)=soilm(c,m)
     enddo
     call drv_writevar_bin(ftn,tempvar)
  enddo
!----------------------------------------------------------------------
! Root Zone Soil Moisture (kg/m2)
! Calculation of root zone soil moisture 
!----------------------------------------------------------------------
  do t=1,lis%d%glbnch
     soilmr(t)=0.
     do m=1,nlevsoi
        soilmr(t)=soilmr(t)+clm(t)%rootfr(m)*clm(t)%h2osoi_liq(m)
     enddo
  enddo
  call drv_writevar_bin(ftn,soilmr)
!----------------------------------------------------------------------
! Total Soil Wetness 
! Calculation of Total column soil wetness and root zone soil wetness
! soilwtc = (vertically averaged soilm - wilting point)/
!          (vertically averaged layer porosity - wilting point)
! where average soilm is swetint, the wilting point is swetwilt,
! and avgwatsat is average porosity.
! totaldepth represents the total depth of all of the layers
!----------------------------------------------------------------------
  do t=1,lis%d%glbnch
     swetint(t)=0.
     swetintr(t)=0.
     totaldepth(t)=0.
     avgwatsat(t)=0.
     do m=1,nlevsoi
        avgwatsat(t)=avgwatsat(t)+clm(t)%dz(m)*clm(t)%watsat(m)
        totaldepth(t)=totaldepth(t)+clm(t)%dz(m)
        swetint(t)=swetint(t)+clm(t)%h2osoi_liq(m)  
        swetintr(t)=swetintr(t)+clm(t)%rootfr(m)*clm(t)%h2osoi_liq(m) 
     enddo
     avgwatsat(t)=avgwatsat(t)/totaldepth(t)
     swetint(t)=(swetint(t)/denh2o)/totaldepth(t)     
     swetintr(t)=(swetintr(t)/denh2o)/totaldepth(t) 
     soilwtc(t)=swetint(t)/avgwatsat(t)
  enddo
  call drv_writevar_bin(ftn,soilwtc)
!----------------------------------------------------------------------
! Evaporation Components
! Vegetation Transpiration (kg/m2s)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,cantrn)
!----------------------------------------------------------------------
! Canopy Water Evaporation (kg/m2s)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn, clm%totqflx_ecanop/2.501E6)
!----------------------------------------------------------------------
! Bare Soil Evaporation (kg/m2s)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,bare)
!----------------------------------------------------------------------
! Canopy Interception (kg/m2)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn, clm%canopint)
!----------------------------------------------------------------------
! Aerodynamic Conductance (m/s)
!----------------------------------------------------------------------
  call drv_writevar_bin(ftn,clm%acond)

  if(lis%o%wfor .eq. 1) then 
!----------------------------------------------------------------------
! Forcing Data Write Option On
!----------------------------------------------------------------------
! Wind (m/s)
!----------------------------------------------------------------------
     call drv_writevar_bin(ftn,sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v))
!----------------------------------------------------------------------
! Rainf (kg/m2s)
!----------------------------------------------------------------------
     call drv_writevar_bin(ftn,clm%forc_rain)
!----------------------------------------------------------------------
! Snowf (kg/m2s)
!----------------------------------------------------------------------
     call drv_writevar_bin(ftn,clm%forc_snow)
!----------------------------------------------------------------------
! Tair (K)
!----------------------------------------------------------------------
     call drv_writevar_bin(ftn,clm%forc_t)
!----------------------------------------------------------------------
! Qair (kg/kg)
!----------------------------------------------------------------------
     call drv_writevar_bin(ftn,clm%forc_q)
!----------------------------------------------------------------------
! PSurf (Pa)
!----------------------------------------------------------------------
     call drv_writevar_bin(ftn,clm%forc_pbot)
!----------------------------------------------------------------------
! SWdown (W/m2)
!----------------------------------------------------------------------
     call drv_writevar_bin(ftn,clm%forc_solad(1)*100.0/35.0)
!----------------------------------------------------------------------
! LWdown (W/m2)
!----------------------------------------------------------------------
     call drv_writevar_bin(ftn,clm%forc_lwrad)
  end if
!EOC
end subroutine clm2_binout
 
