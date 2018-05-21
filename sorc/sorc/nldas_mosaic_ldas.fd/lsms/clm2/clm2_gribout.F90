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
! !ROUTINE: clm2_gribout.F90
!
! !DESCRIPTION:  
!  LIS CLM2 data writer: Write CLM2 output in Grib format
!
! !REVISION HISTORY:
! 02 Dec  2003: Sujay Kumar; Initial Version
! 
! !INTERFACE:
subroutine clm2_gribout(ftn)
! !USES:
  use lis_module
  use lisdrv_module, only : lis, gindex
  use drv_output_mod, only : drv_writevar_grib
  use clm_varcon, only : denh2o, denice, hvap, hsub, hfus, istwet 
  use clm_varpar, only : nlevsoi
  use clm_varmap, only : patchvec
  use clm_varctl, only : clmdrv
  use clm_varder
  use time_manager, only : tick

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
  real :: snowmeltflux(lis%d%glbnch)
  
  real :: cantrn(lis%d%glbnch)
  real :: bare(lis%d%glbnch)
  real :: snowevp(lis%d%glbnch)
  real :: potevp(lis%d%glbnch)
  real :: vmean,vstdev,vmin,vmax
  real :: tempvar(lis%d%glbnch)
  integer :: ftn,c,t,m,i,r,k
  logical*1 :: lismask(lis%d%lnc,lis%d%lnr)
  character*8 :: today, yesterday
  character*1 :: tod(8), yes(8)
  character(len=100) :: temp1
  real*8 :: dummytime
  real  :: dummygmt
  integer:: ss1,ts,mn1,hr1,da1,mo1,yr1,ts1,doy1
  integer :: kpds(25)
  real :: interval

!BOC
  soilmtc=0.0
  delsoilmoist = 0.0
  delswe = 0.0
  soilmr=0.0
  soilwtc=0.0
  
  interval = clmdrv%writeintc2
  hr1=lis%t%hr
  da1=lis%t%da
  mo1=lis%t%mo
  yr1=lis%t%yr
  mn1=lis%t%mn
  ss1=0
  ts1=-3600*24
  dummygmt=1.0
  dummytime=1.0

  write(unit=temp1,fmt='(i4,i2,i2)')yr1,mo1,da1
  read(unit=temp1,fmt='(8a1)')tod
  do i=1,8
     if(tod(i).eq.(' '))tod(i)='0'
  enddo
  today=tod(1)//tod(2)//tod(3)//tod(4)//tod(5) &
       //tod(6)//tod(7)//tod(8)

  call tick(dummytime,doy1,dummygmt,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  write(unit=temp1,fmt='(i4,i2,i2)')yr1,mo1,da1
  read(unit=temp1,fmt='(8a1)')yes
  do i=1,8
     if(yes(i).eq.(' '))yes(i)='0'
  enddo
  yesterday=yes(1)//yes(2)//yes(3)//yes(4)//yes(5) &
       //yes(6)//yes(7)//yes(8)

  do i=1,25
     kpds(i)=0
  enddo
  kpds(1)=221               !id for gsfc products
  kpds(2)=221               !id for clm2 model (change value for other models)
  kpds(4)=192               !bms flag... don't worry about this.
  kpds(12)=0                !assume output time minute always = 0
  kpds(13)=1                !forecast time unit (hours)
  kpds(17)=int((clmdrv%writeintc2*3600.0)/lis%t%ts) !number of time steps in
  !averaged/accum variables
  kpds(18)=0                !grib version -- left as 0 in ncep products
  kpds(19)=1                !version number of kpds.tbl for lisas.
  kpds(20)=0                !none missing from averages/accumulations (always4)
  kpds(23)=221              !gsfc id#
  kpds(24)=0                !does not apply to lisas output
  kpds(25)=0

  open (unit = 69, file = './src/tables/KPDS_completeclm.tbl')
  do k = 1, 42
     read(69,*)
  end do

  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        if(gindex(c,r).gt.0.5) then
           lismask(c,r)=.true.
        else
           lismask(c,r)=.false.
        endif
     enddo
  enddo

!----------------------------------------------------------------------
! Net Surface Shortwave (Absorbed) Radiation (W/m2)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%totfsa=clm%totfsa/float(clm%count)
  call drv_writevar_grib(ftn,clm%totfsa,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Net Surface Longwave Radiation (W/m2)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%toteflx_lwrad_net=-1.0*clm%toteflx_lwrad_net/float(clm%count)
  call drv_writevar_grib(ftn,clm%toteflx_lwrad_net,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Latent Heat Flux (W/m2)       
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%toteflx_lh_tot=clm%toteflx_lh_tot/float(clm%count)
  call drv_writevar_grib(ftn,clm%toteflx_lh_tot,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Sensible Heat Flux (W/m2)       
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%toteflx_sh_tot=clm%toteflx_sh_tot/float(clm%count)
  call drv_writevar_grib(ftn,clm%toteflx_sh_tot,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Ground Heat Flux (W/m2)       
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%toteflx_soil_grnd=clm%toteflx_soil_grnd/float(clm%count)
  call drv_writevar_grib(ftn,clm%toteflx_soil_grnd,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! General Water Balance Components (Time Averaged)
! Snowfall (kg/m2s)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%totsnow = clm%totsnow/float(clm%count)
  call drv_writevar_grib(ftn,clm%totsnow,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Rainfall (kg/m2s)       
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%totrain = clm%totrain/float(clm%count)
  call drv_writevar_grib(ftn,clm%totrain,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Total Evaporation (kg/m2s)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%totqflx_evap = clm%totqflx_evap/float(clm%count)
  call drv_writevar_grib(ftn,clm%totqflx_evap,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Surface Runoff (kg/m2s)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%totqflx_surf = clm%totqflx_surf/float(clm%count)
!  clm%totqflx_surf = (clm%totqflx_surf/float(clm%count))*interval*3600.0
  call drv_writevar_grib(ftn,clm%totqflx_surf,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Subsurface Runoff (kg/m2s)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  clm%totqflx_drain = clm%totqflx_drain/float(clm%count)
!  clm%totqflx_drain = (clm%totqflx_drain/float(clm%count))*interval*3600.0
  call drv_writevar_grib(ftn,clm%totqflx_drain,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Snowmelt (kg/m2s)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  snowmelt=clm%totqflx_snomelt/float(clm%count)
!  snowmelt=(clm%qflx_snomelt/float(clm%count))*interval*3600.0
  call drv_writevar_grib(ftn,snowmelt,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Snowmelt Heat flux(W/m2)
!----------------------------------------------------------------------
!  call readkpdsclm(69,kpds)
!  snowmeltflux=clm%totqflx_snomelt/float(clm%count)
!  call drv_writevar_grib(ftn,snowmeltflux,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Total soil moisture (liquid+ice) in each layer
! Calculation of total column soil moisture        
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
  call readkpdsclm(69,kpds)
  delsoilmoist = (soilmtc-clm%soilmtc_prev)/float(clm%count)
  call drv_writevar_grib(ftn,delsoilmoist,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Change in Snow water equivalent (kg/m2)     
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  delswe = (clm%h2osno-clm%h2osno_prev)/float(clm%count)
  call drv_writevar_grib(ftn,delswe,kpds,lismask,interval,today,yesterday)
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
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,snowtemp,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Canopy Temperature(K)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,clm%t_veg,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Bare Soil Surface Temperature(K)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,clm%t_grnd,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Average Surface Temperature(K)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,asurft,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Effective Radiative Surface Temperature (K)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,clm%t_rad,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Surface Albedo
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,clm%surfalb,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Snow Water Equivalent (kg/m2)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,clm%h2osno,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------  
! Subsurface State Variables
! Average Layer Soil Temperature (K)       
!----------------------------------------------------------------------
  do m=1,nlevsoi
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,clm%t_soisno(m),kpds,lismask,interval,today,yesterday)
  enddo
!----------------------------------------------------------------------
! Average Layer Soil Moisture (kg/m2)
!----------------------------------------------------------------------
  do m=1,nlevsoi
     do c=1,lis%d%glbnch
        tempvar(c)=soilm(c,m)
     enddo
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,tempvar,kpds,lismask,interval,today,yesterday)
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
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,soilmr,kpds,lismask,interval,today,yesterday)
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
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,soilwtc,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Evaporation Components
! Vegetation Transpiration (kg/m2s)
!----------------------------------------------------------------------
  cantrn=(clm%totqflx_tran_veg/float(clm%count))
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,cantrn,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Canopy Water Evaporation (kg/m2s)
!----------------------------------------------------------------------
  clm%totqflx_ecanop=clm%totqflx_ecanop/float(clm%count)
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,clm%totqflx_ecanop/2.501E6,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Potential Evaporation (kg/m2s)
!----------------------------------------------------------------------
!  potevp=lis%d%udef
!----------------------------------------------------------------------
! Bare Soil Evaporation (kg/m2s)
!----------------------------------------------------------------------
  bare=(clm%totqflx_evap_grnd/float(clm%count))
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,bare,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Canopy Interception
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,clm%canopint,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Snowpack Evaporation (kg/m2s)
!----------------------------------------------------------------------
!  snowevp=(clm%totqflx_sub_snow/float(clm%count))
!  call readkpdsclm(69,kpds)
!  call drv_writevar_grib(ftn,snowevp,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Aerodynamic Conductance (m/s)
!----------------------------------------------------------------------
  call readkpdsclm(69,kpds)
  call drv_writevar_grib(ftn,clm%acond,kpds,lismask,interval,today,yesterday)

  if(lis%o%wfor .eq. 1) then 
!----------------------------------------------------------------------
! Forcing Data Write Option On
!----------------------------------------------------------------------
! Wind (m/s)
!----------------------------------------------------------------------
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v), &
                            kpds,lismask,interval,today,yesterday)
!    call readkpdsclm(69,kpds)
!     call drv_writevar_grib(ftn,clm%forc_u,kpds,lismask,interval,today,yesterday)
!    call readkpdsclm(69,kpds)
!     call drv_writevar_grib(ftn,clm%forc_v,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Rainf (kg/m2s) - ** NOTE - CHANGE TO ACPCP-Convective Precip Field
!----------------------------------------------------------------------
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,clm%forc_rain,kpds,lismask,interval,today,yesterday)
!     clm%totrain = clm%totrain/float(clm%count) 
!     call drv_writevar_grib(ftn,clm%totrain,kpds,lismask,interval,today,yesterday)
!     call drv_writevar_grib(ftn,(clm%forc_rain)*interval*3600.0,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Snowf (kg/m2s)
!----------------------------------------------------------------------
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,clm%forc_snow,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Tair (K)
!----------------------------------------------------------------------
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,clm%forc_t,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! Qair (kg/kg)
!----------------------------------------------------------------------
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,clm%forc_q,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! PSurf (Pa)
!----------------------------------------------------------------------
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,clm%forc_pbot,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! SWdown (W/m2)
!----------------------------------------------------------------------
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,clm%forc_solad(1)*100.0/35.0,kpds,lismask,interval,today,yesterday)
!----------------------------------------------------------------------
! LWdown (W/m2)
!----------------------------------------------------------------------
     call readkpdsclm(69,kpds)
     call drv_writevar_grib(ftn,clm%forc_lwrad,kpds,lismask,interval,today,yesterday)
  end if
  close(69)

!EOC
end subroutine clm2_gribout
  
