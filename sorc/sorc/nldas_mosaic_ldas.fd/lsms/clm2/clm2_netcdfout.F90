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
! !ROUTINE: clm2_netcdfout.F90
!
! !DESCRIPTION:  
!  LIS CLM2 data writer: Write CLM2 output in netcdf format
!
! !REVISION HISTORY:
! 02 Dec  2003: Sujay Kumar; Initial Version
! 
! !INTERFACE:
subroutine clm2_netcdfout(check,ftn)
! !USES:
  use lisdrv_module , only : lis
#if ( defined USE_NETCDF )
  use netcdf
#endif
  use drv_output_mod, only : drv_writevar_netcdf, drv_writevar_netcdf3d
  use clm_varcon, only : denh2o, denice, hvap, hsub, hfus, istwet 
  use clm_varpar, only : nlevsoi
  use clm_varmap, only : patchvec
  use clm_varder
  use clm_varctl, only : clmdrv
!EOP
  implicit none
  logical, intent(in) :: check
  real :: gtmp(lis%d%lnc,lis%d%lnr)
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
  integer :: ftn,c,t,m,i,dim1
  integer :: dimID(2),dimid1(3),days1(12),days2(12)
  data days1 /31,28,31,30,31,30,31,31,30,31,30,31/
  data days2 /31,29,31,30,31,30,31,31,30,31,30,31/
  logical :: leap
!BOC
  soilmtc=0.0
  delsoilmoist = 0.0
  delswe = 0.0
  soilmr=0.0
  soilwtc=0.0 
#if ( defined USE_NETCDF )
  if(.not.check) then 
     if((mod(lis%t%yr,4).eq.0.and.mod(lis%t%yr,100).ne.0) &    
          .or.(mod(lis%t%yr,400).eq.0)) then            
        leap = .true.                   
     else
        leap = .false.
     endif
     if(lis%t%mo .eq. lis%t%emo .and. &
          lis%t%yr .eq. lis%t%eyr) then
        dim1 = ((lis%t%eda-lis%t%da -1)*24+&
             lis%t%ehr+(24-lis%t%hr))/clmdrv%writeintc2
     else
        if(leap) then 
           dim1 = ((24-lis%t%hr)+&
                (days2(lis%t%mo)-lis%t%da)*24)/clmdrv%writeintc2
        else
           dim1 = ((24-lis%t%hr)+&
                (days1(lis%t%mo)-lis%t%da)*24)/clmdrv%writeintc2
        endif
     endif
     print*, 'dim1 ',dim1
     call check_nc(nf90_def_dim(ftn, 'land',lis%d%glbnch, dimID(1)))
     call check_nc(nf90_def_dim(ftn, 'time',dim1, dimID(2)))

     call check_nc(nf90_def_dim(ftn, 'land1',lis%d%glbnch, dimID1(1)))
     call check_nc(nf90_def_dim(ftn, 'soil1',10,dimid1(2)))
     call check_nc(nf90_def_dim(ftn, 'time1',dim1, dimID1(3)))

     call check_nc(nf90_def_var(ftn, "Swnet", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(1)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(1),"units","W/m2"))

     call check_nc(nf90_def_var(ftn, "Lwnet", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(2)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(2),"units","W/m2"))

     call check_nc(nf90_def_var(ftn, "Qle", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(3)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(3),"units","W/m2"))

     call check_nc(nf90_def_var(ftn, "Qh", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(4)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(4),"units","W/m2"))

     call check_nc(nf90_def_var(ftn, "Qg", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(5)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(5),"units","W/m2"))

     call check_nc(nf90_def_var(ftn, "Snowf", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(6)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(6),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "Rainf", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(7)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(7),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "Evap", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(8)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(8),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "Qs", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(9)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(9),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "Qsb", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(10)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(10),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "Qsm", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(11)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(11),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "DelSoilMoist", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(12)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(12),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "DelSWE", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(13)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(13),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "SnowT", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(14)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(14),"units","K"))

     call check_nc(nf90_def_var(ftn, "VegT", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(15)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(15),"units","K"))

     call check_nc(nf90_def_var(ftn, "BareSoilT", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(16)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(16),"units","K"))

     call check_nc(nf90_def_var(ftn, "AvgSurfT", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(17)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(17),"units","K"))

     call check_nc(nf90_def_var(ftn, "RadT", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(18)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(18),"units","K"))

     call check_nc(nf90_def_var(ftn, "Albedo", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(19)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(19),"units","-"))

     call check_nc(nf90_def_var(ftn, "SWE", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(20)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(20),"units","kg/m2"))

     call check_nc(nf90_def_var(ftn, "SoilTemp", nf90_float, &
          dimids = dimid1, varID = clmdrv%varid(21)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(21),"units","K"))    

     call check_nc(nf90_def_var(ftn, "SoilMoist", nf90_float, &
          dimids = dimid1, varID = clmdrv%varid(22)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(22),"units","kg/m2"))    

     call check_nc(nf90_def_var(ftn, "RootMoist", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(23)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(23),"units","kg/m2"))

     call check_nc(nf90_def_var(ftn, "SoilWet", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(24)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(24),"units","-"))

     call check_nc(nf90_def_var(ftn, "TVeg", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(25)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(25),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "ECanop", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(26)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(26),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "ESoil", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(27)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(27),"units","kg/m2s"))

     call check_nc(nf90_def_var(ftn, "CanopInt", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(28)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(28),"units","kg/m2"))

     call check_nc(nf90_def_var(ftn, "Acond", nf90_float, &
          dimids = dimID, varID = clmdrv%varid(29)))
     call check_nc(nf90_put_att(ftn,clmdrv%varid(29),"units","m/s"))

     call check_nc(nf90_enddef(ftn))
  endif
  dim1 = clmdrv%numout

!----------------------------------------------------------------------
! Net Surface Shortwave (absorbed) Radiation (W/m2)
!----------------------------------------------------------------------
  clm%totfsa=clm%totfsa/float(clm%count)
  call drv_writevar_netcdf(ftn,clm%totfsa, dim1,clmdrv%varid(1)) 
!----------------------------------------------------------------------
! Net Surface Longwave Radiation (W/m2)
!----------------------------------------------------------------------
  clm%toteflx_lwrad_net=-1.0*clm%toteflx_lwrad_net/float(clm%count)
  call drv_writevar_netcdf(ftn,clm%toteflx_lwrad_net, dim1,clmdrv%varid(2))
!----------------------------------------------------------------------
! Latent Heat Flux (W/m2)       
!----------------------------------------------------------------------
  clm%toteflx_lh_tot=clm%toteflx_lh_tot/float(clm%count)
  call drv_writevar_netcdf(ftn,clm%toteflx_lh_tot, dim1,clmdrv%varid(3))
!----------------------------------------------------------------------
! Sensible Heat Flux (W/m2)       
!----------------------------------------------------------------------
  clm%toteflx_sh_tot=clm%toteflx_sh_tot/float(clm%count)
  call drv_writevar_netcdf(ftn,clm%toteflx_sh_tot, dim1,clmdrv%varid(4))
!----------------------------------------------------------------------
! Ground Heat Flux (W/m2)       
!----------------------------------------------------------------------
  clm%toteflx_soil_grnd=clm%toteflx_soil_grnd/float(clm%count) 
  call drv_writevar_netcdf(ftn,clm%toteflx_soil_grnd, dim1,clmdrv%varid(5))
!----------------------------------------------------------------------
! General Water Balance Components (Time Averaged)
! Snowfall (kg/m2s)
!----------------------------------------------------------------------
  clm%totsnow = clm%totsnow/float(clm%count) 
  call drv_writevar_netcdf(ftn,clm%totsnow, dim1,clmdrv%varid(6))
!----------------------------------------------------------------------
! Rainfall (kg/m2s)       
!----------------------------------------------------------------------
  clm%totrain = clm%totrain/float(clm%count) 
  call drv_writevar_netcdf(ftn,clm%totrain, dim1,clmdrv%varid(7))
!----------------------------------------------------------------------
! Total Evaporation (kg/m2s)
!----------------------------------------------------------------------
  clm%totqflx_evap = clm%totqflx_evap/float(clm%count) 
  call drv_writevar_netcdf(ftn,clm%totqflx_evap, dim1,clmdrv%varid(8))
!----------------------------------------------------------------------
! Surface Runoff (kg/m2s)
!----------------------------------------------------------------------
  clm%totqflx_surf = clm%totqflx_surf/float(clm%count) 
  call drv_writevar_netcdf(ftn,clm%totqflx_surf, dim1,clmdrv%varid(9))
!----------------------------------------------------------------------
! Subsurface Runoff (kg/m2s)
!----------------------------------------------------------------------
  clm%totqflx_drain = clm%totqflx_drain/float(clm%count) 
  call drv_writevar_netcdf(ftn,clm%totqflx_drain, dim1,clmdrv%varid(10))
!----------------------------------------------------------------------
! Snowmelt (kg/m2s)
!----------------------------------------------------------------------
  snowmelt=clm%totqflx_snomelt/float(clm%count) 
  call drv_writevar_netcdf(ftn,snowmelt, dim1,clmdrv%varid(11))
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
  call drv_writevar_netcdf(ftn,delsoilmoist, dim1,clmdrv%varid(12))
!----------------------------------------------------------------------
! Change in Snow water equivalent (kg/m2)     
!----------------------------------------------------------------------
  delswe = (clm%h2osno-clm%h2osno_prev)/float(clm%count)
  call drv_writevar_netcdf(ftn,delswe, dim1,clmdrv%varid(13))
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
  call drv_writevar_netcdf(ftn,snowtemp, dim1,clmdrv%varid(14))
!----------------------------------------------------------------------
! Canopy Temperature(K)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,clm%t_veg, dim1,clmdrv%varid(15))
!----------------------------------------------------------------------
! Bare Soil Surface Temperature(K)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,clm%t_grnd, dim1,clmdrv%varid(16))
!----------------------------------------------------------------------
! Average Surface Temperature(K)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,asurft, dim1,clmdrv%varid(17))
!----------------------------------------------------------------------
! Effective Radiative Surface Temperature (K)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,clm%t_rad, dim1,clmdrv%varid(18))
!----------------------------------------------------------------------
! Surface Albedo
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,clm%surfalb, dim1,clmdrv%varid(19))
!----------------------------------------------------------------------
! Snow Water Equivalent (kg/m2)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,clm%h2osno, dim1,clmdrv%varid(20))
!----------------------------------------------------------------------  
! Subsurface State Variables
! Average Layer Soil Temperature (K)       
!----------------------------------------------------------------------
  do m=1,nlevsoi
     call drv_writevar_netcdf3d(ftn,clm%t_soisno(m),dim1,m,clmdrv%varid(22))
  enddo
!----------------------------------------------------------------------  
! Subsurface State Variables
! Average Layer Soil Moisture (kg/m2)       
!----------------------------------------------------------------------
  do m=1,nlevsoi
     do c=1,lis%d%glbnch
        tempvar(c)=soilm(c,m)
     enddo
     call drv_writevar_netcdf3d(ftn,tempvar, dim1,m,clmdrv%varid(21))
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
  call drv_writevar_netcdf(ftn,soilmr, dim1,clmdrv%varid(27))
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
  call drv_writevar_netcdf(ftn,soilwtc, dim1,clmdrv%varid(23))
!----------------------------------------------------------------------
! Evaporation Components
! Vegetation Transpiration (kg/m2s)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,cantrn, dim1,clmdrv%varid(25))
!----------------------------------------------------------------------
! Canopy Water Evaporation (kg/m2s)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,clm%totqflx_ecanop/2.501E6,&
       dim1,clmdrv%varid(24))
!----------------------------------------------------------------------
! Bare Soil Evaporation (kg/m2s)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,bare, dim1,clmdrv%varid(26))
!----------------------------------------------------------------------
! Canopy Interception (kg/m2)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,clm%canopint, dim1, clmdrv%varid(28))
!----------------------------------------------------------------------
! Aerodynamic Conductance (m/s)
!----------------------------------------------------------------------
  call drv_writevar_netcdf(ftn,clm%acond, dim1,clmdrv%varid(29))
#endif
!EOC
end subroutine clm2_netcdfout
  
