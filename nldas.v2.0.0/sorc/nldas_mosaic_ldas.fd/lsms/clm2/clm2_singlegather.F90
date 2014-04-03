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
! !ROUTINE: clm2_singlegather.F90
!
! !DESCRIPTION:
!  Gather single variable for output 
!
! !REVISION HISTORY:
! Apr 2003 ; Sujay Kumar, Initial Code
!
! !INTERFACE:
subroutine clm2_singlegather(index, var)
! !USES:
  use lisdrv_module, only : lis
  use clm_varcon, ONLY : istwet,denh2o
  use clm_varpar, ONLY : nlevsoi
  use clm_varder
  use tile_spmdMod
  use clm2pardef_module
!EOP
  implicit none
  integer :: index
  real :: var(lis%d%glbnch)
  real :: snowmelt(di_array(iam))
  real :: snowt(di_array(iam))
  real :: soilmr(di_array(iam))
  real :: snowtemp(di_array(iam))
  real :: totaldepth(di_array(iam))
  real :: avgwatsat(di_array(iam))
  real :: swetint(di_array(iam))
  real :: asurft(di_array(iam))
  real :: var_temp(di_array(iam))
  integer :: t,i,m
  integer ierr
!BOC
  if ( index == 1 ) then 
     var_temp = clm%totfsa/float(clm%count)
  elseif ( index == 2 ) then 
     var_temp = -1.0*clm%toteflx_lwrad_net/float(clm%count)
  elseif ( index == 3 ) then 
     var_temp = clm%toteflx_lh_tot/float(clm%count)
  elseif ( index == 4 ) then 
     var_temp = clm%toteflx_sh_tot/float(clm%count)
  elseif ( index == 5 ) then
     var_temp = clm%toteflx_soil_grnd/float(clm%count)
  elseif ( index == 6 ) then 
     var_temp = clm%totsnow/float(clm%count)
  elseif ( index == 7 ) then
     var_temp = clm%totrain/float(clm%count)
  elseif ( index == 8 ) then 
     var_temp = clm%totqflx_evap/float(clm%count)
  elseif ( index == 9 ) then 
     var_temp = clm%totqflx_surf/float(clm%count)
  elseif ( index == 10 ) then 
     var_temp = clm%totqflx_drain/float(clm%count)
  elseif ( index == 11 ) then 
     var_temp = clm%totqflx_snomelt/float(clm%count)
  elseif ( index == 14 ) then  !SnowT
     do t=1,di_array(iam)
        snowtemp(t)=0.
        if ( clm(t)%itypwat /= istwet ) then
           if ( clm(t)%snl < 0 ) then
              totaldepth(t)=0.
              do i=clm(t)%snl+1,0    ! Compute total depth of snow layers
                 totaldepth(t)=totaldepth(t)+clm(t)%dz(i)
              enddo
              do i=clm(t)%snl+1,0    ! Compute snow temperature
                 snowtemp(t)=snowtemp(t)+(clm(t)%t_soisno(i)*clm(t)%dz(i))
              enddo
              snowtemp(t)=snowtemp(t)/totaldepth(t)
           endif
           if (snowtemp(t).eq.0)snowtemp(t)=lis%d%udef
        endif
     enddo
     var_temp = snowtemp
  elseif ( index == 15 ) then !VegT
     var_temp = clm%t_veg
  elseif ( index == 16 ) then !BareSoilT
     var_temp = clm%t_grnd
  elseif ( index == 17 ) then !AvgSurfT
     do t=1,di_array(iam)
        snowt(t) = 0.0
        if ( clm(t)%itypwat /= istwet ) then
           if ( clm(t)%snl < 0 ) then
              snowt(t) = clm(t)%t_soisno(clm(t)%snl+1)
           endif
        endif
        if ( snowt(t) == 0.0 ) then
           snowt(t) = lis%d%udef
        endif
        
        if ( snowt(t) /= lis%d%udef ) then
           asurft(t)=clm(t)%frac_sno*snowt(t)+ & 
                clm(t)%frac_veg_nosno*clm(t)%t_veg+  & 
                (1-(clm(t)%frac_sno+clm(t)%frac_veg_nosno))* & 
                clm(t)%t_grnd
        else
           asurft(t)=clm(t)%frac_veg_nosno*clm(t)%t_veg+ & 
                (1-clm(t)%frac_veg_nosno)*clm(t)%t_grnd
        endif
     enddo
     var_temp = asurft
  elseif ( index == 18 ) then !AvgSurfT
     var_temp = clm%t_rad
  elseif ( index == 19 ) then  !Albedo
     var_temp = clm%surfalb
  elseif ( index == 20 ) then  !SWE
     var_temp = clm%h2osno
  elseif ( index == 21 ) then  !SoilTemp1
     var_temp = clm%t_soisno(1)
  elseif ( index == 22 ) then  !SoilTemp2
     var_temp = clm%t_soisno(2)
  elseif ( index == 23 ) then  !SoilTemp3
     var_temp = clm%t_soisno(3)
  elseif ( index == 24 ) then  !SoilTemp4
     var_temp = clm%t_soisno(4)
  elseif ( index == 25 ) then  !SoilTemp5
     var_temp = clm%t_soisno(5)
  elseif ( index == 26 ) then  !SoilTemp6
     var_temp = clm%t_soisno(6)
  elseif ( index == 27 ) then  !SoilTemp7
     var_temp = clm%t_soisno(7)
  elseif ( index == 28 ) then  !SoilTemp8
     var_temp = clm%t_soisno(8)
  elseif ( index == 29 ) then  !SoilTemp9
     var_temp = clm%t_soisno(9)
  elseif ( index == 30 ) then  !SoilTemp10
     var_temp = clm%t_soisno(10)
  elseif ( index == 31 ) then  !SoilMoist1
     var_temp = clm%h2osoi_liq(1)+clm%h2osoi_ice(1)
  elseif ( index == 32 ) then  !SoilMoist2
     var_temp = clm%h2osoi_liq(2)+clm%h2osoi_ice(2)
  elseif ( index == 33 ) then  !SoilMoist3
     var_temp = clm%h2osoi_liq(3)+clm%h2osoi_ice(3)
  elseif ( index == 34 ) then  !SoilMoist4
     var_temp = clm%h2osoi_liq(4)+clm%h2osoi_ice(4)
  elseif ( index == 35 ) then  !SoilMoist5
     var_temp = clm%h2osoi_liq(5)+clm%h2osoi_ice(5)
  elseif ( index == 36 ) then  !SoilMoist6
     var_temp = clm%h2osoi_liq(6)+clm%h2osoi_ice(6)
  elseif ( index == 37 ) then  !SoilMoist7
     var_temp = clm%h2osoi_liq(7)+clm%h2osoi_ice(7)
  elseif ( index == 38 ) then  !SoilMoist8
     var_temp = clm%h2osoi_liq(8)+clm%h2osoi_ice(8)
  elseif ( index == 39 ) then  !SoilMoist9
     var_temp = clm%h2osoi_liq(9)+clm%h2osoi_ice(9)
  elseif ( index == 40 ) then  !SoilMoist10
     var_temp = clm%h2osoi_liq(10)+clm%h2osoi_ice(10)
  elseif ( index == 41 ) then !RootMoist
     do t=1,di_array(iam)
        soilmr(t) = 0.0
        do m=1,nlevsoi
           soilmr(t) = soilmr(t)+clm(t)%rootfr(m)*clm(t)%h2osoi_liq(m)
        enddo
     enddo
     var_temp = soilmr
  elseif ( index == 42 ) then !Soilwet
     do t=1,di_array(iam)
        swetint(t) = 0.0
        avgwatsat(t) = 0.0
        totaldepth(t) = 0.0
        do m=1,nlevsoi
           avgwatsat(t)=avgwatsat(t)+clm(t)%dz(m)*clm(t)%watsat(m)
           totaldepth(t)=totaldepth(t)+clm(t)%dz(m)
           swetint(t)=swetint(t)+clm(t)%h2osoi_liq(m)  
        enddo
        avgwatsat(t) = avgwatsat(t)/totaldepth(t)
        swetint(t) = (swetint(t)/denh2o)/totaldepth(t)
        var_temp(t) = 100*swetint(t)/avgwatsat(t)
     enddo
  elseif ( index == 43 ) then  !TVeg
     var_temp = clm%totqflx_tran_veg/float(clm%count)
  elseif ( index == 44 ) then  !ECanop
     var_temp = clm%totqflx_ecanop/float(clm%count)
     var_temp = var_temp / 2.501E6
  elseif ( index == 45 ) then  !ESoil
     var_temp = clm%totqflx_evap_grnd/float(clm%count)
  elseif ( index == 46 ) then  !Canopint
     var_temp = clm%canopint
  elseif ( index == 47 ) then !ACond
     var_temp = clm%acond
! Forcing
  elseif ( index == 48 ) then !Wind
     var_temp = sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v)
  elseif ( index == 49 ) then !Rainfall rate
     var_temp = clm%forc_rain
  elseif ( index == 50 ) then !Snowfall rate
     var_temp = clm%forc_snow
  elseif ( index == 51 ) then !Air temperature
     var_temp = clm%forc_t
  elseif ( index == 52 ) then !Specific humidity
     var_temp = clm%forc_q
  elseif ( index == 53 ) then !Surface pressure
     var_temp = clm%forc_pbot
  elseif ( index == 54 ) then !Shortwave down
     var_temp = clm%forc_solad(1)*100.0/35.0
  elseif ( index == 55 ) then !Longwave down
     var_temp = clm%forc_lwrad
  endif
#if (defined SPMD)  
  call MPI_GATHERV(var_temp(1:di_array(iam)),di_array(iam), & 
       MPI_REAL,var,di_array,displs,MPI_REAL, & 
       0,MPI_COMM_WORLD, ierr)
#else
   var = var_temp
#endif 
!EOC
end subroutine clm2_singlegather

