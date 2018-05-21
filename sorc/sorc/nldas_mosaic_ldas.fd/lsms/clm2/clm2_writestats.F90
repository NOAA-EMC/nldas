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
! !ROUTINE: clm2_writestats.F90
!
! !DESCRIPTION:  
!  LIS CLM2 data writer: Write CLM2 stats
!
! !REVISION HISTORY:
! 02 Dec  2003: Sujay Kumar; Initial Version
! 
! !INTERFACE:
subroutine clm2_writestats(ftn_stats)
! !USES:
  use lisdrv_module, only : lis
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
  integer :: ftn_stats,c,t,m,i
!BOC
  soilmtc=0.0
  delsoilmoist = 0.0
  delswe = 0.0
  soilmr=0.0
  soilwtc=0.0
  

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

  delsoilmoist = (soilmtc-clm%soilmtc_prev)/float(clm%count)

  delswe = (clm%h2osno-clm%h2osno_prev)/float(clm%count)
  do t=1,lis%d%glbnch
     snowt(t)=0.
     if (clm(t)%itypwat/=istwet)then 
        if(clm(t)%snl < 0)then
           snowt(t)=clm(t)%t_soisno(clm(t)%snl+1)
        endif
     endif
     if(snowt(t)==0.)snowt(t)=lis%d%udef  
  enddo

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
  
  cantrn=(clm%totqflx_tran_veg/float(clm%count))
  bare=(clm%totqflx_evap_grnd/float(clm%count))
  snowevp=(clm%totqflx_sub_snow/float(clm%count))
  potevp=lis%d%udef

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

  do m=1,nlevsoi
     do c=1,lis%d%glbnch
        tempvar(c)=soilm(c,m)
     enddo
  enddo
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
  
  do t=1,lis%d%glbnch
     soilmr(t)=0.
     do m=1,nlevsoi
        soilmr(t)=soilmr(t)+clm(t)%rootfr(m)*clm(t)%h2osoi_liq(m)
     enddo
  enddo

  call stats(clm%totfsa,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
       vmax)
  write(ftn_stats,999)'SWnet(W/m2): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%toteflx_lwrad_net,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,999)'LWnet (W/m2): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%toteflx_lh_tot,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,999)'Qle (W/m2):    ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%toteflx_sh_tot,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,999)'Qh (W/m2):    ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%toteflx_soil_grnd,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,999)'Qg (W/m2): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%totsnow,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, &
       vmax)
  write(ftn_stats,998)'Snowf (kg/m2s): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%totrain,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
       vmax)
  write(ftn_stats,998)'Rainf (kg/m2s): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%totqflx_evap,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,998)'Evap (kg/m2s): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%totqflx_surf,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,998)'Qs (kg/m2s): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%totqflx_drain,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,998)'Qsb (kg/m2s): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%totqflx_snomelt/2.5e6,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,998)'Qsm (kg/m2s):  ', & 
       vmean,vstdev,vmin,vmax
  call stats(delsoilmoist,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,998)'DelSoilMoist (kg/m2):  ', & 
       vmean,vstdev,vmin,vmax
  call stats(delswe,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,998)'DelSWE (kg/m2): ', & 
       vmean,vstdev,vmin,vmax
  call stats(snowtemp,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999)'SnowT (K): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%t_veg,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
       vmax)
  write(ftn_stats,999)'VegT (K): ',vmean,vstdev,vmin,vmax 
  call stats(clm%t_grnd,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
       vmax)
  write(ftn_stats,999)'BaresoilT (K): ',vmean,vstdev,vmin,vmax
  call stats(asurft,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999)'AvgSurfT (K): ',vmean,vstdev,vmin,vmax
  call stats(clm%t_rad,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, &
       vmax)
  write(ftn_stats,999)'RadT (K): ',vmean,vstdev,vmin,vmax
  call stats(clm%surfalb,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
       vmax)
  write(ftn_stats,999)'Albedo (-): ',vmean,vstdev,vmin,vmax
  call stats(clm%h2ocan,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
       vmax)
  write(ftn_stats,998)'SWE (kg/m2):        ',vmean,vstdev,vmin,vmax
  
  do m=1,nlevsoi
     do c=1,lis%d%glbnch
        tempvar(c)=soilm(c,m)
     enddo
     call stats(tempvar,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)       
     write(ftn_stats,995)'SoilMoist',m,' (kg/m2): ',vmean,vstdev,vmin,vmax
  enddo
  do m=1,nlevsoi
     call stats(clm%t_soisno(m),lis%d%udef,lis%d%glbnch,&
          vmean,vstdev,vmin,vmax)       
     write(ftn_stats,995)'SoilTemp',m,'(K): ',vmean,vstdev,vmin,vmax
  enddo
  
  call stats(soilwtc,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)       
  write(ftn_stats,999)'SoilWet (-): ',vmean,vstdev,vmin,vmax

  call stats(clm%totqflx_ecanop/2.501E6,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'ECanop(kg/m2s): ',vmean,vstdev,vmin,vmax

  call stats(cantrn,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,998)'TVeg (kg/m2s):  ', & 
       vmean,vstdev,vmin,vmax
  call stats(bare,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,998)'ESoil (kg/m2s): ', & 
       vmean,vstdev,vmin,vmax
  call stats(soilmr,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,999)'RootMoist (kg/m2): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%canopint,lis%d%udef,lis%d%glbnch, & 
       vmean,vstdev,vmin,vmax)
  write(ftn_stats,999)'CanopInt (kg/m2): ', & 
       vmean,vstdev,vmin,vmax
  call stats(clm%acond,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
       vmax)
  write(ftn_stats,998)'ACond (m/s): ', & 
       vmean,vstdev,vmin,vmax
  if(lis%o%wfor.eq.1) then 
     call stats(sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v), &
          lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,999)'Wind(m/s): ', & 
          vmean,vstdev,vmin,vmax
     call stats(clm%forc_rain,lis%d%udef,lis%d%glbnch, & 
          vmean,vstdev,vmin,vmax)
     write(ftn_stats,998)'Rainf(kg/m2s): ', & 
          vmean,vstdev,vmin,vmax
     call stats(clm%forc_snow,lis%d%udef,lis%d%glbnch, & 
          vmean,vstdev,vmin,vmax)
     write(ftn_stats,998)'Snowf(kg/m2s): ', & 
          vmean,vstdev,vmin,vmax
     call stats(clm%forc_t,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
          vmax)
     write(ftn_stats,999)'Tair(K): ', & 
         vmean,vstdev,vmin,vmax
       call stats(clm%forc_q,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
            vmax)
       write(ftn_stats,999)'Qair(kg/kg): ', & 
            vmean,vstdev,vmin,vmax
       call stats(clm%forc_pbot,lis%d%udef,lis%d%glbnch,vmean, & 
            vstdev,vmin,vmax)
       write(ftn_stats,999)'PSurf(Pa): ',& 
            vmean,vstdev,vmin,vmax
       call stats(clm%forc_solad(1)*100.0/35.0,lis%d%udef,lis%d%glbnch, & 
            vmean,vstdev,vmin,vmax)
       write(ftn_stats,999)'SWdown(W/m2): ', & 
            vmean,vstdev,vmin,vmax
       call stats(clm%forc_lwrad,lis%d%udef,lis%d%glbnch,vmean, & 
            vstdev,vmin,vmax)
       write(ftn_stats,999)'LWdown(W/m2): ', & 
            vmean,vstdev,vmin,vmax
    endif
995 format (1x,a10,I1,a9,4f14.3)
999 format (1x,a15,4f14.3)
998 format (1x,a15,4e14.3)
!EOC
  end subroutine clm2_writestats
  
