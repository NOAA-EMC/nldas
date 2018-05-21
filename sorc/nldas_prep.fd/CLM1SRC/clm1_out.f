!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! clm_out.f: 
!
! DESCRIPTION:
!  CLM output writer.
!
! REVISION HISTORY:
! 29 Oct. 1999: Jon Radakovich; Initial code
! 27 Sep. 2000: Brian Cosgrove; Major revisions to enable CLM to 
!               output ALMA/LDAS variables
! 27 Sep. 2000: Added abitrary root zone cutoff value of .05 so 
!               that root zone is considered only those levels
!               with a rooting value >= .05
! 19 Mar  2000: Cosgrove; Changed code to allow LDAS GRIB output.  Removed
!               calls to LATS4D that had been used for grib output and 
!               added a call to griboutclm which outputs GRIB CLM data.
! 12 Apr  2002: Jambor; Modified output directory structure of HDF and
!               binary files to match that of GRIB output. 
! 21 May  2002: Cosgrove; fixed calculation of snowmelt variable
!               Previous version had been multiplying by ldas%ts
!               twice for hdf/binary output
!=========================================================================
      subroutine clm1_out (ldas, tile, grid, clm_tile, clm1)
 
!=== Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use tile_module      ! LDAS non-model-specific tile variables
      use grid_module      ! LDAS non-model-specific grid variables
      use clm1type          ! 1-D CLM variables
      use drv_tilemodule       ! 1-D CLM variables
      use clm1_varcon, ONLY : denh2o, denice, hvap, hsub, hfus, istwet 
      use clm1_varpar, ONLY : nlevsoi
      implicit none
      type (ldasdec)      ::  ldas
      type (tiledec)      ::  tile(ldas%nch)
      type (griddec)      ::  grid(ldas%nc,ldas%nr)
      type (clm_tiledec)  ::  clm_tile(ldas%nch)
      type (clm11d)        ::  clm1(ldas%nch)

!=== Local variables =====================================================
      integer :: t,c,r,m,i,j,flag,tt

!=== Temporary transfer variables 
      real :: tempvar(ldas%nch)
      real :: tempvarts(ldas%nc,ldas%nr,1:ldas%maxt)
      real :: tempvarb(ldas%nc,ldas%nr)
      real :: tempvarts3(ldas%nc,ldas%nr,1:nlevsoi)

!=== Variables used for LDAS standard output
      real :: depth,factor
      real :: snowflux(ldas%nch)
      real :: snowmelt(ldas%nch)
      real :: snowtemp(ldas%nch)
      real :: totaldepth(ldas%nch)
      real :: snowt(ldas%nch)
      real :: asurft(ldas%nch)
      real :: soilm(ldas%nch,1:nlevsoi)
      real :: soilmtc(ldas%nch)
      real :: soilmr(ldas%nch)
      real :: soilm1m(ldas%nch)
      real :: soilwtc(ldas%nch)
      real :: soilwr(ldas%nch)
      real :: swetwilt(ldas%nch)
      real :: avgwatsat(ldas%nch)
      real :: swetint(ldas%nch)
      real :: swetintr(ldas%nch)

      real :: cantrn(ldas%nch)
      real :: bare(ldas%nch)
      real :: snowevp(ldas%nch)
      real :: potevp(ldas%nch)
      real :: ccond(ldas%nch)

      real :: g_t_soisno(ldas%nc,ldas%nr,1:nlevsoi)
      real :: g_h2osoi_liq(ldas%nc,ldas%nr,1:nlevsoi)
      real :: g_h2osoi_ice(ldas%nc,ldas%nr,1:nlevsoi)
      real :: g_soilm(ldas%nc,ldas%nr,1:nlevsoi)

      real :: g_fsa(ldas%nc,ldas%nr)
      real :: g_eflx_lwrad_net(ldas%nc,ldas%nr)
      real :: g_eflx_lh_tot(ldas%nc,ldas%nr)
      real :: g_eflx_sh_tot(ldas%nc,ldas%nr)
      real :: g_eflx_soil_grnd(ldas%nc,ldas%nr)
      real :: g_snowflux(ldas%nc,ldas%nr)
      real :: g_solisbd(ldas%nc,ldas%nr)
      real :: g_forc_lwrad(ldas%nc,ldas%nr)
      real :: g_snow(ldas%nc,ldas%nr)
      real :: g_rain(ldas%nc,ldas%nr)
      real :: g_qflx_evap(ldas%nc,ldas%nr)
      real :: g_qflx_surf(ldas%nc,ldas%nr)
      real :: g_qflx_drain(ldas%nc,ldas%nr)
      real :: g_snowmelt(ldas%nc,ldas%nr)
      real :: g_snowtemp(ldas%nc,ldas%nr)
      real :: g_t_veg(ldas%nc,ldas%nr)
      real :: g_t_grnd(ldas%nc,ldas%nr)
      real :: g_asurft(ldas%nc,ldas%nr)
      real :: g_t_rad(ldas%nc,ldas%nr)
      real :: g_surfalb(ldas%nc,ldas%nr)
      real :: g_h2osno(ldas%nc,ldas%nr)
      real :: g_h2ocan(ldas%nc,ldas%nr)
      real :: g_soilmtc(ldas%nc,ldas%nr)
      real :: g_soilmr(ldas%nc,ldas%nr)
      real :: g_soilm1m(ldas%nc,ldas%nr)
      real :: g_soilwtc(ldas%nc,ldas%nr)
      real :: g_soilwr(ldas%nc,ldas%nr)
      real :: g_qflx_ecanop(ldas%nc,ldas%nr)
      real :: g_cantrn(ldas%nc,ldas%nr)
      real :: g_bare(ldas%nc,ldas%nr)
      real :: g_snowevp(ldas%nc,ldas%nr)
      real :: g_potevp(ldas%nc,ldas%nr)
      real :: g_acond(ldas%nc,ldas%nr)
      real :: g_ccond(ldas%nc,ldas%nr)
      real :: g_tlai(ldas%nc,ldas%nr)
      real :: g_snowdp(ldas%nc,ldas%nr)
      real :: g_frac_sno(ldas%nc,ldas%nr)
      real :: g_snowa(ldas%nc,ldas%nr)
      real :: g_forc_t(ldas%nc,ldas%nr)
      real :: g_forc_q(ldas%nc,ldas%nr)
      real :: g_forc_u(ldas%nc,ldas%nr)
      real :: g_forc_v(ldas%nc,ldas%nr)
      real :: g_forc_pbot(ldas%nc,ldas%nr)



!=== Variables used for statistical summary
      real vmean,vstdev,vmin,vmax

      character*80 fileng,filent,filent3,mkfyrmo,cdir,namet,nameg,namet3
      character*80 filengb
      character*1  fname(80),fbase(40),fsubst(80),fmkdir(80)
      character*1  ftime(8),fyrmodir(80),fsubsg(80),fsubst3(80)
      character*1  fcd(3),frm(3),flats(13),ftimeb(10),fsubgb(8)
      character*1 ftimec(4)

!=== Variables used for writing output in HDF format
      integer,parameter :: nvarsg=46,nvarst=45,nvarst3=9
      character*80 :: titleg,titlet,source,contact
      character*12 :: levunits
      character*80 :: vnameg(nvarsg),vtitleg(nvarsg),vunitsg(nvarsg)
      character*80 :: vnamet(nvarst),vtitlet(nvarst),vunitst(nvarst)
      character*80 :: vnamet3(nvarst3),vtitlet3(nvarst3)
      character*80 :: vunitst3(nvarst3)
      integer :: kmvarg(nvarsg),kmvart(nvarst),kmvart3(nvarst3)
      integer :: kbeg,timinc
      real :: lat(ldas%nr),lon(ldas%nc)
      real :: levsg(nlevsoi),levst(ldas%maxt)
      real :: valid_rangeg(2,nvarsg),packing_rangeg(2,nvarsg)
      real :: valid_ranget(2,nvarst),packing_ranget(2,nvarst)
      real :: valid_ranget3(2,nvarst3),packing_ranget3(2,nvarst3)
      integer :: prec,kount

      data titleg /"CLM Grid-space Output based on LDAS Forcing"/
      data titlet /"CLM Tile-space Output based on LDAS Forcing"/
      data source /"NASA GSFC/Data Assimilation Office"/
      data contact /"houser@dao.gsfc.nasa.gov"/
      data levunits /"sigma_level"/

      data vtitlet /"Vegetation Type of Tile",
     &  "Fraction of Grid Covered by Tile",
     &  "Net Surface Shortwave Radiation (W/m2)",
     &  "Net Surface Longwave Radiation (W/m2)",
     &  "Latent Heat Flux (W/m2)",
     &  "Sensible Heat Flux (W/m2)",
     &  "Ground Heat Flux (W/m2)",
     &  "Snow Phase Change Heat Flux (W/m2)",
     &  "Downward Surface Shortwave Radiation (W/m2)",
     &  "Downward Surface Longwave Radiation (W/m2)",
     &  "Snowfall (kg/m2)",
     &  "Rainfall (kg/m2)",
     &  "Total Evaporation (kg/m2)",
     &  "Surface Runoff (kg/m2)",
     &  "Subsurface Runoff (kg/m2)",
     &  "Snowmelt (kg/m2)",
     &  "Snow Temperature (K)",
     &  "Canopy Temperature (K)",
     &  "Bare Soil Surface Temperature (K)",
     &  "Average Surface Temperature (K)",
     &  "Effective Radiative Surface Temperature (K)",
     &  "Surface Albedo, All Wavelengths (%)",
     &  "Snowpack Water Equivalent (kg/m2)",
     &  "Plant Canopy Surface Water Storage (kg/m2)",
     &  "Total Column Soil Moisture (kg/m2)",
     &  "Root Zone Soil Moisture (kg/m2)",
     &  "Top 1-meter Soil Moisture (kg/m2)",
     &  "Total Soil Column Wetness (%)",
     &  "Root Zone Wetness (%)",
     &  "Canopy Surface Water Evaporation (W/m2)",
     &  "Canopy Transpiration (W/m2)",
     &  "Bare Soil Evaporation (W/m2)",
     &  "Snow Evaporation (W/m2)",
     &  "Potential Evaporation (W/m2)",
     &  "Aerodynamic Conductance (m/s)",
     &  "Canopy Conductance (m/s)",
     &  "Leaf Area Index",
     &  "Snow Depth (m)",
     &  "Snow Cover (%)",
     &  "Snow Albedo (%)",
     &  "Two Meter Temperature (K)",
     &  "Two Meter Humidity (kg/kg)",
     &  "Ten Meter U Wind (m/s)",
     &  "Ten Meter V Wind (m/s)",
     &  "Surface Pressure (mb)" /


      DATA VTITLEG /
     &  "Net Surface Shortwave Radiation (W/m2)",
     &  "Net Surface Longwave Radiation (W/m2)",
     &  "Latent Heat Flux (W/m2)",
     &  "Sensible Heat Flux (W/m2)",
     &  "Ground Heat Flux (W/m2)",
     &  "Snow Phase Change Heat Flux (W/m2)",
     &  "Downward Surface Shortwave Radiation (W/m2)",
     &  "Downward Surface Longwave Radiation (W/m2)",
     &  "Snowfall (kg/m2)",
     &  "Rainfall (kg/m2)",
     &  "Total Evaporation (kg/m2)",
     &  "Surface Runoff (kg/m2)",
     &  "Subsurface Runoff (kg/m2)",
     &  "Snowmelt (kg/m2)",
     &  "Snow Temperature (K)",
     &  "Canopy Temperature (K)",
     &  "Bare Soil Surface Temperature (K)",
     &  "Average Surface Temperature (K)",
     &  "Effective Radiative Surface Temperature (K)",
     &  "Surface Albedo, All Wavelengths (%)",
     &  "Snowpack Water Equivalent (kg/m2)",
     &  "Plant Canopy Surface Water Storage (kg/m2)",
     &  "Soil Temperature (K)",
     &  "Total Column Soil Moisture (kg/m2)",
     &  "Root Zone Soil Moisture (kg/m2)",
     &  "Top 1-meter Soil Moisture (kg/m2)",
     &  "Soil Moisture (kg/m2)",
     &  "Liquid Soil Moisture (kg/m2)",
     &  "Total Soil Column Wetness (%)",
     &  "Root Zone Wetness (%)",
     &  "Canopy Surface Water Evaporation (W/m2)",
     &  "Canopy Transpiration (W/m2)",
     &  "Bare Soil Evaporation (W/m2)",
     &  "Snow Evaporation (W/m2)",
     &  "Potential Evaporation (W/m2)",
     &  "Aerodynamic Conductance (m/s)",
     &  "Canopy Conductance (m/s)",
     &  "Leaf Area Index",
     &  "Snow Depth (m)",
     &  "Snow Cover (%)",
     &  "Snow Albedo (%)",
     &  "Two Meter Temperature (K)",
     &  "Two Meter Humidity (kg/kg)",
     &  "Ten Meter U Wind (m/s)",
     &  "Ten Meter V Wind (m/s)",
     &  "Surface Pressure (mb)" /

      data vtitlet3 /"Soil Layer Temperature (1)",
     & "Soil Layer Temperature (2)",
     & "Soil Layer Temperature (3)",
     & "Soil Layer Moisture (1)",
     & "Soil Layer Moisture (2)",
     & "Soil Layer Moisture (3)",
     & "Liquid in Soil Layer (1)",
     & "Liquid in Soil Layer (2)",
     & "Liquid in Soil Layer (3)" /
                  
      data vunitsg /
     &  "watts per meter squared",
     &  "watts per meter squared",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "K",
     &  "K",
     &  "K",
     &  "K",
     &  "K",
     &  "%",
     &  "kg/m2",
     &  "kg/m2",
     &  "K",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "%",
     &  "%",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "m/s",
     &  "m/s",
     &  "Unitless",
     &  "m",
     &  "%",
     &  "%",
     &  "K",
     &  "kg/kg",
     &  "m/s",
     &  "m/s",
     &  "mb" /

      data vunitst /
     &  "Type",
     &  "Fraction",
     &  "watts per meter squared",
     &  "watts per meter squared",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "K",
     &  "K",
     &  "K",
     &  "K",
     &  "K",
     &  "%",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "%",
     &  "%",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "m/s",
     &  "m/s",
     &  "Unitless",
     &  "m",
     &  "%",
     &  "%",
     &  "K",
     &  "kg/kg",
     &  "m/s",
     &  "m/s",
     &  "mb" /


      data vunitst3 /"Kelvin","Kelvin","Kelvin",
     & "kg/m2","kg/m2","kg/m2","kg/m2",
     & "kg/m2","kg/m2"/

      data vnameg / 
     &  "nswrs",
     &  "nlwrs",
     &  "lhtfl",
     &  "shtfl",
     &  "gflux",
     &  "snohf",
     &  "dswrf",
     &  "dlwrf",
     &  "asnow",
     &  "arain",
     &  "evp",
     &  "ssrun",
     &  "bgrun",
     &  "snom",
     &  "snowt",
     &  "vegt",
     &  "baret",
     &  "avsft",
     &  "radt",
     &  "albdo",
     &  "weasd",
     &  "cwat",
     &  "soilt",
     &  "soilmc",
     &  "soilmr",
     &  "soilmt1",
     &  "soilm",
     &  "lsoil",
     &  "mstavc",
     &  "mstavr",
     &  "evcw",
     &  "trans",
     &  "evbs",
     &  "sbsno",
     &  "pevpr",
     &  "acond",
     &  "ccond",
     &  "lai",
     &  "snod",
     &  "snoc",
     &  "salbd",
     &  "tmp2m",
     &  "humid",
     &  "uwind",
     &  "vwind",
     &  "sfcprs"/


      data vnamet /
     &  "vegt",
     &  "fgrd",
     &  "nswrs",
     &  "nlwrs",
     &  "lhtfl",
     &  "shtfl",
     &  "gflux",
     &  "snohf",
     &  "dswrf",
     &  "dlwrf",
     &  "asnow",
     &  "arain",
     &  "evp",
     &  "ssrun",
     &  "bgrun",
     &  "snom",
     &  "snowt",
     &  "vegt",
     &  "baret",
     &  "avsft",
     &  "radt",
     &  "albdo",
     &  "weasd",
     &  "cwat", 
     &  "soilmc",
     &  "soilmr",
     &  "soilmt1",
     &  "mstavc",
     &  "mstavr",
     &  "evcw",
     &  "trans",
     &  "evbs",
     &  "sbsno",
     &  "pevpr",
     &  "acond",
     &  "ccond",
     &  "lai",
     &  "snod",
     &  "snoc",
     &  "salbd",
     &  "tmp2m",
     &  "humid",
     &  "uwind",
     &  "vwind",
     &  "sfcprs"/


      data vnamet3 /"soilt1","soilt2","soilt3","soilm1","soilm2",
     & "soilm3","lsoil1","lsoil2","lsoil3" /

      data kmvarg /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,
     &  0,0,10,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/

      character*40 file
      character*80 name
    
!=== End Variable List =========================================================

!=== Total arrays hold a running total of each output variable for time
!=== averaging, between output writes
      clm1%count=clm1%count+1
      clm1%totfsa=clm1%totfsa+clm1%fsa
      clm1%toteflx_lwrad_net=clm1%toteflx_lwrad_net+  
     &            clm1%eflx_lwrad_net
      clm1%toteflx_lh_tot=clm1%toteflx_lh_tot+clm1%eflx_lh_tot
      clm1%toteflx_sh_tot=clm1%toteflx_sh_tot+clm1%eflx_sh_tot
      do tt=1,ldas%nch
      clm1(tt)%toteflx_soil_grnd=clm1(tt)%toteflx_soil_grnd+
     &  clm1(tt)%eflx_soil_grnd
	enddo
        clm1%totqflx_snomelt=clm1%totqflx_snomelt+
     &    clm1%qflx_snomelt*ldas%ts
      clm1%totsolisbd=clm1%totsolisbd+
     & (clm1%forc_solad(1)+clm1%forc_solad(2)+
     &  clm1%forc_solai(1)+clm1%forc_solai(2))
      clm1%totforc_lwrad=clm1%totforc_lwrad+clm1%forc_lwrad

      clm1%totsnow=clm1%totsnow+(clm1%forc_snow*ldas%ts)
      clm1%totrain=clm1%totrain+(clm1%forc_rain*ldas%ts)
      clm1%totqflx_evap=clm1%totqflx_evap+(clm1%qflx_evap_tot*ldas%ts)

      clm1%totqflx_surf=clm1%totqflx_surf+
     &     (ldas%ts*(clm1%qflx_surf+clm1%qflx_qrgwl-clm1%qflx_qirr))
      clm1%totqflx_drain=clm1%totqflx_drain+
     &    (ldas%ts*(clm1%qflx_drain))

      clm1%totqflx_ecanop=clm1%totqflx_ecanop+          
     &  (-hvap*(clm1%qflx_evap_veg-clm1%qflx_tran_veg))
      clm1%totqflx_tran_veg=clm1%totqflx_tran_veg+clm1%qflx_tran_veg*
     &                   (-hvap)
      clm1%totqflx_evap_grnd=clm1%totqflx_evap_grnd+clm1%qflx_evap_grnd*
     &                    (-hvap)
      clm1%totqflx_sub_snow=clm1%totqflx_sub_snow+clm1%qflx_sub_snow*
     &              (-hsub)


!=== Test to see if output writing interval has been reached
      if(mod(ldas%gmt,ldas%writeintc1).eq.0)then
       ldas%numoutc1=ldas%numoutc1+1    !Counts number of output times

!=== Initialize grid space average arrays 
       do m=1,nlevsoi
        do c=1,ldas%nc
         do r=1,ldas%nr
          g_t_soisno(c,r,m)=0.
          g_h2osoi_liq(c,r,m)=0.
          g_h2osoi_ice(c,r,m)=0.
	  g_soilm(c,r,m)=0.
         enddo
        enddo
       enddo 


       	g_fsa=0.
	g_eflx_lwrad_net=0.
	g_eflx_lh_tot=0.
	g_eflx_sh_tot=0.
	g_eflx_soil_grnd=0.
	g_snowflux=0.
	g_solisbd=0.
	g_forc_lwrad=0.
	g_snow=0.
	g_rain=0.
	g_qflx_evap=0.
	g_qflx_surf=0.
	g_qflx_drain=0.
	g_snowmelt=0.
	g_snowtemp=0.
	g_t_veg=0.
	g_t_grnd=0.
        g_asurft=0.
	g_t_rad=0.
        g_surfalb=0.
        g_h2osno=0.
	g_h2ocan=0.
	g_soilmtc=0.
	g_soilmr=0.
	g_soilm1m=0.
	g_soilwtc=0.
	g_soilwr=0.
	g_qflx_ecanop=0.
	g_cantrn=0.
	g_bare=0.
	g_snowevp=0.
	g_potevp=0.
	g_acond=0.
	g_ccond=0.
	g_tlai=0.
	g_snowdp=0.
	g_frac_sno=0.
	g_snowa=0.
	g_forc_t=0.
	g_forc_q=0.
	g_forc_u=0.
	g_forc_v=0.
	g_forc_pbot=0.

!=== Perform time averaging for the total arrays where necessary
!=== and find grid averages

!=== LDAS Energy Balance Components (time averaged)
       clm1%totfsa=(-1)*clm1%totfsa/float(clm1%count)
       call t2gr(clm1%totfsa,g_fsa,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       clm1%toteflx_lwrad_net=clm1%toteflx_lwrad_net/float(clm1%count)
       call t2gr(clm1%toteflx_lwrad_net,g_eflx_lwrad_net,ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row)
       clm1%toteflx_lh_tot=clm1%toteflx_lh_tot/float(clm1%count)
       call t2gr(clm1%toteflx_lh_tot,g_eflx_lh_tot,ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row)
	clm1%toteflx_sh_tot=clm1%toteflx_sh_tot/float(clm1%count)
       call t2gr(clm1%toteflx_sh_tot,g_eflx_sh_tot,ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row)
	clm1%toteflx_soil_grnd=
     &  (-1)*clm1%toteflx_soil_grnd/float(clm1%count)
       call t2gr(clm1%toteflx_soil_grnd,g_eflx_soil_grnd,ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row)
       snowflux=(clm1%totqflx_snomelt/float(clm1%count))*(-hfus)
       call t2gr(snowflux,g_snowflux,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
	clm1%totsolisbd=clm1%totsolisbd/float(clm1%count)
       call t2gr(clm1%totsolisbd,g_solisbd,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
	clm1%totforc_lwrad=clm1%totforc_lwrad/float(clm1%count)
       call t2gr(clm1%totforc_lwrad,g_forc_lwrad,ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row)

!=== LDAS Water Balance Components (accumulated)
       call t2gr(clm1%totsnow,g_snow,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%totrain,g_rain,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%totqflx_evap,g_qflx_evap,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%totqflx_surf,g_qflx_surf,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%totqflx_drain,g_qflx_drain,ldas%nc,ldas%nr,
     &  ldas%nch,tile%fgrd,tile%col,tile%row)
       snowmelt=clm1%totqflx_snomelt
       call t2gr(snowmelt,g_snowmelt,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

!=== LDAS Surface State Variables (Instantaneous)
!Snow Temperature Calculation
       do t=1,ldas%nch
        snowtemp(t)=0.
         if (clm1(t)%itypwat/=istwet)then
          if(clm1(t)%snl < 0)then
           totaldepth(t)=0.
           do i=clm1(t)%snl+1,0    ! Compute total depth of snow layers
            totaldepth(t)=totaldepth(t)+clm1(t)%dz(i)
           enddo

           do i=clm1(t)%snl+1,0    ! Compute snow temperature
            snowtemp(t)=snowtemp(t)+(clm1(t)%t_soisno(i)*clm1(t)%dz(i))
           enddo
           snowtemp(t)=snowtemp(t)/totaldepth(t)
         endif
         if(snowtemp(t).eq.0)snowtemp(t)=ldas%udef
        endif
       enddo

       call t2gr(snowtemp,g_snowtemp,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%t_veg,g_t_veg,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%t_grnd,g_t_grnd,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

!=== Average Surface Temperature Calculation
       do t=1,ldas%nch
! SnowT is the snow surface temperature, i.e. top layer t_soisno  
        snowt(t)=0.
        if (clm1(t)%itypwat/=istwet)then 
         if(clm1(t)%snl < 0)then
          snowt(t)=clm1(t)%t_soisno(clm1(t)%snl+1)
         endif
        endif
       if(snowt(t)==0.)snowt(t)=ldas%udef  !SnowT is undefined when there is no snow
       enddo

! AvgSurfT is the average surface temperature which depends on
! the snow temperature, bare soil temperature and canopy temperature
       do t=1,ldas%nch
        if(snowt(t).ne.ldas%udef)then
         asurft(t)=clm1(t)%frac_sno*snowt(t)+
     &             clm1(t)%frac_veg_nosno*clm1(t)%t_veg+ 
     &             (1-(clm1(t)%frac_sno+clm1(t)%frac_veg_nosno))*
     &             clm1(t)%t_grnd
        else
         asurft(t)=clm1(t)%frac_veg_nosno*clm1(t)%t_veg+
     &             (1-clm1(t)%frac_veg_nosno)*clm1(t)%t_grnd
        endif
       enddo
       call t2gr(asurft,g_asurft,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%t_rad,g_t_rad,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%surfalb,g_surfalb,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%h2osno,g_h2osno,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%h2ocan,g_h2ocan,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

!=== LDAS Subsurface State Variables
       do m=1,nlevsoi
        do t=1,ldas%nch
         tempvar(t)=clm1(t)%t_soisno(m)
        enddo
        call t2gr(tempvar,g_t_soisno(:,:,m),ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row) 
       enddo

       do m=1,nlevsoi
        do t=1,ldas%nch
         tempvar(t)=(clm1(t)%h2osoi_liq(m))
        enddo
        call t2gr(tempvar,g_h2osoi_liq(:,:,m),ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row)
       enddo

       do m=1,nlevsoi
        do t=1,ldas%nch
         tempvar(t)=(clm1(t)%h2osoi_ice(m))
        enddo
        call t2gr(tempvar,g_h2osoi_ice(:,:,m),ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row)
       enddo
   
! Total soil moisture (liquid+ice) in each layer
       do m=1,nlevsoi 
        do t=1,ldas%nch
         soilm(t,m)=clm1(t)%h2osoi_liq(m)+clm1(t)%h2osoi_ice(m)
        enddo
       enddo

       do m=1,nlevsoi
        do t=1,ldas%nch
         tempvar(t)=soilm(t,m)
        enddo
        call t2gr(tempvar,g_soilm(:,:,m),ldas%nc,
     &  ldas%nr,ldas%nch,tile%fgrd,tile%col,tile%row)
       enddo

       soilmtc=0.0
       soilmr=0.0
       soilm1m=0.0
       soilwtc=0.0
       soilwr=0.0

! Calculation of total column soil moisture        
       do m=1,nlevsoi
        do t=1,ldas%nch
         soilmtc(t)=soilmtc(t)+soilm(t,m)
        enddo
       enddo
       call t2gr(soilmtc,g_soilmtc,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

! Calculation of root zone soil moisture 
      do t=1,ldas%nch
        soilmr(t)=0.
        do m=1,nlevsoi
          soilmr(t)=soilmr(t)+clm1(t)%rootfr(m)*clm1(t)%h2osoi_liq(m)
        enddo
       enddo
       call t2gr(soilmr,g_soilmr,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

! Calculation of top 1-m soil moisture
       do t=1,ldas%nch
        depth=0.0
        flag=0
        do i=1,nlevsoi
          if (flag.ne.1) then
            soilm1m(t)=0.0
            depth=depth+clm1(t)%dz(i)
            if (depth.ge.1.0) then
              factor=(1.0-((depth-1.0)/clm1(t)%dz(i)))
              if (depth.le.1.00001) then
                do j=1,i
                  soilm1m(t)=soilm1m(t)+soilm(t,j)
                enddo
                flag=1
              else
                do j=1,i-1
                  soilm1m(t)=soilm1m(t)+soilm(t,j)
                enddo
                soilm1m(t)=soilm1m(t)+(factor*soilm(t,i))
                flag=1
              endif
            endif
          endif
        enddo
        if (flag.eq.0) then
          do i=1,nlevsoi
            soilm1m(t)=soilm1m(t)+soilm(t,i)
          enddo
          soilm1m(t)=soilm1m(t)/depth
        endif
       enddo
       call t2gr(soilm1m,g_soilm1m,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

! Calculation of Total column soil wetness and root zone soil wetness
!soilwtc = (vertically averaged soilm - wilting point)/
!          (vertically averaged layer porosity - wilting point)
!where average soilm is swetint, the wilting point is swetwilt,
!and avgwatsat is average porosity.
!totaldepth represents the total depth of all of the layers
       do t=1,ldas%nch
        swetwilt(t)=0.
        swetint(t)=0.
        swetintr(t)=0.
        totaldepth(t)=0.
        avgwatsat(t)=0.
        do m=1,nlevsoi
         swetwilt(t)=swetwilt(t) + clm1(t)%dz(m)*(clm1(t)%watsat(m)*
     &    ((-1)*clm1(t)%smpmax/clm1(t)%sucsat(m))**(-1/clm1(t)%bsw(m)))
         avgwatsat(t)=avgwatsat(t)+clm1(t)%dz(m)*clm1(t)%watsat(m)
         totaldepth(t)=totaldepth(t)+clm1(t)%dz(m)
! Total column soil moisture
         swetint(t)=swetint(t)+clm1(t)%h2osoi_liq(m)  
! Root zone soil moisture
         swetintr(t)=swetintr(t)+clm1(t)%rootfr(m)*clm1(t)%h2osoi_liq(m) 
        enddo
        swetwilt(t)=swetwilt(t)/totaldepth(t)
        avgwatsat(t)=avgwatsat(t)/totaldepth(t)
        swetint(t)=(swetint(t)/denh2o)/totaldepth(t)     
        swetintr(t)=(swetintr(t)/denh2o)/totaldepth(t) 
        soilwtc(t)=100.*swetint(t)/avgwatsat(t)
        soilwr(t)=100.*swetintr(t)/avgwatsat(t)
       enddo

       call t2gr(soilwtc,g_soilwtc,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

       call t2gr(soilwr,g_soilwr,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

!=== LDAS Evaporation Components
       clm1%totqflx_ecanop=clm1%totqflx_ecanop/float(clm1%count)
       call t2gr(clm1%totqflx_ecanop,g_qflx_ecanop,ldas%nc,ldas%nr,
     &  ldas%nch,tile%fgrd,tile%col,tile%row)
       cantrn=(clm1%totqflx_tran_veg/float(clm1%count))*(-hvap)
       call t2gr(cantrn,g_cantrn,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       bare=(clm1%totqflx_evap_grnd/float(clm1%count))*(-hvap)
       call t2gr(bare,g_bare,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       snowevp=(clm1%totqflx_sub_snow/float(clm1%count))*(-hsub)
       call t2gr(snowevp,g_snowevp,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       potevp=ldas%udef
       call t2gr(potevp,g_potevp,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%acond,g_acond,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       ccond=ldas%udef
       call t2gr(ccond,g_ccond,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%tlai,g_tlai,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

!=== LDAS Cold Season Processes
       call t2gr(clm1%snowdp,g_snowdp,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%frac_sno,g_frac_sno,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%snoalb,g_snowa,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)

!=== LDAS Forcing Components
       call t2gr(clm1%forc_t,g_forc_t,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%forc_q,g_forc_q,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%forc_u,g_forc_u,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%forc_v,g_forc_v,ldas%nc,ldas%nr,ldas%nch,
     &  tile%fgrd,tile%col,tile%row)
       call t2gr(clm1%forc_pbot/100.,g_forc_pbot,ldas%nc,ldas%nr,
     &  ldas%nch,tile%fgrd,tile%col,tile%row)

!=== Generate directory structure and file names for clm1 Output 
 91    FORMAT(A7,I3,A1)
 92    FORMAT(80A1)
 93    FORMAT(A80)
 94    FORMAT(I4,I2,I2)
 95    FORMAT(8A1)
 96    FORMAT(A40)
 97    FORMAT(A4,I3,A4)
 98    FORMAT(A4,I3,A5,I4,A1,I4,I2,I2)
100    FORMAT(A9)
101    FORMAT(A8)
102    FORMAT(A3)
103    FORMAT(A80,A12,A26,A4,A26,A35)
104    FORMAT(A183)
105    FORMAT(A3,A80)
106    FORMAT(A83)
107    FORMAT(40A1)
108    FORMAT(A13)
109    FORMAT(I4,I2,I2,I2)
110    FORMAT(10A1)
111    FORMAT(I4)

       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
       WRITE(90,94,REC=1)LDAS%YR,LDAS%MO,LDAS%DA
       READ(90,95,REC=1)FTIME
       DO I=1,8
        IF(FTIME(I).EQ.(' '))FTIME(I)='0'
       ENDDO

       WRITE(90,111,REC=1)LDAS%YR
       READ(90,95,REC=1)FTIMEC
       DO I=1,4
        IF(FTIMEC(I).EQ.(' '))FTIMEC(I)='0'
       ENDDO

       WRITE(90,91,REC=1)'/LDAS.E',LDAS%EXPCODE,'.'
       READ(90,92,REC=1) (FNAME(I),I=1,11)
       DO I=1,11
        IF(FNAME(I).EQ.(' '))FNAME(I)='0'
       ENDDO

       WRITE(90,96,REC=1) LDAS%ODIR
       READ(90,107,REC=1) (FBASE(I),I=1,40)
       C=0
       DO I=1,40
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO

       WRITE(90,98,REC=1)'/EXP',LDAS%EXPCODE,'/CLM/',
     &  LDAS%YR,'/',LDAS%YR,LDAS%MO,LDAS%DA
       READ(90,92,REC=1) (FYRMODIR(I),I=1,25)
       DO I=1,25
        IF(FYRMODIR(I).EQ.(' '))FYRMODIR(I)='0'
       ENDDO

       WRITE(90,100,REC=1)'mkdir -p '
       READ(90,92,REC=1)(FMKDIR(I),I=1,9)

       WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &  (FYRMODIR(I),I=1,25)
       READ(90,93,REC=1)MKFYRMO

       CLOSE(90)

!== Make the directories for the clm1 output files              
       call system(mkfyrmo)

!=== Generate file name for binary output
       if(ldas%wbin.eq.1)then
        open(90,file='temp',form='formatted',access='direct',recl=80)
        write(90,109,rec=1)ldas%yr,ldas%mo,ldas%da,ldas%hr
        read(90,110,rec=1)ftimeb
        do i=1,10
         if(ftimeb(i).eq.(' '))ftimeb(i)='0'
        enddo
        write(90,101,rec=1)'.clm1gbin'
        read(90,92,rec=1) (fsubgb(i),i=1,8)

        write(90,92,rec=1)(fbase(i),i=1,c),(fyrmodir(i),i=1,25),
     &       (fname(i),i=1,11),(ftimeb(i),i=1,10),(fsubgb(i),i=1,8 ) 
        read(90,93,rec=1)filengb
        close(90)

       endif

!=== Open HDF daily output file     
       if((ldas%numoutc1.eq.1.and.ldas%whdf.eq.1))then
        open(90,file='temp',form='formatted',access='direct',recl=80)
        write(90,101,rec=1)'.clmgrid'
        read(90,92,rec=1) (fsubsg(i),i=1,8)

        write(90,92,rec=1)(fbase(i),i=1,c),(fyrmodir(i),i=1,25), 
     &       (fname(i),i=1,11),(ftime(i),i=1,8),(fsubsg(i),i=1,8 )
        read(90,93,rec=1)fileng

        if(ldas%wtil.eq.1)then
         write(90,101,rec=1)'.clmtile'
         read(90,92,rec=1) (fsubst(i),i=1,8)
         write(90,100,rec=1)'.clmtile3'
         read(90,92,rec=1) (fsubst3(i),i=1,9)

         write(90,92,rec=1)(fbase(i),i=1,c), (fname(i),i=1,29),
     &                     (ftime(i),i=1,8),(fsubst(i),i=1,8 ) 
         read(90,93,rec=1)filent

         write(90,92,rec=1)(fbase(i),i=1,c), (fname(i),i=1,29),
     &                     (ftime(i),i=1,8),(fsubst3(i),i=1,9 ) 
         read(90,93,rec=1)filent3
        endif
        close(90)
        close(89)

!=== clm1 Output File Name Generation Complete

!=== Set up HDF Input Parameters ===============================================
        prec = 0     !Data precision (0=32 bit, 1=64 bit) 
        timinc= 10000*int(ldas%writeintc1)

!=== Set up HDF grid specific input parameters
        do c=1,2
         do r=1,nvarsg 
          valid_rangeg(c,r)=ldas%udef      !Do not use packing algorithm
          packing_rangeg(c,r)=ldas%udef
         enddo
        enddo

        do c=1,ldas%nc
         do r=1,ldas%nr
          lat(r)=grid(c,r)%lat
          lon(c)=grid(c,r)%lon
         enddo
        enddo

        do m=1,nlevsoi
         levsg(m)=clm1(1)%dz(m)
        enddo

!=== Create HDF file for grid space output===========================================

        if(ldas%wtil.eq.1)then 
!=== Set up HDF 3D tile specific parameters
         do c=1,2
          do r=1,nvarst3 
           valid_ranget3(c,r)=ldas%udef      !Do not use packing algorithm
           packing_ranget3(c,r)=ldas%udef
          enddo
         enddo

         do m=1,nvarst3
          kmvart3(m)=nlevsoi
         enddo

!=== Set up HDF tile specific input parameters
         do c=1,2
          do r=1,nvarst 
           valid_ranget(c,r)=ldas%udef      !Do not use packing algorithm
           packing_ranget(c,r)=ldas%udef
          enddo
         enddo

         do m=1,nvarst
          kmvart(m)=ldas%maxt
         enddo

         do m=1,ldas%maxt
          levst(m)=m
         enddo

!=== Create HDF file for 3D tile space output========================================

!=== Create HDF file for tile space output===========================================

        endif
       endif       !End HDF/Grib daily output file setup

       if(ldas%whdf.eq.1)then

        if(ldas%wtil.eq.1)then
!=== Write 3D tile output in HDF format with the subroutine gfio_putvar
         kbeg=1                 !First level to write; if 2D grid kbeg=0
         kount=nlevsoi       !number of levels to write

!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

!=== Write to HDF
         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &     tile(t)%pveg.eq.1)then
           do m=1,nlevsoi
            tempvarts3(tile(t)%col,tile(t)%row,m)=clm1(t)%t_soisno(m)
           enddo
          endif
         enddo

!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &    tile(t)%pveg.eq.2)then
           do m=1,nlevsoi
            tempvarts3(tile(t)%col,tile(t)%row,m)=clm1(t)%t_soisno(m)
           enddo
          endif
         enddo

!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &     tile(t)%pveg.eq.3)then
           do m=1,nlevsoi
            tempvarts3(tile(t)%col,tile(t)%row,m)=clm1(t)%t_soisno(m)
           enddo
          endif
         enddo

!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

!=== Write to HDF
         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &     tile(t)%pveg.eq.1)then
           do m=1,nlevsoi
            tempvarts3(tile(t)%col,tile(t)%row,m)=soilm(t,m)
           enddo
          endif
         enddo
!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

!=== Write to HDF
         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &     tile(t)%pveg.eq.2)then
           do m=1,nlevsoi
            tempvarts3(tile(t)%col,tile(t)%row,m)=soilm(t,m)
           enddo
          endif
         enddo

!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

!=== Write to HDF
         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &     tile(t)%pveg.eq.3)then
           do m=1,nlevsoi
            tempvarts3(tile(t)%col,tile(t)%row,m)=soilm(t,m)
           enddo
          endif
         enddo

!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

!=== Write to HDF
         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &     tile(t)%pveg.eq.1)then
           do m=1,nlevsoi 
            tempvarts3(tile(t)%col,tile(t)%row,m)=clm1(t)%h2osoi_liq(m)
           enddo
          endif
         enddo

!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

!=== Write to HDF      
         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &     tile(t)%pveg.eq.2)then
           do m=1,nlevsoi
            tempvarts3(tile(t)%col,tile(t)%row,m)=clm1(t)%h2osoi_liq(m)

           enddo
          endif
         enddo

!=== Initialize temporary tile array
         do m=1,nlevsoi
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts3(c,r,m)=ldas%udef
           enddo
          enddo 
         enddo

!=== Write to HDF 
         do t=1,ldas%nch
          if(grid(tile(t)%col,tile(t)%row)%imask.ge.1.and.
     &     tile(t)%pveg.eq.3)then
           do m=1,nlevsoi
            tempvarts3(tile(t)%col,tile(t)%row,m)=clm1(t)%h2osoi_liq(m)
           enddo
          endif
         enddo

!=== End 3D tile space output
         if((24-ldas%gmt).le.ldas%writeintc1.or.
     &     ldas%endtime.eq.1)then
         endif

!=== Write tile output in HDF format with the subroutine gfio_putvar
         kbeg=1                 !First level to write; if 2D grid kbeg=0
         kount=ldas%maxt        !Number of levels to write

!=== Initialize temporary tile array
         do m=1,ldas%maxt
          do r=1,ldas%nr
           do c=1,ldas%nc
            tempvarts(c,r,m)=ldas%udef
           enddo
          enddo
         enddo

         do t=1,ldas%nch
          tempvarts(tile(t)%col,tile(t)%row,tile(t)%pveg)=
     &     float(tile(t)%vegt)
         enddo
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(2),kbeg,kount,tile%fgrd,                   !Fraction of grid covered by tile
     &    ldas,grid,tile)
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(3),kbeg,kount,
     &  clm1%totfsa,              !Net Surface Shortwave Radiation (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(4),kbeg,kount,
     &  clm1%toteflx_lwrad_net,            !Net Surface Longwave Radiation (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(6),kbeg,kount,
     &  clm1%toteflx_lh_tot,             !Latent Heat Flux (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(6),kbeg,kount,
     &  clm1%toteflx_sh_tot,              !Sensible Heat Flux (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(7),kbeg,kount,
     &  clm1%toteflx_soil_grnd,              !Ground Heat Flux (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(8),kbeg,kount,
     &  snowflux,           !Snow Phase Change Heat Flux (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(9),kbeg,kount,
     &  clm1%totsolisbd,            !Downward Surface Shortwave Radiation (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(10),kbeg,kount,
     &  clm1%totforc_lwrad,               !Downward Surface Longwave Radiation (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(11),kbeg,kount,
     &  clm1%totsnow,              !Snowfall (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(12),kbeg,kount,
     &  clm1%totrain,              !Rainfall (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(13),kbeg,kount,
     &  clm1%totqflx_evap,            !Total Evaporation (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(14),kbeg,kount,
     &  clm1%totqflx_surf,              !Surface Runoff (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(15),kbeg,kount,
     &  clm1%totqflx_drain,            !Subsurface Runoff (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(16),kbeg,kount,
     &  snowmelt,          !Snowmelt (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(17),kbeg,kount,
     &  snowtemp,               !Snow Temperature (K)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(18),kbeg,kount,
     &  clm1%t_veg,                   !Canopy Temperature (K)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(19),kbeg,kount,
     &  clm1%t_grnd,                  !Bare Soil Temperature (K)
     &    ldas,grid,tile)
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(20),kbeg,kount,
     &  asurft,                     !Average Surface Temperature (K)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(21),kbeg,kount,
     &  clm1%t_rad,                   !Effective Radiative Surface Temperature (K)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(22),kbeg,kount,
     &  clm1%surfalb,                 !Surface Albedo, All Wavelengths (%)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(23),kbeg,kount,
     &  clm1%h2osno,                    !Snowpack Water Equivalent (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(24),kbeg,kount,
     &  clm1%h2ocan,                   !Plant Canopy Surface Water Storage (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(25),kbeg,kount,
     &  soilmtc,                !Total Column Soil Moisture (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(26),kbeg,kount,
     &  soilmr,                 !Root Zone Soil Moisture (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(27),kbeg,kount,
     &  soilm1m,             !Top 1-meter Soil Moisture (kg/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(28),kbeg,kount,
     &  soilwtc,                !Total Soil Column Wetness (%)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(29),kbeg,kount,
     &  soilwr,              !Root Zone Wetness (%)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(30),kbeg,kount,
     &  clm1%totqflx_ecanop,            !Canopy Surface Water Evaporation (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(31),kbeg,kount,
     &  cantrn,            !Canopy Transpiration (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(32),kbeg,kount,
     &  bare,              !Bare Soil Evaporation (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(33),kbeg,kount,
     &  snowevp,           !Snow Evaporation (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(34),kbeg,kount,
     &  potevp,            !Potential Evaporation (W/m2)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(35),kbeg,kount,
     &  clm1%acond,                  !Aerodynamic Conductance (m/s)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(36),kbeg,kount,
     &  ccond,                  !Canopy Conductance (m/s)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(37),kbeg,kount,
     &  clm1%tlai,                    !Leaf Area Index
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(38),kbeg,kount,
     &  clm1%snowdp,                 !Snow Depth (m)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(39),kbeg,kount,
     &  clm1%frac_sno,                   !Snow Cover (%)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(40),kbeg,kount,
     &  clm1%snoalb,                  !Snow Albedo (%)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(41),kbeg,kount,
     &  clm1%forc_t,                     !Two Meter Temperature (K)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(42),kbeg,kount,
     &  clm1%forc_q,                     !Two Meter Humidity (kg/kg)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(43),kbeg,kount,
     &  clm1%forc_u,                     !Ten Meter U Wind (m/s)
     &    ldas,grid,tile) 
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(44),kbeg,kount,
     &  clm1%forc_v,                     !Ten Meter V Wind (m/s)
     &    ldas,grid,tile)
         call tshdf(ldas%fidtc,ldas%rctc,vnamet(45),kbeg,kount,
     &  clm1%forc_pbot/100.,              !Surface Pressure (mb)
     &    ldas,grid,tile)


         if((24-ldas%gmt).le.ldas%writeintc1.or.
     &     ldas%endtime.eq.1)then
         endif
        endif   !End tile space output

!=== Write grid output in HDF format with the subroutine gfio_putvar
        kbeg=0     ! For 2D fields
        kount=1

        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(1),kbeg,kount,
     &   g_fsa,      !Net Surface Shortwave  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(2),kbeg,kount,
     &   g_eflx_lwrad_net,    !Net Longwave  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(3),kbeg,kount,
     &   g_eflx_lh_tot,      !Latent Heat Flux  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(4),kbeg,kount,
     &   g_eflx_sh_tot,      !Sensible Heat Flux  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(5),kbeg,kount,
     &   g_eflx_soil_grnd,      !Ground Heat Flux  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(6),kbeg,kount,g_snowflux,   !Snow Phase Change Heat Flux [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(7),kbeg,kount,g_solisbd,   !Downward Surface Shortwave  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(8),kbeg,kount,
     &   g_forc_lwrad,        !Downward Surface Longwave
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(9),kbeg,kount,g_snow,       !Snowfall  [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(10),kbeg,kount,g_rain,      !Rainfall  [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(11),kbeg,kount,
     &   g_qflx_evap,    !Evaporation [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(12),kbeg,kount,
     &   g_qflx_surf,      !Surface Runoff [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(13),kbeg,kount,
     &  g_qflx_drain,      !Subsurface Runoff [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(14),kbeg,kount,
     &  g_snowmelt,  !Snowmelt [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(15),kbeg,kount,
     &  g_snowtemp,  !Snow Temperature [K]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(16),kbeg,kount,g_t_veg,      !Canopy Temperature [K]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(17),kbeg,kount,g_t_grnd,        !Bare Soil Temperature [K]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(18),kbeg,kount,g_asurft,        !Average Surface Temperature [K]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(19),kbeg,kount,g_t_rad,      !Effective Radiative Temperature [K]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(20),kbeg,kount,g_surfalb,    !Surface Albedo [%]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(21),kbeg,kount,g_h2osno,       !Snowpack Water Equivalent [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(22),kbeg,kount,g_h2ocan,      !Plant Canopy Sfc Water Storage [kg/m2]
     &   ldas,grid)

        kbeg=1                 !First level to write; if 2D grid kbeg=0
        kount=nlevsoi       !number of levels to write

        do m=1,nlevsoi
         do c=1,ldas%nc
          do r=1,ldas%nr
           if(grid(c,r)%imask.eq.0) g_t_soisno(c,r,m)=ldas%udef
          enddo
         enddo
        enddo

        kbeg=0     ! For 2D fields
        kount=1

        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(24),kbeg,kount,g_soilmtc,      !Total Column Soil Moisture  [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(25),kbeg,kount,g_soilmr,       !Root Zone Soil Moisture  [kg/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(26),kbeg,kount,g_soilm1m,      !Top 1-meter Soil Moisture  [kg/m2]
     &   ldas,grid)

        kbeg=1                 !First level to write; if 2D grid kbeg=0
        kount=nlevsoi       !number of levels to write
        do m=1,nlevsoi
         do c=1,ldas%nc
          do r=1,ldas%nr
           if(grid(c,r)%imask.eq.0) g_soilm(c,r,m)=ldas%udef
          enddo
         enddo
        enddo

        do m=1,nlevsoi
         do c=1,ldas%nc
          do r=1,ldas%nr
           if(grid(c,r)%imask.eq.0) g_h2osoi_liq(c,r,m)=ldas%udef
          enddo
         enddo
        enddo

        kbeg=0     ! For 2D fields
        kount=1

        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(29),kbeg,kount,g_soilwtc,      !Total Column Soil Wetness  [%]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(30),kbeg,kount,g_soilwr,       !Root Zone Soil Wetness  [%]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(31),kbeg,kount,
     &   g_qflx_ecanop,       !Canopy Sfc Water Evaporation  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(32),kbeg,kount,g_cantrn,       !Canopy Transpiration  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(33),kbeg,kount,g_bare,         !Bare Soil Evaporation  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(34),kbeg,kount,g_snowevp,      !Snow Evaporation  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(35),kbeg,kount,g_potevp,       !Potential Evaporation  [W/m2]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(36),kbeg,kount,g_acond,        !Aerodynamic Conductance  [m/s]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(37),kbeg,kount,g_ccond,        !Canopy Conductance  [m/s]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(38),kbeg,kount,g_tlai,          !Leaf Area Index  [Unitless]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(39),kbeg,kount,g_snowdp,         !Snow Depth  [M]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(40),kbeg,kount,
     &   g_frac_sno,        !Snow Cover [%]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(41),kbeg,kount,g_snowa,        !Snow Albedo [%]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(42),kbeg,kount,g_forc_t,        !Two Meter Temperature [K]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(43),kbeg,kount,g_forc_q,        !Two Meter Humidity [kg/kg]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(44),kbeg,kount,g_forc_u,        !Ten Meter U Wind [m/s]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(45),kbeg,kount,g_forc_v,        !Ten Meter V Wind [m/s]
     &   ldas,grid)
        call gshdf(ldas%fidgc,ldas%rcgc,vnameg(46),kbeg,kount,
     &   g_forc_pbot,      !Surface Pressure [mb]
     &   ldas,grid)


         if((24-ldas%gmt).le.ldas%writeintc1.or.
     &     ldas%endtime.eq.1)then
         ldas%numoutc1=0    !Reset output time counter
        endif

       endif

!=== Write grid output in GRIB

        IF(LDAS%WGRB.EQ.1)THEN
        CALL GRIBOUTclm1 (LDAS,GRID,TILE,clm1,FBASE,FYRMODIR)
        ENDIF


!=== Write grid output in binary

       if(ldas%wbin.eq.1)then
        open(57,file=filengb,form='unformatted')

        write(57) g_fsa
        write(57) g_eflx_lwrad_net
        write(57) g_eflx_lh_tot
        write(57) g_eflx_sh_tot
        write(57) g_eflx_soil_grnd
        write(57) g_snowflux
        write(57) g_solisbd
        write(57) g_forc_lwrad
        write(57) g_snow
        write(57) g_rain
        write(57) g_qflx_evap
        write(57) g_qflx_surf
        write(57) g_qflx_drain
        write(57) g_snowmelt
        write(57) g_snowtemp
        write(57) g_t_veg
        write(57) g_t_grnd
        write(57) g_asurft
        write(57) g_t_rad
        write(57) g_surfalb
        write(57) g_h2osno
        write(57) g_h2ocan

        do m=1,nlevsoi
         do c=1,ldas%nc
          do r=1,ldas%nr
           tempvarb(c,r)=g_t_soisno(c,r,m)
          enddo
         enddo
         write(57)tempvarb
        enddo

        write(57) g_soilmtc
        write(57) g_soilmr
        write(57) g_soilm1m

        do m=1,nlevsoi
         do c=1,ldas%nc
          do r=1,ldas%nr
           tempvarb(c,r)=g_soilm(c,r,m)
          enddo
         enddo
         write(57)tempvarb
        enddo

        do m=1,nlevsoi
         do c=1,ldas%nc
          do r=1,ldas%nr
           tempvarb(c,r)=g_h2osoi_liq(c,r,m)
          enddo
         enddo
         write(57)tempvarb
        enddo

        write(57) g_soilwtc
        write(57) g_soilwr
        write(57) g_qflx_ecanop
        write(57) g_cantrn
        write(57) g_bare
        write(57) g_snowevp
        write(57) g_potevp
        write(57) g_acond
        write(57) g_ccond
        write(57) g_tlai
        write(57) g_snowdp
        write(57) g_frac_sno
        write(57) g_snowa
        write(57) g_forc_t
        write(57) g_forc_q
        write(57) g_forc_u
        write(57) g_forc_v
        write(57) g_forc_pbot

        close(57)
       endif

!=== Write statistical output
      if(ldas%clm1open.eq.0)then
       file='clm1stats.dat'
       call openfile(name,ldas%odir,ldas%expcode,file)
       if(ldas%startcode.eq.1)then
        open(60,file=name,form='formatted',status='unknown',
     1   position='append')
       else
        open(60,file=name,form='formatted',status='replace')
       endif
       ldas%clm1open=1
      ENDIF

       write(60,996)'       Statistical Summary of clm1 Output for:  ',
     & ldas%mo,'/',ldas%da,'/',ldas%yr,ldas%hr,':',
     & ldas%mn,':',ldas%ss
996    format(a47,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
       write(60,*)
       write(60,997)
!997    format(t28,'Mean',t42,'StDev',t56,'Min',t70,'Max') 
997    format(t36,'Mean',t50,'StDev',t64,'Min',t78,'Max')
       do t=1,ldas%nch
        tempvar(t)=clm1(t)%t_soisno(5)
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'T_SOISNO05(K):        ',vmean,vstdev,vmin,vmax
       do t=1,ldas%nch
        tempvar(t)=clm1(t)%t_soisno(3)
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'T_SOISNO03 (K):        ',vmean,vstdev,vmin,vmax
       do t=1,ldas%nch
        tempvar(t)=clm1(t)%t_soisno(1)
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'T_SOISNO01 (K):        ',vmean,vstdev,vmin,vmax

       do t=1,ldas%nch
        tempvar(t)=(clm1(t)%h2osoi_liq(5))/
     &   (denh2o*clm1(t)%dz(5))
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'H2OSOI_LIQ05 (Vol):    ',
     &    vmean,vstdev,vmin,vmax
       do t=1,ldas%nch
        tempvar(t)=(clm1(t)%h2osoi_liq(3))/
     &   (denh2o*clm1(t)%dz(3))
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'H2OSOI_LIQ03 (Vol):    ',
     &    vmean,vstdev,vmin,vmax
       do t=1,ldas%nch
        tempvar(t)=(clm1(t)%h2osoi_liq(1))/
     &   (denh2o*clm1(t)%dz(1))
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'H2OSOI_LIQ01 (Vol):    ',
     &    vmean,vstdev,vmin,vmax

       do t=1,ldas%nch
        tempvar(t)=(clm1(t)%h2osoi_ice(5))/
     &   (denice*clm1(t)%dz(5))
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'H2OSOI_ICE05 (Vol):    ',
     &    vmean,vstdev,vmin,vmax
       do t=1,ldas%nch
        tempvar(t)=(clm1(t)%h2osoi_ice(3))/
     &   (denice*clm1(t)%dz(3))
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'H2OSOI_ICE03 (Vol):    ',
     &    vmean,vstdev,vmin,vmax
       do t=1,ldas%nch
        tempvar(t)=(clm1(t)%h2osoi_ice(1))/
     &   (denice*clm1(t)%dz(1))
       enddo
       call stats(tempvar,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'H2OSOI_ICE01 (Vol):    ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%totfsa,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'FSA (W/m2):            ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%toteflx_lwrad_net,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'EFLX_LWRAD_NET (W/m2): ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%toteflx_lh_tot,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'EFLX_LH_TOT (W/m2):    ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%toteflx_sh_tot,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'EFLX_SH_TOT (W/m2):    ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%toteflx_soil_grnd,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'EFLX_SOIL_GRND (W/m2): ',
     &    vmean,vstdev,vmin,vmax
       call stats(snowflux,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SNOWFLUX (W/m2):       ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%totsolisbd,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'SOLISBD (W/m2):        ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%totforc_lwrad,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'FORC_LWRAD (W/m2):     ',
     &    vmean,vstdev,vmin,vmax
      call stats(clm1%totsnow,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,998)'SNOW (kg/m2):          ',
     &    vmean,vstdev,vmin,vmax
      call stats(clm1%totrain,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,998)'RAIN (kg/m2):          ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%totqflx_evap,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,998)'QFLX_EVAP (kg/m2):     ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%totqflx_surf,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,998)'QFLX_SURF (kg/m2):     ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%totqflx_drain,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,998)'QFLX_DRAIN (kg/m2):    ',
     &    vmean,vstdev,vmin,vmax
       call stats(snowmelt,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,998)'QFLX_SNOMELT (kg/m2):  ',
     &    vmean,vstdev,vmin,vmax
       call stats(snowtemp,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SNOWTEMP (K):          ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%t_veg,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'T_VEG (K):             ',vmean,vstdev,vmin,vmax
       call stats(clm1%t_grnd,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'T_GRND (K):            ',vmean,vstdev,vmin,vmax
       call stats(asurft,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'ASURFT (K):            ',vmean,vstdev,vmin,vmax
       call stats(clm1%t_rad,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'T_RAD (K):             ',vmean,vstdev,vmin,vmax
      call stats(clm1%surfalb,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SURFALB (-):           ',vmean,vstdev,vmin,vmax
       call stats(clm1%h2ocan,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,998)'H2OCAN (kg/m2):        ',vmean,vstdev,vmin,vmax
       call stats(clm1%h2osno,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,998)'H2OSNO (kg/m2):        ',vmean,vstdev,vmin,vmax

       call stats(soilmtc,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SOILMTC (kg/m2):       ',vmean,vstdev,vmin,vmax
       call stats(soilmr,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SOILMR (kg/m2):        ',vmean,vstdev,vmin,vmax
       call stats(soilm1m,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SOILM1M (kg/m2):       ',vmean,vstdev,vmin,vmax
       call stats(soilwtc,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SOILWTC (%):           ',vmean,vstdev,vmin,vmax
       call stats(soilwr,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SOILWR (%):            ',vmean,vstdev,vmin,vmax
       call stats(clm1%totqflx_ecanop,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'QFLX_ECANOP (W/m2):    ',
     &    vmean,vstdev,vmin,vmax
       call stats(cantrn,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'QFLX_TRAN_VEG (W/m2):  ',
     &    vmean,vstdev,vmin,vmax
       call stats(bare,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'QFLX_EVAP_GRND (W/m2): ',
     &    vmean,vstdev,vmin,vmax
       call stats(snowevp,ldas%udef,ldas%nch,
     &    vmean,vstdev,vmin,vmax)
       write(60,999)'QFLX_SUB_SNO (W/m2):   ',
     &    vmean,vstdev,vmin,vmax
       call stats(potevp,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'POTEVP (W/m2):         ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%acond,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,998)'ACOND (m/s):           ',
     &    vmean,vstdev,vmin,vmax
       call stats(ccond,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'CCOND (m/s):           ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%tlai,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'TLAI (-):              ',
     &    vmean,vstdev,vmin,vmax
       call stats(clm1%snowdp,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SNOWDP (m):            ',vmean,vstdev,vmin,vmax
       call stats(clm1%frac_sno,ldas%udef,
     &    ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'FRAC_SNO (%):          ',vmean,vstdev,vmin,vmax
       call stats(clm1%snoalb,ldas%udef,ldas%nch,vmean,vstdev,vmin,vmax)
       write(60,999)'SNOALB (%):            ',vmean,vstdev,vmin,vmax

999    format (1x,a25,4f14.3)
998    format (1x,a25,4e14.3)
       write(60,*)
       write(60,*)
      
!=== Reinitialize total arrays to zero
      clm1%totfsa=0.                   
      clm1%toteflx_lwrad_net=0.
      clm1%toteflx_lh_tot=0.
      clm1%toteflx_sh_tot=0.
      clm1%toteflx_soil_grnd=0.
      clm1%totqflx_snomelt=0.
      clm1%totsolisbd=0.
      clm1%totforc_lwrad=0.
      clm1%totsnow=0.
      clm1%totrain=0.
      clm1%totqflx_evap=0.
      clm1%totqflx_surf=0.
      clm1%totqflx_drain=0.
      clm1%totqflx_ecanop=0.
      clm1%totqflx_tran_veg=0.
      clm1%totqflx_evap_grnd=0.
      clm1%totqflx_sub_snow=0.
 
      clm1%count=0 

!=== Reinitialize grid space arrays to zero
       do m=1,nlevsoi
        do c=1,ldas%nc
         do r=1,ldas%nr
          g_t_soisno(c,r,m)=0.
          g_h2osoi_liq(c,r,m)=0.
          g_h2osoi_ice(c,r,m)=0.
          g_soilm(c,r,m)=0.
         enddo
        enddo
       enddo


        g_fsa=0.
        g_eflx_lwrad_net=0.
        g_eflx_lh_tot=0.
        g_eflx_sh_tot=0.
        g_eflx_soil_grnd=0.
        g_snowflux=0.
        g_solisbd=0.
        g_forc_lwrad=0.
        g_snow=0.
        g_rain=0.
        g_qflx_evap=0.
        g_qflx_surf=0.
        g_qflx_drain=0.
        g_snowmelt=0.
        g_snowtemp=0.
        g_t_veg=0.
        g_t_grnd=0.
        g_asurft=0.
        g_t_rad=0.
        g_surfalb=0.
        g_h2osno=0.
        g_h2ocan=0.
        g_soilmtc=0.
        g_soilmr=0.
        g_soilm1m=0.
        g_soilwtc=0.
        g_soilwr=0.
        g_qflx_ecanop=0.
        g_cantrn=0.
        g_bare=0.
        g_snowevp=0.
        g_potevp=0.
        g_acond=0.
        g_ccond=0.
        g_tlai=0.
        g_snowdp=0.
        g_frac_sno=0.
        g_snowa=0.
        g_forc_t=0.
        g_forc_q=0.
        g_forc_u=0.
        g_forc_v=0.
        g_forc_pbot=0.


      endif
      end subroutine clm1_out

