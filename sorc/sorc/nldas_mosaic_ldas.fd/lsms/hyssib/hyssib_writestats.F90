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
! !ROUTINE: hyssib_writestats
!
! !DESCRIPTION:
!  LIS HY-SSiB data writer: Writes HY-SSiB output in grid space
!
! !REVISION HISTORY:
! 02 Dec 2003: Sujay Kumar, Initial Version
! 21 Apr 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
subroutine hyssib_writestats(ftn_stats)
! !USES:
  use lisdrv_module, only : lis
  use hyssib_varder
  
  implicit none
! !ARGUMENTS:
  integer ::ftn_stats
!EOP
  real :: vmean,vstdev,vmin,vmax
  real :: rainf_in(lis%d%glbnch)
  real :: snowf_in(lis%d%glbnch)
  integer :: t
!BOC
  do t=1,lis%d%glbnch
     if (hyssib(t)%forcing(1).lt.273.15) then
        rainf_in(t) = 0.0
        snowf_in(t) = hyssib(t)%forcing(8)
     else
        rainf_in(t) = hyssib(t)%forcing(8)
        snowf_in(t) = 0.0
     endif
  enddo

!-----------------------------------------------------------------------
! General Energy Balance Components
!-----------------------------------------------------------------------

  call stats(hyssib%swnet,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SWnet(W/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%lwnet,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'LWnet(W/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qle,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'Qle(W/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qh,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'Qh(W/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qg,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'Qg(W/m2)',vmean,vstdev,vmin,vmax


  call stats(hyssib%qf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'Qf(W/m2)',vmean,vstdev,vmin,vmax


  call stats(hyssib%qv,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'Qv(W/m2)',vmean,vstdev,vmin,vmax


  call stats(hyssib%qtau,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'Qtau(N/m2)',vmean,vstdev,vmin,vmax


  call stats(hyssib%qa,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'Qa(W/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%delsurfheat,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'DelSurfHeat(W/m2)',vmean,vstdev,vmin,vmax


  call stats(hyssib%delcoldcont,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'DelColdCont(W/m2)',vmean,vstdev,vmin,vmax

!-----------------------------------------------------------------------
! General Water Balance Components
!-----------------------------------------------------------------------

  call stats(hyssib%snowf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Snowf(kg/m2s)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%rainf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Rainf(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%evap,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Evap(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qs,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Qs(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qrec,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Qrec(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qsb,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Qsb(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qsm,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Qsm(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qfz,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Qfz(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%qst,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'Qst(kg/m2s)',vmean,vstdev,vmin,vmax
      
  call stats(hyssib%delsoilmoist,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'DelSoilMoist(kg/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%delswe,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'DelSWE(kg/m2)',vmean,vstdev,vmin,vmax


  call stats(hyssib%delintercept,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'DelIntercept(kg/m2)',vmean,vstdev,vmin,vmax

!-----------------------------------------------------------------------
! Surface State Variables
!----------------------------------------------------------------------
  if (hyssibdrv%STATEVAR_AVG.eq.1) then
     do t = 1,lis%d%glbnch
        if (hyssib(t)%snowtcount.gt.0) then
           hyssib(t)%snowt = hyssib(t)%snowt/float(hyssib(t)%snowtcount)
        else
           hyssib(t)%snowt = lis%d%udef
        endif
        if (hyssib(t)%albedocount.gt.0) then
           hyssib(t)%albedo = hyssib(t)%albedo/float(hyssib(t)%albedocount)
        else
           hyssib(t)%albedo = lis%d%udef
        endif
        if (hyssib(t)%sliqfraccount.gt.0) then
           hyssib(t)%sliqfrac = hyssib(t)%sliqfrac/float(hyssib(t)%sliqfraccount)
        else
           hyssib(t)%sliqfrac = lis%d%udef
        endif
     enddo
     hyssib%vegtc = hyssib%vegtc/float(hyssib%count)
     hyssib%baresoilt = hyssib%baresoilt/float(hyssib%count)
     hyssib%avgsurft = hyssib%avgsurft/float(hyssib%count)
     hyssib%radteff = hyssib%radteff/float(hyssib%count)
     hyssib%swe = hyssib%swe/float(hyssib%count)
     hyssib%sweveg = hyssib%sweveg/float(hyssib%count)
  endif
  
  call stats(hyssib%snowt,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SnowT(K)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%vegtc,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'VegT(K)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%baresoilt,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'BaresoilT(K)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%avgsurft,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'AvgSurfT(K)',vmean,vstdev,vmin,vmax

  call stats(hyssib%radteff,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'RadT(K)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%albedo,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'Albedo(-)',vmean,vstdev,vmin,vmax

  call stats(hyssib%swe,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SWE(kg/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%sweveg,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SWEVeg(kg/m2)',vmean,vstdev,vmin,vmax

!-----------------------------------------------------------------------
! Subsurface State Variables
!-----------------------------------------------------------------------
  if (hyssibdrv%STATEVAR_AVG.eq.1) then
     hyssib%soilmoist1 = hyssib%soilmoist1/float(hyssib%count)
     hyssib%soilmoist2 = hyssib%soilmoist2/float(hyssib%count)
     hyssib%soilmoist3 = hyssib%soilmoist3/float(hyssib%count)
     hyssib%soiltemp = hyssib%soiltemp/float(hyssib%count)
     hyssib%soilwet = hyssib%soilwet/float(hyssib%count)
  endif
  
  call stats(hyssib%soilmoist1,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SoilMoist1(kg/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%soilmoist2,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SoilMoist2(kg/m2)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%soilmoist3,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SoilMoist3(kg/m2)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%soiltemp,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SoilTemp(K)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%soilwet,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SoilWet(-)',vmean,vstdev,vmin,vmax
  
!-----------------------------------------------------------------------
! Evaporation Components
!-----------------------------------------------------------------------
  call stats(hyssib%potevap,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'PotEvap(kg/m2s)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%ecanop,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'ECanop(kg/m2s)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%tveg,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'TVeg(kg/m2s)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%esoil,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'ESoil(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%rootmoist,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'RootMoist(kg/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%canopint,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'CanopInt(kg/m2)',vmean,vstdev,vmin,vmax

  call stats(hyssib%subsnow,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'SubSnow(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%subsurf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,998) 'SubSurf(kg/m2s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%acond,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'ACond(m/s)',vmean,vstdev,vmin,vmax

  call stats(hyssib%snowfrac,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SnowFrac(-)',vmean,vstdev,vmin,vmax
  
!-----------------------------------------------------------------------
! Cold Season Processes
!-----------------------------------------------------------------------
  call stats(hyssib%snowdepth,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SnowDepth(m)',vmean,vstdev,vmin,vmax
  
  call stats(hyssib%sliqfrac,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
  write(ftn_stats,999) 'SliqFrac(-)',vmean,vstdev,vmin,vmax
  
!-----------------------------------------------------------------------
! Forcing Variables
!-----------------------------------------------------------------------
  if (lis%o%wfor.eq.1) then
     call stats(sqrt(hyssib%forcing(5)*hyssib%forcing(5)+ & 
          hyssib%forcing(6)*hyssib%forcing(6)), & 
          lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,999) 'Wind(m/s)',vmean,vstdev,vmin,vmax
     
     call stats(rainf_in,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,998) 'Rainf(kg/m2s)',vmean,vstdev,vmin,vmax
     
     call stats(snowf_in,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,998) 'Snowf(kg/m2s)',vmean,vstdev,vmin,vmax
     
     call stats(hyssib%forcing(1),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,999) 'Tair(K)',vmean,vstdev,vmin,vmax
     
     call stats(hyssib%forcing(2),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,999) 'Qair(kg/kg)',vmean,vstdev,vmin,vmax
     
     call stats(hyssib%forcing(7),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,999) 'PSurf(Pa)',vmean,vstdev,vmin,vmax
     
     call stats(hyssib%forcing(3),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,999) 'SWdown (W/m2)',vmean,vstdev,vmin,vmax
     
     call stats(hyssib%forcing(4),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
     write(ftn_stats,999) 'LWdown(W/m2)',vmean,vstdev,vmin,vmax
  endif
  
998 FORMAT(1X,A18,4E14.3)
999 FORMAT(1X,A18,4F14.3)
  return
!EOC
end subroutine hyssib_writestats

