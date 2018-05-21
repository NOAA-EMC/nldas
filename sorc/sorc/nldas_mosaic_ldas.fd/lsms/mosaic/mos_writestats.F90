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
! !ROUTINE: mos_writevars.F90
!
! !DESCRIPTION:  
!  LIS MOS data writer: Writes mos output
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine mos_writestats(ftn_stats)
! !USES:
  use lisdrv_module, only : lis
  use mos_varder
 
  implicit none
  
  integer :: ftn_stats,t
!EOP
  real :: vmean,vstdev,vmin,vmax
  real :: rainf(lis%d%glbnch)
  real :: snowf(lis%d%glbnch)
!BOC
  do t=1,lis%d%glbnch
     if(mos(t)%forcing(1) < 273.15) then
        rainf(t) = 0.0
        snowf(t) = mos(t)%forcing(8)
     else
        rainf(t) = mos(t)%forcing(8)
        snowf(t) = 0.0
     endif
  enddo
!---------------------------------------------------------------------------
! General Energy Balance Components
!---------------------------------------------------------------------------
   call stats(mos%swnet,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'SWnet(W/m2)', &
        vmean,vstdev,vmin,vmax
   call stats(mos%lwnet,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'LWnet(W/m2)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%qle,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'Qle(W/m2)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%qh,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'Qh(W/m2)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%qg,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'Qg(W/m2)',&
        vmean,vstdev,vmin,vmax
!---------------------------------------------------------------------------
! General Water Balance Components
!---------------------------------------------------------------------------
   call stats(mos%snowf,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'Snowf(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%rainf,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'Rainf(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%evap,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'Evap(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%qs,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'Qs(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%qsb,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'Qsb(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats((mos%water1+ &
        mos%water2+ & 
        mos%water3  & 
        -mos%soilm_prev)/float(mos%count), & 
        lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, vmax)
   write(ftn_stats,999) 'DelSoilMoist(kg/m2s)', & 
        vmean,vstdev,vmin,vmax
   call stats((mos%snow-mos%swe_prev)/float(mos%count), & 
        lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, vmax)
   write(ftn_stats,999) 'DelSWE(kg/m2s)', & 
        vmean,vstdev,vmin,vmax
!---------------------------------------------------------------------------
! Surface State Variables
!---------------------------------------------------------------------------
   call stats(mos%avgsurft,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'AvgSurfT(K)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%albedo,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'Albedo(-)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%swe,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'SWE(kg/m2)',&
        vmean,vstdev,vmin,vmax
!---------------------------------------------------------------------------
! Subsurface State Variables
!---------------------------------------------------------------------------
   call stats(mos%soilmoist1,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'SoilMoist1(kg/m2)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%soilmoist2,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'SoilMoist2(kg/m2)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%soilmoist3,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,999) 'SoilMoist3(kg/m2)',&
        vmean,vstdev,vmin,vmax
   write(ftn_stats,998) 'SoilWet(-)',&
        vmean,vstdev,vmin,vmax
!---------------------------------------------------------------------------
! Evaporation Components
!---------------------------------------------------------------------------
   call stats(mos%ecanop,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'ECanop(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%tveg,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'TVeg(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%esoil,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'ESoil(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%rootmoist,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'RootMoist(kg/m2)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%canopint,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'Canopint(kg/m2s)',&
        vmean,vstdev,vmin,vmax
   call stats(mos%acond,lis%d%udef,lis%d%glbnch,vmean, & 
        vstdev,vmin, vmax)
   write(ftn_stats,998) 'Acond(m/s)',&
        vmean,vstdev,vmin,vmax
   if(lis%o%wfor.eq.1) then
      call stats(sqrt(mos%forcing(5)*mos%forcing(5)+ & 
           mos%forcing(6)*mos%forcing(6)), & 
           lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, vmax)
      write(ftn_stats,999) 'Wind(m/s)', & 
           vmean,vstdev,vmin,vmax
      call stats(rainf, & 
           lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, vmax)
      write(ftn_stats,998) 'Rainf(kg/m2s)', & 
           vmean,vstdev,vmin,vmax
      call stats(snowf, & 
           lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, vmax)
      write(ftn_stats,998) 'Snowf(kg/m2s)', & 
           vmean,vstdev,vmin,vmax
      call stats(mos%forcing(1),lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'Tair(K)', & 
           vmean,vstdev,vmin,vmax
      call stats(mos%forcing(2),lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'Qair(kg/kg)', & 
           vmean,vstdev,vmin,vmax
      call stats(mos%forcing(7),lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'PSurf(Pa)', & 
           vmean,vstdev,vmin,vmax
      call stats(mos%forcing(3),lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'SWdown (W/m2)', & 
           vmean,vstdev,vmin,vmax
      call stats(mos%forcing(4),lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'LWdown(W/m2)', & 
           vmean,vstdev,vmin,vmax
   endif
998    FORMAT(1X,A18,4E14.3)
999    FORMAT(1X,A18,4F14.3)
!EOC
 end subroutine mos_writestats
 
