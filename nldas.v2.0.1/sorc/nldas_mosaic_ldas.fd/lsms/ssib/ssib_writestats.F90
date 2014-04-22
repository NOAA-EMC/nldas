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
! !ROUTINE: ssib_tileout.F90
!
! !DESCRIPTION:
!  LIS HY-SSiB data writer: Writes HY-SSiB output in tile space
!
! !REVISION HISTORY:
! 02 Dec 2003: Sujay Kumar, Initial Version
! 22 May 2004: David Mocko, Conversion from HY-SSiB to SSiB
!
! !INTERFACE:
subroutine ssib_writestats(ftn_stats)
! !USES:
  use lisdrv_module, only : lis, tile
  use ssib_varder

      implicit none
! !ARGUMENTS:
      integer :: ftn_stats
!EOP
      real :: gtmp(lis%d%glbngrid)
      real :: vmean,vstdev,vmin,vmax
      real :: rainf(lis%d%glbnch)
      real :: snowf(lis%d%glbnch)
      integer :: t
!BOC
      do t=1,lis%d%glbnch
         if (ssib(t)%forcing(1).lt.273.15) then
            rainf(t) = 0.0
            snowf(t) = ssib(t)%forcing(8)
         else
            rainf(t) = ssib(t)%forcing(8)
            snowf(t) = 0.0
         endif
      enddo

!-----------------------------------------------------------------------
! General Energy Balance Components
!-----------------------------------------------------------------------
      call stats(ssib%swnet,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SWnet(W/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%lwnet,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'LWnet(W/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%qle,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'Qle(W/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%qh,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'Qh(W/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%qg,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'Qg(W/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%qf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'Qf(W/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%qtau,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'Qtau(W/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%delsurfheat,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'DelSurfHeat(J/m2)',vmean,vstdev,vmin,vmax

!-----------------------------------------------------------------------
! General Water Balance Components
!-----------------------------------------------------------------------
      call stats(ssib%snowf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'Snowf(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%rainf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'Rainf(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%evap,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'Evap(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%qs,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'Qs(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%qsb,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'Qsb(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%qsm,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'Qsm(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%delsoilmoist,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'DelSoilMoist(kg/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%delswe,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'DelSWE(kg/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%delintercept,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'DelIntercept(kg/m2)',vmean,vstdev,vmin,vmax

!-----------------------------------------------------------------------
! Surface State Variables
!-----------------------------------------------------------------------
      if (ssibdrv%STATEVAR_AVG.eq.1) then
         do t = 1,lis%d%glbnch
            if (ssib(t)%albedocount.gt.0) then
               ssib(t)%albedo = ssib(t)%albedo/float(ssib(t)%albedocount)
            else
               ssib(t)%albedo = lis%d%udef
            endif
         enddo
         ssib%vegtc = ssib%vegtc/float(ssib%count)
         ssib%baresoilt = ssib%baresoilt/float(ssib%count)
         ssib%avgsurft = ssib%avgsurft/float(ssib%count)
         ssib%radteff = ssib%radteff/float(ssib%count)
         ssib%swe = ssib%swe/float(ssib%count)
         ssib%sweveg = ssib%sweveg/float(ssib%count)
      endif

      call stats(ssib%vegtc,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'VegT(K)',vmean,vstdev,vmin,vmax

      call stats(ssib%baresoilt,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'BaresoilT(K)',vmean,vstdev,vmin,vmax

      call stats(ssib%avgsurft,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'AvgSurfT(K)',vmean,vstdev,vmin,vmax

      call stats(ssib%radteff,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'RadT(K)',vmean,vstdev,vmin,vmax

      call stats(ssib%albedo,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'Albedo(-)',vmean,vstdev,vmin,vmax

      call stats(ssib%swe,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SWE(kg/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%sweveg,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SWEVeg(kg/m2)',vmean,vstdev,vmin,vmax

!-----------------------------------------------------------------------
! Subsurface State Variables
!-----------------------------------------------------------------------
      if (ssibdrv%STATEVAR_AVG.eq.1) then
         ssib%soilmoist1 = ssib%soilmoist1/float(ssib%count)
         ssib%soilmoist2 = ssib%soilmoist2/float(ssib%count)
         ssib%soilmoist3 = ssib%soilmoist3/float(ssib%count)
         ssib%soiltemp = ssib%soiltemp/float(ssib%count)
         ssib%soilwet = ssib%soilwet/float(ssib%count)
      endif

      call stats(ssib%soilmoist1,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SoilMoist1(kg/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%soilmoist2,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SoilMoist2(kg/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%soilmoist3,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SoilMoist3(kg/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%soiltemp,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SoilTemp(K)',vmean,vstdev,vmin,vmax

      call stats(ssib%soilwet,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SoilWet(-)',vmean,vstdev,vmin,vmax

!-----------------------------------------------------------------------
! Evaporation Components
!-----------------------------------------------------------------------
      call stats(ssib%ecanop,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'ECanop(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%tveg,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'TVeg(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%esoil,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'ESoil(kg/m2s)',vmean,vstdev,vmin,vmax

      call stats(ssib%rootmoist,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'RootMoist(kg/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%canopint,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,998) 'CanopInt(kg/m2)',vmean,vstdev,vmin,vmax

      call stats(ssib%acond,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'ACond(m/s)',vmean,vstdev,vmin,vmax

      call stats(ssib%snowfrac,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
      write(ftn_stats,999) 'SnowFrac(-)',vmean,vstdev,vmin,vmax

!-----------------------------------------------------------------------
! Cold Season Processes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Forcing Variables
!-----------------------------------------------------------------------
      if (lis%o%wfor.eq.1) then
         call stats(sqrt(ssib%forcing(5)*ssib%forcing(5)+ & 
                    ssib%forcing(6)*ssib%forcing(6)), & 
                    lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
         write(ftn_stats,999) 'Wind(m/s)',vmean,vstdev,vmin,vmax

         call stats(rainf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
         write(ftn_stats,998) 'Rainf(kg/m2s)',vmean,vstdev,vmin,vmax

         call stats(snowf,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
         write(ftn_stats,998) 'Snowf(kg/m2s)',vmean,vstdev,vmin,vmax

         call stats(ssib%forcing(1),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
         write(ftn_stats,999) 'Tair(K)',vmean,vstdev,vmin,vmax

         call stats(ssib%forcing(2),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
         write(ftn_stats,999) 'Qair(kg/kg)',vmean,vstdev,vmin,vmax

         call stats(ssib%forcing(7),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
         write(ftn_stats,999) 'PSurf(Pa)',vmean,vstdev,vmin,vmax

         call stats(ssib%forcing(3),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
         write(ftn_stats,999) 'SWdown (W/m2)',vmean,vstdev,vmin,vmax

         call stats(ssib%forcing(4),lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin,vmax)
         write(ftn_stats,999) 'LWdown(W/m2)',vmean,vstdev,vmin,vmax

      endif

 998  FORMAT(1X,A18,4E14.3)
 999  FORMAT(1X,A18,4F14.3)
      return
!EOC
    end subroutine ssib_writestats

