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
! !ROUTINE: ssib_gridout.F90
!
! !DESCRIPTION:
!  LIS SSiB data writer: Writes SSiB output in grid space
!
! !REVISION HISTORY:
! 02 Dec 2003: Sujay Kumar, Initial Version
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
! 22 May 2004: David Mocko, made compatible with SiB-lings
!
! !INTERFACE:
      subroutine ssib_binout(ftn)
! !USES:
      use lisdrv_module, only : lis
      use drv_output_mod, only : drv_writevar_bin
      use ssib_varder

      implicit none
! !ARGUMENTS:
      integer :: ftn
!EOP
      real :: gtmp(lis%d%lnc,lis%d%lnr)
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
      ssib%swnet = ssib%swnet/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%swnet)

      ssib%lwnet = ssib%lwnet/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%lwnet)

      ssib%qle = ssib%qle/float(ssib%count)
      call drv_writevar_bin(ftn,ssib%qle)

      ssib%qh = ssib%qh/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%qh)

      ssib%qg = ssib%qg/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%qg)

      ssib%qf = ssib%qf/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%qf)

      ssib%qtau = ssib%qtau/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%qtau)

      call drv_writevar_bin(ftn, ssib%delsurfheat)
!-----------------------------------------------------------------------
! General Water Balance Components
!-----------------------------------------------------------------------
      ssib%snowf = ssib%snowf/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%snowf)

      ssib%rainf = ssib%rainf/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%rainf)

      ssib%evap = ssib%evap/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%evap)

      ssib%qs = ssib%qs/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%qs)

      ssib%qsb = ssib%qsb/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%qsb)

      ssib%qsm = ssib%qsm/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%qsm)

      call drv_writevar_bin(ftn, ssib%delsoilmoist)

      call drv_writevar_bin(ftn, ssib%delswe)

      call drv_writevar_bin(ftn, ssib%delintercept)

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

      call drv_writevar_bin(ftn, ssib%vegtc)

      call drv_writevar_bin(ftn, ssib%baresoilt)

      call drv_writevar_bin(ftn, ssib%avgsurft)

      call drv_writevar_bin(ftn, ssib%radteff)

      call drv_writevar_bin(ftn, ssib%albedo)

      call drv_writevar_bin(ftn, ssib%swe)

      call drv_writevar_bin(ftn, ssib%sweveg)

!-----------------------------------------------------------------------
! Subsurface State Variables
!-----------------------------------------------------------------------
      ssib%soilmoist1 = ssib%soilmoist1/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%soilmoist1)

      ssib%soilmoist2 = ssib%soilmoist2/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%soilmoist2)

      ssib%soilmoist3 = ssib%soilmoist3/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%soilmoist3)

      ssib%soiltemp = ssib%soiltemp/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%soiltemp)

      ssib%soilwet = ssib%soilwet/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%soilwet)

!-----------------------------------------------------------------------
! Evaporation Components
!-----------------------------------------------------------------------
      ssib%ecanop = ssib%ecanop/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%ecanop)

      ssib%tveg = ssib%tveg/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%tveg)

      ssib%esoil = ssib%esoil/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%esoil)

      ssib%rootmoist = ssib%rootmoist/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%rootmoist)

      ssib%canopint = ssib%canopint/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%canopint)

      ssib%acond = ssib%acond/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%acond)

      ssib%snowfrac = ssib%snowfrac/float(ssib%count)
      call drv_writevar_bin(ftn, ssib%snowfrac)

!-----------------------------------------------------------------------
! Cold Season Processes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Forcing Variables
!-----------------------------------------------------------------------
      if (lis%o%wfor.eq.1) then
         call drv_writevar_bin(ftn, sqrt(ssib%forcing(5)*ssib%forcing(5)+ & 
                        ssib%forcing(6)*ssib%forcing(6)))

         call drv_writevar_bin(ftn, rainf)

         call drv_writevar_bin(ftn, snowf)

         call drv_writevar_bin(ftn, ssib%forcing(1))

         call drv_writevar_bin(ftn, ssib%forcing(2))

         call drv_writevar_bin(ftn, ssib%forcing(7))

         call drv_writevar_bin(ftn, ssib%forcing(3))

         call drv_writevar_bin(ftn, ssib%forcing(4))

      endif

 998  FORMAT(1X,A18,4E14.3)
 999  FORMAT(1X,A18,4F14.3)
      return
!EOC
    end subroutine ssib_binout

