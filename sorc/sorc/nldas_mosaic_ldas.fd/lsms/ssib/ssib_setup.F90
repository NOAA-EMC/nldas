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
! !ROUTINE: ssib_setup.F90
!
! !DESCRIPTION:
!  Complete the setup routines for SSiB
!
! !REVISION HISTORY:
!  4 Nov 1999: Paul Houser, Initial Code
! 28 Apr 2002: Kristi Arsenault, Modified to SSiB LSM 2.5 code to LDAS
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
! 22 May 2004: David Mocko, made compatible with SiB-lings
! 
! !INTERFACE:
      subroutine ssib_setup
! !USES:
      use lisdrv_module, only : lis,tile
      use ssib_varder
      use spmdMod, only : masterproc, npes
!EOP
      implicit none

      integer :: t, n
!BOC
!=== End Variable List ===================================================
#if ( ! defined OPENDAP )
      if ( masterproc ) then
#endif

      call setssibp
      call ssib_gfrac
!      call ssib_alb(lis%d, lis%t, tile)  
      call ssib_coldstart

      do t=1,lis%d%nch
         ssib(t)%swnet = 0
         ssib(t)%lwnet = 0
         ssib(t)%qle = 0
         ssib(t)%qh = 0
         ssib(t)%qg = 0
         ssib(t)%qf = 0
         ssib(t)%qtau = 0
         ssib(t)%delsurfheat = 0
         ssib(t)%snowf = 0
         ssib(t)%rainf = 0
         ssib(t)%evap = 0
         ssib(t)%qs = 0
         ssib(t)%qsb = 0
         ssib(t)%qsm = 0
         ssib(t)%delsoilmoist = 0
         ssib(t)%delswe = 0
         ssib(t)%delintercept = 0
         ssib(t)%vegtc = 0
         ssib(t)%baresoilt = 0
         ssib(t)%avgsurft = 0
         ssib(t)%radteff = 0
         ssib(t)%albedo = 0
         ssib(t)%swe = 0
         ssib(t)%sweveg = 0
         ssib(t)%soilmoist1 = 0
         ssib(t)%soilmoist2 = 0
         ssib(t)%soilmoist3 = 0
         ssib(t)%soiltemp = 0
         ssib(t)%soilwet = 0
         ssib(t)%ecanop = 0
         ssib(t)%tveg = 0
         ssib(t)%esoil = 0
         ssib(t)%rootmoist = 0
         ssib(t)%canopint = 0
         ssib(t)%acond = 0
         ssib(t)%snowfrac = 0

         ssib(t)%count = 0
         ssib(t)%albedocount = 0
      enddo

#if ( ! defined OPENDAP )
      endif

      if ( npes > 1 ) then
         call ssib_scatter
      endif
#endif
      return
!EOC
      end subroutine ssib_setup

