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
! !ROUTINE : hyssib$_-$setup.F90
!
! !DESCRIPTION:
!  Complete the setup routines for HY-SSiB
!
! !REVISION HISTORY:
!  4 Nov 1999: Paul Houser, Initial Code
! 28 Apr 2002: Kristi Arsenault, Modified to SSiB LSM 2.5 code to LDAS
! 21 Apr 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      subroutine hyssib_setup()
! !USES:
      use lisdrv_module, only : lis,tile
      use hyssib_varder
      use spmdMod, only : masterproc, npes
!EOP
      implicit none

      integer :: t, n
!BOC
!=== End Variable List ===================================================
#if ( ! defined OPENDAP )
      if ( masterproc ) then
#endif

      call sethyssibp()
      call hyssib_gfrac()
      call hyssib_albedo()
      call hyssib_coldstart()

      do t=1,lis%d%nch
         hyssib(t)%swnet = 0
         hyssib(t)%lwnet = 0
         hyssib(t)%qle = 0
         hyssib(t)%qh = 0
         hyssib(t)%qg = 0
         hyssib(t)%qf = 0
         hyssib(t)%qv = 0
         hyssib(t)%qtau = 0
         hyssib(t)%qa = 0
         hyssib(t)%delsurfheat = 0
         hyssib(t)%delcoldcont = 0
         hyssib(t)%snowf = 0
         hyssib(t)%rainf = 0
         hyssib(t)%evap = 0
         hyssib(t)%qs = 0
         hyssib(t)%qrec = 0
         hyssib(t)%qsb = 0
         hyssib(t)%qsm = 0
         hyssib(t)%qfz = 0
         hyssib(t)%qst = 0
         hyssib(t)%delsoilmoist = 0
         hyssib(t)%delswe = 0
         hyssib(t)%delintercept = 0
         hyssib(t)%snowt = 0
         hyssib(t)%vegtc = 0
         hyssib(t)%baresoilt = 0
         hyssib(t)%avgsurft = 0
         hyssib(t)%radteff = 0
         hyssib(t)%albedo = 0
         hyssib(t)%swe = 0
         hyssib(t)%sweveg = 0
         hyssib(t)%soilmoist1 = 0
         hyssib(t)%soilmoist2 = 0
         hyssib(t)%soilmoist3 = 0
         hyssib(t)%soiltemp = 0
         hyssib(t)%soilwet = 0
         hyssib(t)%potevap = 0
         hyssib(t)%ecanop = 0
         hyssib(t)%tveg = 0
         hyssib(t)%esoil = 0
         hyssib(t)%rootmoist = 0
         hyssib(t)%canopint = 0
         hyssib(t)%subsnow = 0
         hyssib(t)%subsurf = 0
         hyssib(t)%acond = 0
         hyssib(t)%snowfrac = 0
         hyssib(t)%snowdepth = 0
         hyssib(t)%sliqfrac = 0

         hyssib(t)%count = 0
         hyssib(t)%snowtcount = 0
         hyssib(t)%albedocount = 0
         hyssib(t)%sliqfraccount = 0
      enddo

#if ( ! defined OPENDAP )
      endif

      if ( npes > 1 ) then
         call hyssib_scatter
      endif
#endif
      return
!EOC
      end subroutine hyssib_setup

