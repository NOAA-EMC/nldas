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
! !MODULE: hyssib_module.f:
!
! !DESCRIPTION:
!  Module for 1-D HY-SSiB land model driver variable specification
!
! !REVISION HISTORY:
! 28 Apr 2002: Kristi Arsenault, added SSiB LSM 2.5 code to LDAS
! 14 Nov 2002: Sujay Kumar, Optimized version for LIS
! 21 Apr 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      MODULE hyssib_module
!EOP
      implicit none

      public hyssibdec
      type hyssibdec
!BOC
      INTEGER :: ts                 !Timestep (seconds)
      INTEGER :: maxt               !Maximum tiles per grid
      INTEGER :: HYSSIBVEG          !UMD to HY-SSiB Vegetation Class Index value
      INTEGER :: COUNT
      INTEGER :: SNOWTCOUNT
      INTEGER :: ALBEDOCOUNT
      INTEGER :: SLIQFRACCOUNT
      INTEGER :: ZOBSOIL(1)  !Zobler Soil Classes (LIS%NCH)

      REAL :: VEGP(20)     !Static vegetation parameters,   dim(HY-SSiB_NVEGP)
      REAL :: VEGIP(11)    !Interpolated monthly parameters,dim(HY-SSiB_NVEGIP)
      REAL :: VEGT         !vegetation type of tile
      REAL :: ALBSF1       !Date 1 Snow-Free Albedo Value
      REAL :: ALBSF2       !Date 2 Snow-Free Albedo Value
      REAL :: ALBSF        !Quarterly Snow-Free Albedo dataset

!=== HY-SSiB-Forcing Variables ==============================================
      REAL :: FORCING(10)        ! TILE FORCING..
!=== HY-SSiB-Constant Variables ==============================================
      REAL :: TEMPBOT
      REAL :: TEMPSTD
!=== HY-SSiB-Initial Variables ==============================================
      REAL :: SGFG
      REAL :: SDENS
!=== HY-SSiB-State Variables ================================================
      REAL :: TC
      REAL :: TG
      REAL :: TSN
      REAL :: TD
      REAL :: WWW(3)
      REAL :: CAPAC(2)
      REAL :: SNOW(2)
!=== HY-SSiB-Output =========================================================
      REAL :: swnet
      REAL :: lwnet
      REAL :: qle
      REAL :: qh
      REAL :: qg
      REAL :: qf
      REAL :: qv
      REAL :: qtau
      REAL :: qa
      REAL :: delsurfheat
      REAL :: delcoldcont
      REAL :: snowf
      REAL :: rainf
      REAL :: evap
      REAL :: qs
      REAL :: qrec
      REAL :: qsb
      REAL :: qsm
      REAL :: qfz
      REAL :: qst
      REAL :: delsoilmoist
      REAL :: delswe
      REAL :: delintercept
      REAL :: snowt
      REAL :: vegtc
      REAL :: baresoilt
      REAL :: avgsurft
      REAL :: radteff
      REAL :: albedo
      REAL :: swe
      REAL :: sweveg
      REAL :: soilmoist1
      REAL :: soilmoist2
      REAL :: soilmoist3
      REAL :: soiltemp
      REAL :: soilwet
      REAL :: potevap
      REAL :: ecanop
      REAL :: tveg
      REAL :: esoil
      REAL :: rootmoist
      REAL :: canopint
      REAL :: subsnow
      REAL :: subsurf
      REAL :: acond
      REAL :: snowfrac
      REAL :: snowdepth
      REAL :: sliqfrac

      end type hyssibdec
!=== End Variable List ===================================================
!EOC
      end module hyssib_module

