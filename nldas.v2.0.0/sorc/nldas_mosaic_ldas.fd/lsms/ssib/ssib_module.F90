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
! !MODULE: ssib_module.f:
!
! !DESCRIPTION:
!  Module for 1-D SSiB land model driver variable specification
!
! !REVISION HISTORY:
! 28 Apr 2002: Kristi Arsenault, added SSiB LSM 2.5 code to LDAS
! 14 Nov 2002: Sujay Kumar, Optimized version for LIS
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
! 22 May 2004: David Mocko, made compatible with SiB-lings
! 
! !INTERFACE:
      MODULE ssib_module
!EOP
      implicit none

      public ssibdec
      type ssibdec
!BOC
      INTEGER :: ts                 !Timestep (seconds)
      INTEGER :: maxt               !Maximum tiles per grid
      INTEGER :: SIBVEG             !UMD to SiB Vegetation Class Index value
      INTEGER :: NSLAY              !Number of SSiB soil layers (4)
      INTEGER :: COUNT
      INTEGER :: ALBEDOCOUNT
      INTEGER :: ZOBSOIL(1)  !Zobler Soil Classes (LIS%NCH)
      INTEGER :: ITRUNKINI
      INTEGER :: ILWINI

      REAL :: VEGP(57)     !Static vegetation parameters,   dim(SSIB_NVEGP)
      REAL :: VEGIP(12)    !Interpolated monthly parameters,dim(SSIB_NVEGIP)
      REAL :: VEGT         !vegetation type of tile
      REAL :: XDTT

!=== SSiB-Forcing Variables ============================================
      REAL :: FORCING(10)        ! TILE FORCING..
!=== SSiB-Initial Variables ============================================
      REAL :: ZWINDINI
      REAL :: ZMETINI
!=== SSiB-State Variables ==============================================
      REAL :: TCINI
      REAL :: TGSINI
      REAL :: TDINI
      REAL :: TAINI
      REAL :: TMINI
      REAL :: HTINI
      REAL :: QAINI
      REAL :: WWWINI(3)
      REAL :: CAPACINI(2)
!==== SSiB-Output ======================================================
      REAL :: swnet
      REAL :: lwnet
      REAL :: qle
      REAL :: qh
      REAL :: qg
      REAL :: qf
      REAL :: qtau
      REAL :: delsurfheat
      REAL :: snowf
      REAL :: rainf
      REAL :: evap
      REAL :: qs
      REAL :: qsb
      REAL :: qsm
      REAL :: delsoilmoist
      REAL :: delswe
      REAL :: delintercept
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
      REAL :: ecanop
      REAL :: tveg
      REAL :: esoil
      REAL :: rootmoist
      REAL :: canopint
      REAL :: acond
      REAL :: snowfrac

      end type ssibdec
!=== End Variable List =================================================
!EOC
      end module ssib_module

