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
! !MODULE: mos_module.F90:
!
! !DESCRIPTION:
!  Module for 1-D MOSAIC land model driver variable specification.
!
! !REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial code
!  11 Feb 2002: Jon Gottschalck; Added AVHRR derived variables
!
! !INTERFACE:
module mos_module
!EOP
  implicit none
  public mosdec
!BOC
  type mosdec

      INTEGER :: ts                 !Timestep (seconds)      
      INTEGER :: NSLAY               !Number of MOSAIC soil layers
      INTEGER :: COUNT                !MOSAIC Output Counter
      INTEGER :: VEGT

      REAL :: VEGP(24)       !Static vegetation parameter values, dim(MOS_NVEGP)
!      REAL :: VEGMP(6,12)    !Monthly vegetation parameter values dim(MOS_NMVEGP,12)
      REAL :: VEGIP(6)      !Interpolated from monthly parameters (MOS_NMVEGP)
      REAL :: SOILP(10)      !Static soil parameter values, dim(NOS_NSOILP)

      REAL :: LAI           !AVHRR Derived LAI
      REAL :: GREEN         !AVHRR Derived GREENNESS
      REAL :: DSAI          !AVHRR Derived DSAI

!=== LDAS-MOSAIC States ===============================================
      REAL :: CT                     !MOSAIC Canopy/Soil Temperature 
      REAL :: QA                     !MOSAIC Canopy Humidity
      REAL :: ICS                    !MOSAIC Interception Canopy Storage
      REAL :: SNOW                   !MOSAIC Snow Depth
      REAL :: SoT                    !MOSAIC Deep Soil Temperaure
      REAL :: SoWet(3)      !MOSAIC Soil Wetness (3 layers)
!=== Analysis and Bias Correction Variables ===========================
      REAL :: DTCANAL                !MOSAIC Change in Temperature based on Analysis
      real :: lai1, lai2
      real :: sai1, sai2
!=== LDAS-MOSAIC OUTPUT States ============================================
      real :: swnet
      real :: lwnet
      real :: qle
      real :: qh
      real :: qg
      real :: swrad
      real :: lwrad
      real :: snowf
      real :: rainf
      real :: evap
      real :: qs
      real :: qsb
      real :: qsm
      real :: avgsurft
      real :: albedo
      real :: swe
      real :: soilmoist1
      real :: soilmoist2
      real :: soilmoist3
      real :: soilwet
      real :: ecanop
      real :: tveg
      real :: esoil 
      real :: rootmoist
      real :: canopint
      real :: acond
      real :: laiout	
      real :: soilm_prev
      real :: swe_prev
      real :: water1
      real :: water2
      real :: water3
      real :: snohf
      real :: sbsno
	real :: snod
	real :: snoc
	real :: mstavr
	real :: ccond
	real :: veg
	real :: soilv1
	real :: soilv2
	real :: soilv3
	real :: soilv4
	real :: soilv5
	real :: soilv6
	real :: soilv7
	real :: soilv8
	real :: soilmtot
	real :: soilm1
	real :: soilmr
        real :: snwfrcout
      real :: forcing(10)
!      REAL :: AC
!      REAL :: CC
      REAL :: LAT, LON
      end type mosdec
!EOC      
      end module mos_module












