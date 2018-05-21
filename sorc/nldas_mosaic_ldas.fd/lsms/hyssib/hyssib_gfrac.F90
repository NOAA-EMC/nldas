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
! !ROUTINE: hyssib_gfrac.F90 
!
! !DESCRIPTION:
!  This subroutine takes vegetation greenness fraction data and the date
!  to interpolate and determine the actual value of the greenness fraction
!  for that date.  This actual value is then returned to the main program.
!  The assumption is that the data point is valid for the 16th of the given
!  month, at 00Z.
!
!  Gustavo will take advantage of this routine to read the other time varying
!  parameters and interpolate them to the appropriate time frame
!
! !REVISION HISTORY:
! 28 Apr 2002: Kristi Arsenault, Added SSiB LSM to LDAS, initial code
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      SUBROUTINE HYSSIB_GFRAC(LD,LT,TILE)
! !USES:
      USE lis_module      ! LIS non-model-specific 1-D variables
      USE hyssib_varder   ! HY-SSIB tile variables
      use time_manager
      use tile_module
#if ( defined OPENDAP )
      use opendap_module
#endif
!EOP
      implicit none

      type (lisdomain) LD
      type (listime) LT
      type (tiledec) :: tile(ld%nch)
!=== Local Variables =====================================================
      INTEGER :: I,T,C,R              ! Loop counters
      INTEGER :: ITYP,IMON,ICG,IWV,ILD,IDP,IBD
!== Gustavo added for SSiB
      PARAMETER (ITYP=13,IMON=12,ICG=2,IWV=3,ILD=2,IDP=3,IBD=2)
      REAL :: RSTPAR_v(ITYP,ICG,IWV), CHIL_v(ITYP,ICG), TOPT_v(ITYP,ICG), &
              TLL_v(ITYP,ICG), TU_v(ITYP,ICG), DEFAC_v(ITYP,ICG),         &
              PH1_v(ITYP,ICG), PH2_v(ITYP,ICG),  ROOTD_v(ITYP,ICG),       &
              BEE_v(ITYP), PHSAT_v(ITYP), SATCO_v(ITYP), POROS_v(ITYP),   &
              ZDEPTH_v(ITYP,IDP), SLOPE_v(ITYP)
      REAL :: GREEN_v(ITYP,IMON,ICG), VCOVER_v(ITYP,IMON,ICG),           &
              ZLT_v(ITYP,IMON,ICG), Z0_v(ITYP,IMON), DD_v(ITYP,IMON),    &
              Z2_v(ITYP,IMON), Z1_v(ITYP,IMON), RDC_v(ITYP,IMON),        &
              RBC_v(ITYP,IMON)

      OPEN(UNIT=11,FILE=hyssibdrv%HYSSIB_VFILE,STATUS='OLD', &
           FORM='UNFORMATTED')
      READ(11) rstpar_v, chil_v, topt_v, tll_v, tu_v, defac_v,   &
               ph1_v, ph2_v, rootd_v, bee_v, phsat_v, satco_v,   &
               poros_v, zdepth_v, slope_v
      READ(11) green_v, vcover_v, zlt_v, z0_v, dd_v, z2_v, z1_v,       &
               rdc_v, rbc_v
      CLOSE(11)

!=== End Variable Definition =============================================

      do i=1,ld%nch
         HYSSIB(I)%VEGIP(1) = Z0_v(TILE(I)%VEGT,LT%MO)
         HYSSIB(I)%VEGIP(2) = Z1_v(TILE(I)%VEGT,LT%MO)
         HYSSIB(I)%VEGIP(3) = Z2_v(TILE(I)%VEGT,LT%MO)
         HYSSIB(I)%VEGIP(4) = DD_v(TILE(I)%VEGT,LT%MO)
         HYSSIB(I)%VEGIP(5) = VCOVER_v(TILE(I)%VEGT,LT%MO,1)
         HYSSIB(I)%VEGIP(6) = VCOVER_v(TILE(I)%VEGT,LT%MO,2)
         HYSSIB(I)%VEGIP(7) = ZLT_v(TILE(I)%VEGT,LT%MO,1)
         HYSSIB(I)%VEGIP(8) = ZLT_v(TILE(I)%VEGT,LT%MO,2)
         HYSSIB(I)%VEGIP(9) = GREEN_v(TILE(I)%VEGT,LT%MO,1)
         HYSSIB(I)%VEGIP(10) = RBC_v(TILE(I)%VEGT,LT%MO)
         HYSSIB(I)%VEGIP(11) = RDC_v(TILE(I)%VEGT,LT%MO)
      enddo

      return
      end

