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
#include "misc.h"
!BOP
!
! !ROUTINE: sethyssibp.f
!
! !DESCRIPTION:
!  This subroutine retrieves HY-SSiB parameters - Significant F90 revisions
!   below this subroutine will be required in the future.
!
! !REVISION HISTORY:
! 28 Apr 2002: Kristi Arsenault, Added SSiB LSM, Initial Code
! 13 Oct 2003: Sujay Kumar, Domain independent modifications
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      SUBROUTINE SETHYSSIBP(LD,LP,tile)
! !USES:
      use lis_module      ! LIS non-model-specific 1-D variables
      use hyssib_varder   ! HY-SSiB tile variables
      use tile_module
#if ( defined OPENDAP )
      use opendap_module
#endif
      implicit none

! !ARGUMENTS:
      type (lisdomain) LD
      type (lisparameters) LP
      type (tiledec) tile(LD%NCH)
!EOP
!=== Local Variables =====================================================
      INTEGER  :: ITYP,ICG,IWV,ILD,IDP,IBD
      PARAMETER (ITYP=13,ICG=2,IWV=3,ILD=2,IDP=3,IBD=2)
      REAL :: RSTPAR_v(ITYP,ICG,IWV), CHIL_v(ITYP,ICG), TOPT_v(ITYP,ICG), &
              TLL_v(ITYP,ICG), TU_v(ITYP,ICG), DEFAC_v(ITYP,ICG),         &
              PH1_v(ITYP,ICG), PH2_v(ITYP,ICG),  ROOTD_v(ITYP,ICG),       &
              BEE_v(ITYP), PHSAT_v(ITYP), SATCO_v(ITYP), POROS_v(ITYP),   &
              ZDEPTH_v(ITYP,IDP), SLOPE_v(ITYP)
      INTEGER :: N,I,J,K,JJ,c,r                  !Loop counters
      REAL :: VALUE(LP%NT,hyssibdrv%HYSSIB_NVEGP)

!=== End Variable Definition =============================================
!=== Parameters throught data

!=== Convert UMD Classes to SIB Classes for Each Tile
      print*,'MSG: sethyssibp -- Calling MAPVEGC to convert UMD to SIB', &
               ' (', iam,')'
      print*,'DBG: sethyssibp -- nch',ld%nch,' (',iam,')'
      print*,'DBG: sethyssibp -- size of hyssib',size(hyssib), &
               ' (',iam,')'

      DO N=1,ld%nch
         CALL HYSSIB_MAPVEGC(TILE(N)%VEGT)
         hyssib(N)%vegt = tile(n)%vegt
      ENDDO                     !N
      print*,'DBG: sethyssibp -- left MAPVEGC',' (',iam,')'

!=== Get Vegetation Parameters for HY-SSiB Model in Tile Space

!=== Read in the HY-SSiB Static Vegetation Parameter Files

      OPEN(UNIT=11,FILE=hyssibdrv%HYSSIB_VFILE,STATUS='OLD', &
           FORM='UNFORMATTED')
      READ(11) rstpar_v, chil_v, topt_v, tll_v, tu_v, defac_v,   &
               ph1_v, ph2_v, rootd_v, bee_v, phsat_v, satco_v,   &
               poros_v, zdepth_v, slope_v
      CLOSE(11)

!=== Assign STATIC vegetation parameters to each tile based on the
!=== type of vegetation present in that tile.
!=== These parameters will be stored in one long array--structured
!=== as follows: Tile 1, all the parameters (1 through numparam)
!=== then Tile 2, all the parameters. 
!=== Then Tile 3, all the parameters etc.

      DO I=1,ld%nch
         HYSSIB(I)%VEGP(1)=CHIL_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(2)=TOPT_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(3)=TLL_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(4)=TU_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(5)=DEFAC_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(6)=PH1_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(7)=PH2_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(8)=ROOTD_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(9)=ROOTD_v(TILE(I)%VEGT,2)
         HYSSIB(I)%VEGP(10)=RSTPAR_v(TILE(I)%VEGT,1,1)
         HYSSIB(I)%VEGP(11)=RSTPAR_v(TILE(I)%VEGT,1,2)
         HYSSIB(I)%VEGP(12)=RSTPAR_v(TILE(I)%VEGT,1,3)
         HYSSIB(I)%VEGP(13)=PHSAT_v(TILE(I)%VEGT)
         HYSSIB(I)%VEGP(14)=POROS_v(TILE(I)%VEGT)
         HYSSIB(I)%VEGP(15)=BEE_v(TILE(I)%VEGT)
         HYSSIB(I)%VEGP(16)=SATCO_v(TILE(I)%VEGT)
         HYSSIB(I)%VEGP(17)=SLOPE_v(TILE(I)%VEGT)
         HYSSIB(I)%VEGP(18)=ZDEPTH_v(TILE(I)%VEGT,1)
         HYSSIB(I)%VEGP(19)=ZDEPTH_v(TILE(I)%VEGT,2)
         HYSSIB(I)%VEGP(20)=ZDEPTH_v(TILE(I)%VEGT,3)
      ENDDO !I

      return
      end

