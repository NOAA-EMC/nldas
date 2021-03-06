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
! !ROUTINE: setssibp.f
!
! !DESCRIPTION:
!  This subroutine retrieves SSiB parameters - Significant F90 revisions
!   below this subroutine will be required in the future.
!
! !REVISION HISTORY:
! 28 Apr 2002: Kristi Arsenault, Added SSiB LSM, Initial Code
! 13 Oct 2003: Sujay Kumar, Domain independent modifications
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
!  5 May 2004: David Mocko, made compatible with SiB-lings
!
! !INTERFACE:
      subroutine setssibp
! !USES:
      use lisdrv_module, only : grid,tile,lis
      use ssib_varder     ! SSiB tile variables
#if ( defined OPENDAP )
      use opendap_module
#endif
!EOP
      implicit none
!=== Local Variables =====================================================
      INTEGER  :: ITYP,ICG,IWV,ILD,IDP,IBD
      PARAMETER (ITYP=13,ICG=2,IWV=3,ILD=2,IDP=3,IBD=2)
      REAL      :: TRAN  (ITYP,ICG,IWV,ILD),REF   (ITYP,ICG,IWV,ILD),&
    		   RSTPAR(ITYP,ICG,IWV)    ,SOREF (ITYP,IWV),&
    		   ZCHIL  (ITYP,ICG),TOPT  (ITYP,ICG),TLL  (ITYP,ICG),&
    		   TU	 (ITYP,ICG),DEFAC (ITYP,ICG),PH1  (ITYP,ICG),&
    		   PH2   (ITYP,ICG),ROOTD (ITYP,ICG),RLMAX(ITYP,ICG),&
    		   ROOTCA(ITYP,ICG),RPLANT(ITYP,ICG),RDRES(ITYP,ICG),&
    		   BEE   (ITYP),PHSAT (ITYP),SATCO (ITYP),POROS(ITYP),&
    		   SLOPE (ITYP),ZDEPTH(ITYP,IDP)
      INTEGER :: N,I,J,K,JJ,c,r                  !Loop counters
!      REAL :: VALUE(LP%NT,ssibdrv%SSIB_NVEGP)

!=== End Variable Definition =============================================
!=== Parameters through data

      	DATA TRAN /  &
         0.5000000E-01,  0.5000000E-01,  0.5000000E-01,  0.5000000E-01,&
     	 0.5000000E-01,  0.5000000E-01,  0.7000000E-01,  0.5000000E-01,&
     	 0.5000000E-01,  0.5000000E-01,  0.1000000E-02,  0.5000000E-01,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.7000000E-01,  0.1000000E-02,  0.7000000E-01,&
     	 0.1000000E-02,  0.7000000E-01,  0.1000000E-02,  0.7000000E-01,&
     	 0.1000000E-02,&
     	 0.2500000E+00,  0.2500000E+00,  0.1500000E+00,  0.1000000E+00,&
     	 0.1000000E+00,  0.2500000E+00,  0.2475000E+00,  0.2500000E+00,&
     	 0.2500000E+00,  0.2500000E+00,  0.1000000E-02,  0.2500000E+00,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.2475000E+00,  0.1000000E-02,  0.2475000E+00,&
     	 0.1000000E-02,  0.2475000E+00,  0.1000000E-02,  0.2475000E+00,&
     	 0.1000000E-02,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.2200000E+00,  0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.2200000E+00,  0.1000000E-02,  0.2200000E+00,&
     	 0.1000000E-02,  0.2200000E+00,  0.1000000E-02,  0.2200000E+00,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.3750000E+00,  0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.3750000E+00,  0.1000000E-02,  0.3750000E+00,&
     	 0.1000000E-02,  0.3750000E+00,  0.1000000E-02,  0.3750000E+00,&
     	 0.1000000E-02,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00 /
      data ref / &
     	 0.1000000E+00,  0.1000000E+00,  0.7000000E-01,  0.7000000E-01,&
     	 0.7000000E-01,  0.1000000E+00,  0.1050000E+00,  0.1000000E+00,&
     	 0.1000000E+00,  0.1000000E+00,  0.1000000E-02,  0.1000000E+00,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.1050000E+00,  0.1000000E-02,  0.1050000E+00,&
     	 0.1000000E-02,  0.1050000E+00,  0.1000000E-02,  0.1050000E+00,&
     	 0.1000000E-02,&
     	 0.4500000E+00,  0.4500000E+00,  0.4000000E+00,  0.3500000E+00,&
     	 0.3500000E+00,  0.4500000E+00,  0.5775000E+00,  0.4500000E+00,&
     	 0.4500000E+00,  0.4500000E+00,  0.1000000E-02,  0.4500000E+00,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.5775000E+00,  0.1000000E-02,  0.5775000E+00,&
     	 0.1000000E-02,  0.5775000E+00,  0.1000000E-02,  0.5775000E+00,&
     	 0.1000000E-02,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,&
     	 0.1600000E+00,  0.1600000E+00,  0.1600000E+00,  0.1600000E+00,&
     	 0.1600000E+00,  0.1600000E+00,  0.3600000E+00,  0.1600000E+00,&
     	 0.1600000E+00,  0.1600000E+00,  0.1000000E-02,  0.1600000E+00,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.3600000E+00,  0.1000000E-02,  0.3600000E+00,&
     	 0.1000000E-02,  0.3600000E+00,  0.1000000E-02,  0.3600000E+00,&
     	 0.1000000E-02,&
     	 0.3900000E+00,  0.3900000E+00,  0.3900000E+00,  0.3900000E+00,&
     	 0.3900000E+00,  0.3900000E+00,  0.5775000E+00,  0.3900000E+00,&
     	 0.3900000E+00,  0.3900000E+00,  0.1000000E-02,  0.3900000E+00,&
     	 0.1000000E-02,&
     	 0.1000000E-02,  0.1000000E-02,  0.1000000E-02,  0.1000000E-02,&
     	 0.1000000E-02,  0.5775000E+00,  0.1000000E-02,  0.5775000E+00,&
     	 0.1000000E-02,  0.5775000E+00,  0.1000000E-02,  0.5775000E+00,&
     	 0.1000000E-02,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00/
     data rstpar / &
     	 0.2335900E+04,  0.9802230E+04,  0.6335955E+04,  0.2869680E+04,&
     	 0.2869680E+04,  0.5665121E+05,  0.2582010E+04,  0.9398942E+05,&
     	 0.9398942E+05,  0.9802230E+04,  0.1000000E+04,  0.7459000E+04,&
     	 0.1000000E+04,&
     	 0.2335900E+04,  0.9802230E+04,  0.6335955E+04,  0.2869680E+04,&
     	 0.2869680E+04,  0.2582010E+04,  0.2582010E+04,  0.2582010E+04,&
     	 0.1000000E+01,  0.2582010E+04,  0.1000000E+04,  0.7459000E+04,&
     	 0.1000000E+04,&
     	 0.1450000E-01,  0.1055000E+02,  0.7120000E+01,  0.3690000E+01,&
     	 0.3690000E+01,  0.1083000E+02,  0.1090000E+01,  0.1000000E-01,&
     	 0.1000000E-01,  0.1055000E+02,  0.1000000E+04,  0.5700000E+01,&
     	 0.1000000E+04,&
     	 0.1450000E-01,  0.1055000E+02,  0.7120000E+01,  0.3690000E+01,&
     	 0.3690000E+01,  0.1090000E+01,  0.1090000E+01,  0.1090000E+01,&
     	 0.1000000E+01,  0.1090000E+01,  0.1000000E+04,  0.5700000E+01,&
     	 0.1000000E+04,&
     	 0.1534900E+03,  0.1800000E+03,  0.2065000E+03,  0.2330000E+03,&
     	 0.2330000E+03,  0.1650000E+03,  0.1100000E+03,  0.8550000E+03,&
     	 0.8550000E+03,  0.1800000E+03,  0.1000000E+04,  0.2520000E+02,&
     	 0.1000000E+04,&
     	 0.1534900E+03,  0.1800000E+03,  0.2065000E+03,  0.2330000E+03,&
     	 0.2330000E+03,  0.1100000E+03,  0.1100000E+03,  0.1100000E+03,&
     	 0.1000000E+01,  0.1100000E+03,  0.1000000E+04,  0.2520000E+02,&
     	 0.1000000E+04/
     data soref / &
     	 0.1100000E+00,  0.1100000E+00,  0.1100000E+00,  0.1100000E+00,&
     	 0.1100000E+00,  0.1100000E+00,  0.1000000E+00,  0.1000000E+00,&
     	 0.3000000E+00,  0.1000000E+00,  0.3000000E+00,  0.1000000E+00,&
     	 0.1000000E+00,&
     	 0.2250000E+00,  0.2250000E+00,  0.2250000E+00,  0.2250000E+00,&
     	 0.2250000E+00,  0.2250000E+00,  0.2000000E+00,  0.2000000E+00,&
     	 0.3500000E+00,  0.2000000E+00,  0.3500000E+00,  0.1500000E+00,&
     	 0.1500000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00/
     data zchil / &
     	 0.1000000E+00,  0.2500000E+00,  0.1300000E+00,  0.1000000E-01,&
     	 0.1000000E-01,  0.1000000E-01, -0.3000000E+00,  0.1000000E-01,&
     	 0.1000000E-01,  0.2000000E+00,  0.1000000E-01, -0.2000000E-01,&
     	 0.1000000E-01,&
     	 0.1000000E+00,  0.2500000E+00,  0.1300000E+00,  0.1000000E-01,&
     	 0.1000000E-01, -0.3000000E+00, -0.3000000E+00, -0.3000000E+00,&
     	 0.1000000E-01,  0.2000000E+00,  0.1000000E-01, -0.2000000E-01,&
     	 0.1000000E-01/
     data topt / &
     	 0.3030000E+03,  0.3000000E+03,  0.2940000E+03,  0.2880000E+03,&
     	 0.2880000E+03,  0.2970000E+03,  0.3130000E+03,  0.3150000E+03,&
     	 0.3150000E+03,  0.3000000E+03,  0.3100000E+03,  0.3000000E+03,&
     	 0.3100000E+03,&
     	 0.3030000E+03,  0.3000000E+03,  0.2940000E+03,  0.2880000E+03,&
     	 0.2880000E+03,  0.3120000E+03,  0.3130000E+03,  0.3130000E+03,&
     	 0.3150000E+03,  0.2890000E+03,  0.3100000E+03,  0.3000000E+03,&
     	 0.3100000E+03/
     data tll / &
     	 0.2730000E+03,  0.2730000E+03,  0.2700000E+03,  0.2680000E+03,&
     	 0.2680000E+03,  0.2730000E+03,  0.2830000E+03,  0.2830000E+03,&
     	 0.2830000E+03,  0.2730000E+03,  0.3000000E+03,  0.2730000E+03,&
     	 0.3000000E+03,&
     	 0.2730000E+03,  0.2730000E+03,  0.2700000E+03,  0.2680000E+03,&
     	 0.2680000E+03,  0.2730000E+03,  0.2830000E+03,  0.2830000E+03,&
     	 0.2830000E+03,  0.2730000E+03,  0.3000000E+03,  0.2730000E+03,&
     	 0.3000000E+03/
     data tu / &
     	 0.3180000E+03,  0.3180000E+03,  0.3150000E+03,  0.3130000E+03,&
     	 0.3130000E+03,  0.3230000E+03,  0.3280000E+03,  0.3230000E+03,&
     	 0.3230000E+03,  0.3230000E+03,  0.3200000E+03,  0.3180000E+03,&
     	 0.3200000E+03,&
     	 0.3180000E+03,  0.3180000E+03,  0.3150000E+03,  0.3130000E+03,&
     	 0.3130000E+03,  0.3230000E+03,  0.3280000E+03,  0.3280000E+03,&
     	 0.3230000E+03,  0.3090000E+03,  0.3200000E+03,  0.3150000E+03,&
     	 0.3200000E+03/
     data defac / &
     	 0.2730000E-01,  0.3570000E-01,  0.3400000E-01,  0.3100000E-01,&
     	 0.3100000E-01,  0.3570000E-01,  0.2380000E-01,  0.2750000E-01,&
     	 0.2750000E-01,  0.2750000E-01,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00,&
     	 0.2730000E-01,  0.3570000E-01,  0.3400000E-01,  0.3100000E-01,&
     	 0.3100000E-01,  0.2380000E-01,  0.2380000E-01,  0.2380000E-01,&
     	 0.2380000E-01,  0.2380000E-01,  0.0000000E+00,  0.0000000E+00,&
     	 0.0000000E+00/
     data ph1 / &
     	 0.1200000E+01,  0.5350000E+01,  0.1920000E+01,  0.3700000E+01,&
     	 0.7800000E+01,  0.1800000E+01,  0.1730000E+01,  0.1920000E+01,&
     	 0.1390000E+01,  0.9600000E+00,  0.3000000E+01,  0.1800000E+01,&
     	 0.5000000E+01,&
     	 0.1200000E+01,  0.5350000E+01,  0.1920000E+01,  0.3700000E+01,&
     	 0.7800000E+01,  0.1800000E+01,  0.1730000E+01,  0.1920000E+01,&
     	 0.1390000E+01,  0.9600000E+00,  0.3000000E+01,  0.1800000E+01,&
     	 0.5000000E+01/
     data ph2 / &
     	 0.6250000E+01,  0.5570000E+01,  0.5730000E+01,  0.5530000E+01,&
     	 0.5660000E+01,  0.5670000E+01,  0.5800000E+01,  0.5610000E+01,&
     	 0.6370000E+01,  0.5370000E+01,  0.6000000E+01,  0.5670000E+01,&
     	 0.6000000E+01,&
     	 0.6250000E+01,  0.5570000E+01,  0.5730000E+01,  0.5530000E+01,&
     	 0.5660000E+01,  0.5670000E+01,  0.5800000E+01,  0.5610000E+01,&
     	 0.6370000E+01,  0.5370000E+01,  0.6000000E+01,  0.5670000E+01,&
     	 0.6000000E+01/
     data rootd / &
     	 0.1000000E+01,  0.1000000E+01,  0.1000000E+01,  0.5000000E+00,&
     	 0.5000000E+00,  0.5000000E+00,  0.5000000E+00,  0.5000000E+00,&
     	 0.5000000E+00,  0.2000000E+00,  0.1000000E+00,  0.1000000E+01,&
     	 0.1000000E+01,&
     	 0.1000000E+01,  0.1000000E+01,  0.1000000E+01,  0.5000000E+00,&
     	 0.5000000E+00,  0.5000000E+00,  0.5000000E+00,  0.5000000E+00,&
     	 0.5000000E+00,  0.2000000E+00,  0.1000000E+00,  0.1000000E+01,&
     	 0.1000000E+01/
     data bee / &
     	 0.7120000E+01,  0.7120000E+01,  0.7120000E+01,  0.7120000E+01,&
     	 0.7120000E+01,  0.7120000E+01,  0.7120000E+01,  0.4050000E+01,&
     	 0.4050000E+01,  0.7120000E+01,  0.4050000E+01,  0.7797000E+01,&
     	 0.4804000E+01/
     data phsat / &
     	-0.8600000E-01, -0.8600000E-01, -0.8600000E-01, -0.8600000E-01,&
     	-0.8600000E-01, -0.8600000E-01, -0.8600000E-01, -0.3500000E-01,&
     	-0.3500000E-01, -0.8600000E-01, -0.3500000E-01, -0.1980000E+00,&
     	-0.1670000E+00/
     data satco / &
     	 0.2000000E-04,  0.2000000E-04,  0.2000000E-04,  0.2000000E-04,&
     	 0.2000000E-04,  0.2000000E-04,  0.2000000E-04,  0.1760000E-03,&
     	 0.1760000E-03,  0.2000000E-04,  0.1760000E-03,  0.3500000E-05,&
     	 0.7620000E-04/
     data poros / &
     	 0.4200000E+00,  0.4200000E+00,  0.4200000E+00,  0.4200000E+00,&
     	 0.4200000E+00,  0.4200000E+00,  0.4200000E+00,  0.4352000E+00,&
     	 0.4352000E+00,  0.4200000E+00,  0.4352000E+00,  0.4577000E+00,&
     	 0.4352000E+00/
     data slope / &
     	 0.1736000E+00,  0.1736000E+00,  0.1736000E+00,  0.1736000E+00,&
     	 0.1736000E+00,  0.1736000E+00,  0.1736000E+00,  0.8720000E-01,&
     	 0.8720000E-01,  0.1736000E+00,  0.8720000E-01,  0.3420000E+00,&
     	 0.8720000E-01/
     data zdepth / &
     	 0.2000000E-01,  0.2000000E-01,  0.2000000E-01,  0.2000000E-01,&
     	 0.2000000E-01,  0.2000000E-01,  0.2000000E-01,  0.2000000E-01,&
     	 0.2000000E-01,  0.2000000E-01,  0.2000000E-01,  0.2000000E-01,&
     	 0.1000000E+01,&
     	 0.1480000E+01,  0.1480000E+01,  0.1480000E+01,  0.1480000E+01,&
     	 0.1480000E+01,  0.1480000E+01,  0.4700000E+00,  0.4700000E+00,&
     	 0.4700000E+00,  0.1700000E+00,  0.1700000E+00,  0.1480000E+01,&
     	 0.1000000E+01,&
     	 0.2000000E+01,  0.2000000E+01,  0.2000000E+01,  0.2000000E+01,&
     	 0.2000000E+01,  0.2000000E+01,  0.1000000E+01,  0.1000000E+01,&
     	 0.1000000E+01,  0.1000000E+01,  0.3000000E+00,  0.2000000E+01,&
     	 0.1000000E+01/

!=== Convert UMD Classes to SIB Classes for Each Tile
      print*,'MSG: setssibp -- Calling MAPVEGC to convert UMD to SIB', & 
               ' (', iam,')'
      print*,'DBG: setssibp -- nch',lis%d%nch,' (',iam,')'
      print*,'DBG: setssibp -- size of ssib',size(ssib), &
               ' (',iam,')'

      DO N=1,lis%d%nch
         CALL SSIB_MAPVEGC(TILE(N)%VEGT)
         ssib(N)%vegt = tile(n)%vegt
      ENDDO                     !N
       print*,'DBG: setssibp -- left MAPVEGC',' (',iam,')'

!=== Get Vegetation Parameters for SSiB Model in Tile Space

!=== Read in the SSiB Static Vegetation Parameter Files

!      OPEN(UNIT=11,FILE=ssibdrv%SSIB_VFILE,STATUS='OLD')

!      DO J=1,ssibdrv%SSIB_NVEGP
!        READ(11,*)(VALUE(I,J),I=1,LP%NT)
!      ENDDO 
!      CLOSE(11)

!=== Assign STATIC vegetation parameters to each tile based on the
!=== type of vegetation present in that tile.
!=== These parameters will be stored in one long array--structured
!=== as follows: Tile 1, all the parameters (1 through numparam)
!=== then Tile 2, all the parameters. 
!=== Then Tile 3, all the parameters etc.

      DO I=1,lis%d%nch
         SSIB(I)%VEGP(1)=TRAN(TILE(I)%VEGT,1,1,1)
         SSIB(I)%VEGP(2)=TRAN(TILE(I)%VEGT,2,1,1)
         SSIB(I)%VEGP(3)=TRAN(TILE(I)%VEGT,1,2,1)
         SSIB(I)%VEGP(4)=TRAN(TILE(I)%VEGT,2,2,1)
         SSIB(I)%VEGP(5)=TRAN(TILE(I)%VEGT,1,3,1)
         SSIB(I)%VEGP(6)=TRAN(TILE(I)%VEGT,2,3,1)
         SSIB(I)%VEGP(7)=TRAN(TILE(I)%VEGT,1,1,2)
         SSIB(I)%VEGP(8)=TRAN(TILE(I)%VEGT,2,1,2)
         SSIB(I)%VEGP(9)=TRAN(TILE(I)%VEGT,1,2,2)
         SSIB(I)%VEGP(10)=TRAN(TILE(I)%VEGT,2,2,2)
         SSIB(I)%VEGP(11)=TRAN(TILE(I)%VEGT,1,3,2)
         SSIB(I)%VEGP(12)=TRAN(TILE(I)%VEGT,2,3,2)
         SSIB(I)%VEGP(13)=REF(TILE(I)%VEGT,1,1,1)
         SSIB(I)%VEGP(14)=REF(TILE(I)%VEGT,2,1,1)
         SSIB(I)%VEGP(15)=REF(TILE(I)%VEGT,1,2,1)
         SSIB(I)%VEGP(16)=REF(TILE(I)%VEGT,2,2,1)
         SSIB(I)%VEGP(17)=REF(TILE(I)%VEGT,1,3,1)
         SSIB(I)%VEGP(18)=REF(TILE(I)%VEGT,2,3,1)
         SSIB(I)%VEGP(19)=REF(TILE(I)%VEGT,1,1,2)
         SSIB(I)%VEGP(20)=REF(TILE(I)%VEGT,2,1,2)
         SSIB(I)%VEGP(21)=REF(TILE(I)%VEGT,1,2,2)
         SSIB(I)%VEGP(22)=REF(TILE(I)%VEGT,2,2,2)
         SSIB(I)%VEGP(23)=REF(TILE(I)%VEGT,1,3,2)
         SSIB(I)%VEGP(24)=REF(TILE(I)%VEGT,2,3,2)
         SSIB(I)%VEGP(25)=RSTPAR(TILE(I)%VEGT,1,1)
         SSIB(I)%VEGP(26)=RSTPAR(TILE(I)%VEGT,1,2)
         SSIB(I)%VEGP(27)=RSTPAR(TILE(I)%VEGT,1,3)
         SSIB(I)%VEGP(28)=RSTPAR(TILE(I)%VEGT,2,1)
         SSIB(I)%VEGP(29)=RSTPAR(TILE(I)%VEGT,2,2)
         SSIB(I)%VEGP(30)=RSTPAR(TILE(I)%VEGT,2,3)
         SSIB(I)%VEGP(31)=SOREF(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(32)=SOREF(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(33)=SOREF(TILE(I)%VEGT,3)
         SSIB(I)%VEGP(34)=ZCHIL(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(35)=ZCHIL(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(36)=TOPT(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(37)=TOPT(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(38)=TLL(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(39)=TLL(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(40)=TU(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(41)=TU(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(42)=DEFAC(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(43)=DEFAC(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(44)=PH1(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(45)=PH1(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(46)=PH2(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(47)=PH2(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(48)=ROOTD(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(49)=ROOTD(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(50)=BEE(TILE(I)%VEGT)
         SSIB(I)%VEGP(51)=PHSAT(TILE(I)%VEGT)
         SSIB(I)%VEGP(52)=SATCO(TILE(I)%VEGT)
         SSIB(I)%VEGP(53)=POROS(TILE(I)%VEGT)
         SSIB(I)%VEGP(54)=SLOPE(TILE(I)%VEGT)
         SSIB(I)%VEGP(55)=ZDEPTH(TILE(I)%VEGT,1)
         SSIB(I)%VEGP(56)=ZDEPTH(TILE(I)%VEGT,2)
         SSIB(I)%VEGP(57)=ZDEPTH(TILE(I)%VEGT,3)
         SSIB(I)%ZWINDINI   = 58.0
         SSIB(I)%ZMETINI    = 58.0
         SSIB(I)%ITRUNKINI  = 20
         SSIB(I)%ILWINI     = 1
      ENDDO !I
      return
!EOC
      end subroutine setssibp

