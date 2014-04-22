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
! !ROUTINE: ssib_com.F90
!
! !DESCRIPTION:
!  Common block internal to the SSiB main program
!
! !REVISION HISTORY:
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
! 22 May 2004: David Mocko, Added some variables to common block
!
! !INTERFACE:
	module ssib_common
!EOP
!===Prognostic variables	
	REAL :: TC, TS, TGS, TD, CAPAC(2), WWW(3)
!===Constants
	REAL :: PIE, TIMCON, CPAIR, RHOAIR, PSY, HLAT, GRAV, VKC,&
                SNOMEL, STEFAN, TF, CLAI, CW ,XTEM1, XTEM2,&
		GASR
!===VEGETATION AND SOIL PARAMETERS
        REAL :: TRAN(2,3,2), REF(2,3,2), GREEN(2), VCOVER(2), &
	        CHIL(2) , RSTPAR(2,3), TOPT(2), TLL(2),       &       
                TU(2), DEFAC(2), PH1(2), PH2(2),              & 
    		ZLT(2), Z0, D, Z2, Z1, RDC, RBC	,             & 	     
    		RPLANT(2), ROOTD(2), ROOTL(2), RDRES(2),      &    
                rootp(3),ROOTCA(2), RLMAX(2),SOREF(3),        & 
     		BEE, PHSAT, POROS, SATCO, SLOPE,ZDEPTH(3)	     
!===INPUT ATMOSPHERIC AND SITE DATA 
	REAL ::  EM, TM, UM, ZM, PSUR, PPC, PPL, RADN(3,2),   &
                 SUNANG, SWDOWN, RNETM, CLOUD, swup,QM,       &
                 totrad ,pr,us
	INTEGER :: NROOT	 
		 
!===SITE PARAMETERS SPECIFIC TO 1-D MODEL OPERATION 
        REAL ::  G1, G2, G3, ZTZ0, CORB1, CORB2, HA, HT,      &    
                 ZWIND, ZMET, ZLONG, ZLAT, SALB(2,2), RAB(2,3,2)                             
!=== VARIABLES RETURNED TO G.C.M. ( ET(KG), H(W M-2), RUNOFF(M),  
!                                   LW(W M-2), DRAG(KG M-1 S-2) )      
        REAL ::  ETMASS, HFLUX, ROFF, ZLWUP, DRAG		     
!===	VARIABLES CALCULATED FROM ABOVE AND AMBIENT CONDITIONS	  
        REAL ::  RA, RB, RD, RCC, RG,	                      &	      		 
        	 TA, EA, ETC, ETGS, GETC, GETGS ,U2, USTAR,   &	 
        	 ALBEDO(2,3,2), EXTK(2,3,2), RADT(2), THERMK, &	 
        	 RADFAC(2,2,2), RADSAV(12),		      & 
        	 PAR(2), PD(2), RST(2), RSTFAC(2,4), DROP,    &
		 RSOIL, PHL(2), PHROOT(2), ROOTR(2),	      &     		 
        	 PHSOIL(3), SATCAP(2), SNOWW(2),     	      &	 
        	 CCX, CG,DTC, DTG			     		 
!===HEAT FLUXES : C-CANOPY, G-GROUND, T-TRANS, E-EVAP  IN J M-2            
        REAL ::  ECT, ECI, EGT, EGI, EGS, HC, HG, CHF, SHF,   &     
                 ECMASS, EGMASS  
!------------- START COMSCROPS ------------------------
      REAL YSTPAR(2,3), YOPT  (2),YLL  (2),&
                    YU  (2),YEFAC (2),YH1 (2), YH2 (2)
!------------- END OF COMSCROPS ------------------------
!-------whatever is missing on TEMRS1-------------------
      REAL :: BPS,XXT,XXQ,XAKMS
!--------------END TEMRS1--------------------------------      
!===OTHERS
	REAL :: ANGULOSOLAR,DUMMYG(2), UWIND
	REAL :: VWIND, SFCSPD
	REAL :: w1i,w2i,w3i,tprec,cprec,ptot,difrat,vnrat
	REAL :: FRAC(2,2),RSTUN,TOTWB,ENDWB,ERROR,CBAL,GBAL
	REAL :: ZLHS,ZRHS,ectw,eciw,egtw,egiw,egsw
	REAL :: swcan,swgnd,DTT,smelt,q3g,tgeff,scov2

	INTEGER :: T,ilw,itrunk
	REAL :: QA,xadj,D1
        REAL ::  ZINC(3), A2(3), Y1(3)                            
        INTEGER ::  ITER(3)                                            
!===OLD ISTAT related variables (PILPS type experiment; see original ssib)
       REAL :: RNOFFS
       REAL :: SNM
       REAL :: SIBSU
       REAL :: BEDO
       REAL :: SOILDIF
       REAL :: SOILDRA
       REAL :: RNOFFT
       REAL :: RNOFFB
!=== End Variable List ===================================================
      end module ssib_common
	
