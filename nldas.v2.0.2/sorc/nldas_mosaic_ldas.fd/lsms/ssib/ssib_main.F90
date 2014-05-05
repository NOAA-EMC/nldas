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
! !ROUTINE: ssib_main.f
!
! !DESCRIPTION:
! SSIB LAND-SURFACE MODEL, UNCOUPLED 1-D COLUMN: VERSION 2.5 OCT 2001
!
! !REVISION HISTORY:
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
! 25 May 2004: David Mocko, removed duplicate code and provided ALMA outputs
!
! PROGRAM HISTORY LOG (prior to LIS)
! VERSION 1.0  --  01 MAR 1999
! VERSION 1.1  --  08 MAR 1999
! VERSION 2.0  --  27 JUL 1999
! VERSION 2.1  --  23 OCT 2000
! VERSION 2.2  --  07 FEB 2001
! VERSION 2.3  --  07 MAY 2001 = operational Eta implementation
! VERSION 2.4  --  27 JUN 2001 = ops Eta with NO physics changes
! VERSION 2.5  --  18 OCT 2001
!  lis VERSION --  28 APR 2002 = SSIB Main added to lis Driver
!                                (NASA GSFC)
! VERSION 2.5.1--  28 MAY 2002 = Updated changes in SSIB LSM along
!                                with correction to SOILRZ and SOIL1M.
! VERSION 2.6  --  24 JUN 2003 = Updated to SSiB LSM v2.6
!
!   Physics changes:
!     in SUBROUTINE SFLX change CSOIL from 1.26E+6 to 2.00E+6
!     in SUBROUTINE SFLX change ZBOT from -3.0 to -8.0
!     Replaced de-bugged SUBROUTINE TDFCND
!     Removed SUBROUTINE REDPRM and moved the parameters to other
!      locations throughout ssib$_-$main and ssib$_-$physics subroutines.
!    VERSION 2.5.2 --  31 MAY 2002
!      Fix in FUNCTION DEVAP related to FX calculation
!    VERSION 2.6   --  Includes changes to certain parameters and
!                     snow-soil physics
! !INTERFACE:
      subroutine ssib_main
! !USES:
      use lisdrv_module, only : lis,grid,tile
      use ssib_varder      ! SSiB tile variables
      use tile_spmdMod
      use ssib_common
!EOP
      implicit none

! Define all local variables according to comsib
!      real :: ttii

! Mocko parameter dimension
      real :: latco, lonco, currwb, suralbedo, soilwet
      real :: startsm, startint, startswe
      integer :: iyearh, imonthh, idayh, isech
      logical prin

! gustavo parameter dimension
      INTEGER, PARAMETER:: NSOLD=20   ! Maximum number of soil layers
! crops section (define dummy for month)
      real :: dum,nsx,tl
! end crops section
!BOC

!=== Convert lis Timestep varname to SSIB timestep varname (DT) (sec)
      do t = 1, di_array(iam)
         DTT = float(lis%t%TS)

!=== RESET SOME LOCAL ACCUMULATIVE VARS
         RNOFFS = 0.
         SNM = 0.
         SIBSU = 0.
         BEDO = 0
         SOILDIF = 0
         SOILDRA = 0
         RNOFFT = 0
         RNOFFB = 0
         ITER = 1
!===END RESET SOME LOCAL ACCUMULATIVE VARS

         latco = grid(tile(t)%index)%lat
         lonco = grid(tile(t)%index)%lon
         prin = .false.
!         if ((latco.eq.66.5).and.(lonco.eq.-179.5)) prin = .true.

         CALL CONSTS

!==== At this point start to assign SSIB parameters to the respective variables

!***********************************************
!========static parameters
         TRAN(1,1,1) = SSIB(T)%VEGP(1)
         TRAN(2,1,1) = SSIB(T)%VEGP(2)
         TRAN(1,2,1) = SSIB(T)%VEGP(3)
         TRAN(2,2,1) = SSIB(T)%VEGP(4)
         TRAN(1,3,1) = SSIB(T)%VEGP(5)
         TRAN(2,3,1) = SSIB(T)%VEGP(6)
         TRAN(1,1,2) = SSIB(T)%VEGP(7)
         TRAN(2,1,2) = SSIB(T)%VEGP(8)
         TRAN(1,2,2) = SSIB(T)%VEGP(9)
         TRAN(2,2,2) = SSIB(T)%VEGP(10)
         TRAN(1,3,2) = SSIB(T)%VEGP(11)
         TRAN(2,3,2) = SSIB(T)%VEGP(12)
         REF(1,1,1) = SSIB(T)%VEGP(13)
         REF(2,1,1) = SSIB(T)%VEGP(14)
         REF(1,2,1) = SSIB(T)%VEGP(15)
         REF(2,2,1) = SSIB(T)%VEGP(16)
         REF(1,3,1) = SSIB(T)%VEGP(17)
         REF(2,3,1) = SSIB(T)%VEGP(18)
         REF(1,1,2) = SSIB(T)%VEGP(19)
         REF(2,1,2) = SSIB(T)%VEGP(20)
         REF(1,2,2) = SSIB(T)%VEGP(21)
         REF(2,2,2) = SSIB(T)%VEGP(22)
         REF(1,3,2) = SSIB(T)%VEGP(23)
         REF(2,3,2) = SSIB(T)%VEGP(24)
         RSTPAR(1,1) = SSIB(T)%VEGP(25)
         RSTPAR(1,2) = SSIB(T)%VEGP(26)
         RSTPAR(1,3) = SSIB(T)%VEGP(27)
         RSTPAR(2,1) = SSIB(T)%VEGP(28)
         RSTPAR(2,2) = SSIB(T)%VEGP(29)
         RSTPAR(2,3) = SSIB(T)%VEGP(30)
         SOREF(1) = SSIB(T)%VEGP(31)
         SOREF(2) = SSIB(T)%VEGP(32)
         SOREF(3) = SSIB(T)%VEGP(33)
         CHIL(1) = SSIB(T)%VEGP(34)
         CHIL(2) = SSIB(T)%VEGP(35)
         TOPT(1) = SSIB(T)%VEGP(36)
         TOPT(2) = SSIB(T)%VEGP(37)
         TLL(1) = SSIB(T)%VEGP(38)
         TLL(2) = SSIB(T)%VEGP(39)
         TU(1) = SSIB(T)%VEGP(40)
         TU(2) = SSIB(T)%VEGP(41)
         DEFAC(1) = SSIB(T)%VEGP(42)
         DEFAC(2) = SSIB(T)%VEGP(43)
         PH1(1) = SSIB(T)%VEGP(44)
         PH1(2) = SSIB(T)%VEGP(45)
         PH2(1) = SSIB(T)%VEGP(46)
         PH2(2) = SSIB(T)%VEGP(47)
         ROOTD(1) = SSIB(T)%VEGP(48)
         ROOTD(2) = SSIB(T)%VEGP(49)
         BEE = SSIB(T)%VEGP(50)
         PHSAT = SSIB(T)%VEGP(51)
         SATCO = SSIB(T)%VEGP(52)
         POROS = SSIB(T)%VEGP(53)
         SLOPE = SSIB(T)%VEGP(54)
         ZDEPTH(1) = SSIB(T)%VEGP(55)
         ZDEPTH(2) = SSIB(T)%VEGP(56)
         ZDEPTH(3) = SSIB(T)%VEGP(57)

!========monthly parameters
         GREEN(1)  = SSIB(T)%VEGIP(1)
         GREEN(2)  = SSIB(T)%VEGIP(2)
         VCOVER(1) = SSIB(T)%VEGIP(3)
         VCOVER(2) = SSIB(T)%VEGIP(4)
         ZLT(1)    = SSIB(T)%VEGIP(5)
         ZLT(2)    = SSIB(T)%VEGIP(6)
         Z0        = SSIB(T)%VEGIP(7)
         D         = SSIB(T)%VEGIP(8)
         Z1        = SSIB(T)%VEGIP(9)
         Z2        = SSIB(T)%VEGIP(10)
         RBC       = SSIB(T)%VEGIP(11)
         RDC       = SSIB(T)%VEGIP(12)

!======other parameters (state)
         TC       = SSIB(T)%TCINI
         TGS      = SSIB(T)%TGSINI
         TD       = SSIB(T)%TDINI
         TA       = SSIB(T)%TAINI
         TM       = SSIB(T)%TMINI
         HT       = SSIB(T)%HTINI
         QA       = SSIB(T)%QAINI
         WWW(1)   = SSIB(T)%WWWINI(1)
         WWW(2)   = SSIB(T)%WWWINI(2)
         WWW(3)   = SSIB(T)%WWWINI(3)
         CAPAC(1) = SSIB(T)%CAPACINI(1)
         CAPAC(2) = SSIB(T)%CAPACINI(2)
         ZWIND    = SSIB(T)%ZWINDINI
         ZMET     = SSIB(T)%ZMETINI
         ILW      = SSIB(T)%ILWINI
         ITRUNK   = SSIB(T)%ITRUNKINI

         VWIND  = (ssib(t)%FORCING(6))*(ssib(t)%FORCING(6))
         UWIND  = (ssib(t)%FORCING(5))*(ssib(t)%FORCING(5))
         SFCSPD = SQRT( UWIND + VWIND )
         swdown = ssib(t)%FORCING(3)
         rnetm =  ssib(t)%FORCING(4)
         tprec =  ssib(t)%FORCING(8)
         cprec =  ssib(t)%FORCING(9)
         tm =  ssib(t)%FORCING(1)
         um =  SFCSPD
         pr =  ssib(t)%FORCING(7)/100.
         QM =  ssib(t)%FORCING(2)
!         em = em * 98.59 *10. /0.622
         em = QM * pr / (0.622 + QM)

              ! specific humidity [kg/kg]  to vapor pressure [mb]
              !                   measured k83 pressure 98.59 Kpa

         um = amax1(um,0.25)
!         ustarm = mustar/100.
!         swdown = amax1(swdown,0.1)
         ppl = tprec - cprec
         ppc = cprec

         startsm = ((www(1) * zdepth(1)) + (www(2) * zdepth(2)) + &
                    (www(3) * zdepth(3))) * poros
         startint = capac(1)
         startswe = capac(2)

         if (prin) then
            print *,' '
            print *,'Latitude ',grid(tile(t)%index)%lat, &
                  '; Longitude ',grid(tile(t)%index)%lon
            print *,' '
            print *,'Constants: '
            print *,'tf: ',tf
            print *,' '
            print *,'Forcing data: '
            print *,'tm: ',tm
            print *,'pr: ',pr
            print *,'qm: ',qm
            print *,'um: ',um
            print *,'ppl: ',ppl
            print *,'ppc: ',ppc
            print *,'swdown: ',swdown
            print *,'rnetm: ',rnetm
            print *,' '
         endif

!	ttii = WWW(3) + 100.
!      IF (ttii.eq.WWW(3)) THEN 
!          WRITE (*,*) '--> Variable ',WWW(3),t,dtt,ttii
!	  STOP
!       ENDIF !
!======end of assigning parameters
!      write(6,*) 'SOIL WETNESS FRACTION INITIALISATION: ',www(1),www(2),www(3)

         iyearh=lis%t%yr
         imonthh=lis%t%mo
         idayh=lis%t%da
         isech=(lis%t%hr*3600)+(lis%t%mn*60)
         if (prin) print *,'asdf',iyearh,imonthh,idayh,isech

         CALL ZENITHH(grid(tile(t)%index)%lat,grid(tile(t)%index)%lon, &
                      dtt,iyearh,imonthh,idayh,isech,angulosolar,prin)

!=== THE FOLLOWING BREAKS DOWN THE FORCING VARIABLES
!     THE FOLLOWING PROGRAM ONLY NEED TO BE CALLED ONCE FOR A LATITUDE
         IF (int(TILE(t)%VEGT).EQ.12) then
            CALL CROPS(grid(tile(t)%index)%lat,dum,lis%t%doy,NSX,VCOVER)
!,CHIL,ZLT,GREEN,RSTPAR,TOPT,TL,TU,DEFAC,PH2,PH1)
         endif

         RHOAIR = (pr*100.0) / GASR / TM

         SUNANG = angulosolar

         cloud = (1160.*sunang - swdown) / (963. * sunang)
         cloud = amax1(cloud,0.)
         cloud = amin1(cloud,1.)

         difrat = 0.0604 / ( sunang-0.0223 ) + 0.0683
         if ( difrat .lt. 0. ) difrat = 0.
         if ( difrat .gt. 1. ) difrat = 1.

         difrat = difrat + ( 1. - difrat ) * cloud
         vnrat = ( 580. - cloud*464. ) / ( ( 580. - cloud*499. ) &
               + ( 580. - cloud*464. ) )
         FRAC(1,1) = (1.-DIFRAT)*VNRAT
         FRAC(1,2) = DIFRAT*VNRAT
         FRAC(2,1) = (1.-DIFRAT)*(1.-VNRAT)
         FRAC(2,2) = DIFRAT*(1.-VNRAT)

         radn(1,1) = (1.-difrat)*vnrat*swdown
         radn(1,2) = difrat*vnrat*swdown
         radn(2,1) = (1.-difrat)*(1.-vnrat)*swdown
         radn(2,2) = difrat*(1.-vnrat)*swdown
         radn(3,2) = rnetm

         if (prin) then
            print *,'Radsplit: '
            print *,'radn_11: ',radn(1,1)
            print *,'radn_12: ',radn(1,2)
            print *,'radn_21: ',radn(2,1)
            print *,'radn_22: ',radn(2,2)
            print *,'sunang: ',sunang
         endif

         ppl = ppl * DTT
         ppc = ppc * DTT

!=== END OF BREAKING DOWN THE FORCING VARIABLES (SUBROUTINE DRIVER)
!=== ===============================================================
!=== ===============================================================
!=== ============BEGIN OF SIMULATION===============================

         CALL RADAB

         if (prin) then
            print *,'salb: ',salb(1,1),salb(1,2),salb(2,1),salb(2,2)
            print *,' '
            print *,'Initial conditions: '
            print *,'capac_1: ',capac(1)
            print *,'capac_2: ',capac(2)
            print *,'w_1: ',www(1)
            print *,'w_2: ',www(2)
            print *,'w_3: ',www(3)
            print *,'waterinlayer_1: ',(www(1) * poros * zdepth(1))
            print *,'waterinlayer_2: ',(www(2) * poros * zdepth(2))
            print *,'waterinlayer_3: ',(www(3) * poros * zdepth(3))
            print *,' '
            print *,'--------------------------------------------------'
            print *,' '
            print *,'w_1: ',www(1)
            print *,'w_2: ',www(2)
            print *,'w_3: ',www(3)
            print *,'phsat: ',phsat
            print *,'bee: ',bee
            print *,'CALLING ROOT'
         endif

         CALL ROOT1

         if (prin) then
            print *,'phsoil_1: ',phsoil(1)
            print *,'phsoil_2: ',phsoil(2)
            print *,'phsoil_3: ',phsoil(3)

            print *,' '
            print *,'--------------------------------------------------'
            print *,' '
            print *,'stefan: ',stefan
            print *,'vcover_1: ',vcover(1)
            print *,'vcover_2: ',vcover(2)
            print *,'thermk: ',thermk
            print *,'tc: ',tc
            print *,'tgs: ',tgs
            print *,'radn_11: ',radn(1,1)
            print *,'radn_12: ',radn(1,2)
            print *,'radn_21: ',radn(2,1)
            print *,'radn_22: ',radn(2,2)
            print *,'CALLED RADUSE'
            print *,'radt_1: ',radt(1)
            print *,'radt_2: ',radt(2)
            print *,'par: ',par(1)
            print *,'pd: ',pd(1)

            print *,' '
            print *,'--------------------------------------------------'
            print *,' '
            print *,'ityp: ',TILE(t)%VEGT
            print *,'zlt_1: ',zlt(1)
            print *,'zlt_2: ',zlt(2)
            print *,'green: ',green(1)
            print *,'vcover_1: ',vcover(1)
            print *,'vcover_2: ',vcover(2)
            print *,'chil: ',chil(1)
            print *,'rstpar_1: ',rstpar(1,1)
            print *,'rstpar_2: ',rstpar(1,2)
            print *,'rstpar_3: ',rstpar(1,3)
            print *,'sunang: ',sunang
            print *,'par: ',par(1)
            print *,'pd: ',pd(1)
            print *,'CALLING STOMAT'
         endif

         CALL STOMAT1                                        
         RSTUN = RST(1)

! ** WATER BALANCE CHECK
         TOTWB = WWW(1) * POROS * ZDEPTH(1) &				 
     	       + WWW(2) * POROS * ZDEPTH(2)  &	 			
     	       + WWW(3) * POROS * ZDEPTH(3)  &	 			
     	       + CAPAC(1) + CAPAC(2)

         if (prin) then
            print *,'rst_1: ',rst(1)
            print *,'rst_2: ',rst(2)

            print *,' '
            print *,'--------------------------------------------------'
            print *,' '
            print *,'dt: ',dtt
            print *,'snomel: ',snomel
            print *,'clai: ',clai
            print *,'cw: ',cw
            print *,'tf: ',tf
            print *,'poros: ',poros
            print *,'satco: ',satco
            print *,'zdepth_1: ',zdepth(1)
            print *,'zdepth_2: ',zdepth(2)
            print *,'zdepth_3: ',zdepth(3)
            print *,'capac_1: ',capac(1)
            print *,'capac_2: ',capac(2)
            print *,'w_1: ',www(1)
            print *,'w_2: ',www(2)
            print *,'w_3: ',www(3)
            print *,'satcap_1: ',satcap(1)
            print *,'satcap_2: ',satcap(2)
            print *,'tm: ',tm
            print *,'ppl: ',ppl
            print *,'ppc: ',ppc
            print *,'tc: ',tc
            print *,'tgs: ',tgs
            print *,'cg: ',cg
            print *,'CALLING INTERC'
         endif

         CALL INTERC                                                           

         if (prin) then
            print *,'tc: ',tc
            print *,'tgs: ',tgs
            print *,'capac_1: ',capac(1)
            print *,'capac_2: ',capac(2)
            print *,'cc: ',ccx
            print *,'cg: ',cg
            print *,'roff: ',roff
         endif

         currwb = (((www(1) * zdepth(1)) + (www(2) * zdepth(2)) + &
                    (www(3) * zdepth(3))) * poros) + capac(1)   + &
                    capac(2) + roff - ((ppl + ppc) / 1000.0)
         if ((abs(currwb-totwb).gt.0.000001).or.(prin)) then
            print *,' '
            print *,'INTERC CAPAC DIFF: ',(currwb-totwb)
            print *,'point: ',t,'  ',lonco,'  ',latco
            print *,'currtotal: ',currwb
            print *,'starttotal: ',totwb
            print *,'capac_1: ',capac(1)
            print *,'capac_2: ',capac(2)
            print *,'w_1: ',www(1)
            print *,'w_2: ',www(2)
            print *,'w_3: ',www(3)
            print *,'waterinlayer_1: ',(www(1) * poros * zdepth(1))
            print *,'waterinlayer_2: ',(www(2) * poros * zdepth(2))
            print *,'waterinlayer_3: ',(www(3) * poros * zdepth(3))
            print *,'ppl: ',(ppl/1000.)
            print *,'ppc: ',(ppc/1000.)
            print *,'roff: ',roff
            print *,'tm: ',tm
            print *,'tc: ',tc
            print *,'tgs: ',tgs
         endif
         if (abs(currwb-totwb).gt.0.000001) then
            print *,'interc capac error!'
            stop
         endif

         if (prin) then
            print *,'tc: ',tc
            print *,'tgs: ',tgs
            print *,'tm: ',tm
            print *,'qm: ',qm
            print *,'CALLING TEMRES'
         endif

         currwb = (((www(1) * zdepth(1)) + (www(2) * zdepth(2)) + &
                    (www(3) * zdepth(3))) * poros) + capac(1)   + &
                    capac(2) + roff - ((ppl + ppc) / 1000.0)
         if ((abs(currwb-totwb).gt.0.000001).or.(prin)) then
            print *,' '
            print *,'TEMRES CAPAC DIFF: ',(currwb-totwb)
            print *,'point: ',t,'  ',lonco,'  ',latco
            print *,'currtotal: ',currwb
            print *,'starttotal: ',totwb
            print *,'capac_1: ',capac(1)
            print *,'capac_2: ',capac(2)
            print *,'w_1: ',www(1)
            print *,'w_2: ',www(2)
            print *,'w_3: ',www(3)
            print *,'waterinlayer_1: ',(www(1) * poros * zdepth(1))
            print *,'waterinlayer_2: ',(www(2) * poros * zdepth(2))
            print *,'waterinlayer_3: ',(www(3) * poros * zdepth(3))
            print *,'ppl: ',(ppl/1000.)
            print *,'ppc: ',(ppc/1000.)
            print *,'roff: ',roff
            print *,'tm: ',tm
            print *,'tc: ',tc
            print *,'tgs: ',tgs
         endif
         if (abs(currwb-totwb).gt.0.000001) then
            print *,'temres capac error!'
            stop
         endif

      CALL TEMRS1(prin)

         currwb = (((www(1) * zdepth(1)) + (www(2) * zdepth(2)) + &
                    (www(3) * zdepth(3))) * poros) + capac(1)   + &
                    capac(2) + roff - ((ppl + ppc) / 1000.0)
         if ((abs(currwb-totwb).gt.0.000001).or.(prin)) then
            print *,' '
            print *,'TEMRES CAPAC DIFF: ',(currwb-totwb)
            print *,'point: ',t,'  ',lonco,'  ',latco
            print *,'currtotal: ',currwb
            print *,'starttotal: ',totwb
            print *,'capac_1: ',capac(1)
            print *,'capac_2: ',capac(2)
            print *,'w_1: ',www(1)
            print *,'w_2: ',www(2)
            print *,'w_3: ',www(3)
            print *,'waterinlayer_1: ',(www(1) * poros * zdepth(1))
            print *,'waterinlayer_2: ',(www(2) * poros * zdepth(2))
            print *,'waterinlayer_3: ',(www(3) * poros * zdepth(3))
            print *,'ppl: ',(ppl/1000.)
            print *,'ppc: ',(ppc/1000.)
            print *,'roff: ',roff
            print *,'tm: ',tm
            print *,'tc: ',tc
            print *,'tgs: ',tgs
         endif
         if (abs(currwb-totwb).gt.0.000001) then
            print *,'temres capac error!'
            stop
         endif

         if (prin) then
            print *,' '
            print *,'--------------------------------------------------'
            print *,' '
            print *,'tc: ',tc
            print *,'tgs: ',tgs
            print *,'tm: ',tm
            print *,'qm: ',qm
            print *,'CALLING UPDATE'
         endif

         CALL UPDAT1(prin)

         ENDWB = WWW(1) * POROS * ZDEPTH(1)  &				 
               + WWW(2) * POROS * ZDEPTH(2)   &	 			
               + WWW(3) * POROS * ZDEPTH(3)   &	 			
               + CAPAC(1) + CAPAC(2) - (PPL+PPC)/1000. + ETMASS/1000. + ROFF 
         ERROR = TOTWB - ENDWB   

         IF (ABS(ERROR).GT.0.0001) THEN
            WRITE(24,882) ERROR
            print *,' '
            print *,'Water Balance ERROR!'
            print *,'Error: ',ERROR
            print *,'latco: ',latco
            print *,'lonco: ',lonco
            print *,'www_1: ',www(1)
            print *,'www_2: ',www(2)
            print *,'www_3: ',www(3)
            print *,'waterinlayer_1: ',(www(1) * poros * zdepth(1))
            print *,'waterinlayer_2: ',(www(2) * poros * zdepth(2))
            print *,'waterinlayer_3: ',(www(3) * poros * zdepth(3))
            print *,'capac_1: ',capac(1)
            print *,'capac_2: ',capac(2)
            print *,'ppl: ',(ppl/1000.)
            print *,'ppc: ',(ppc/1000.)
            print *,'etmass: ',(etmass/1000.)
            print *,'roff: ',roff
            stop
         endif
 882     FORMAT(1X,'WARNING WATER BALANCE AFTER UPDATE ',f8.6)

         CBAL = RADT(1) - CHF - (ECT+HC+ECI)/DTT                           
         GBAL = RADT(2) - SHF - (EGT+EGI+HG+EGS)/DTT                       
         ZLHS = RADT(1) + RADT(2) - CHF - SHF
         ZRHS = HFLUX + (ECT + ECI + EGT + EGI + EGS)/DTT  

         IF (ABS(ZLHS-ZRHS).GT.1.) WRITE(24,881) ABS(ZLHS-ZRHS)
 881     FORMAT(1X,' WARNING ENERGY BALANCE ',f8.6)

! latent heat flux
         ectw = ect / dtt
         eciw = eci / dtt
         egtw = egt / dtt
         egiw = egi / dtt
         egsw = egs / dtt
         if (prin) then
            print *,'ect: ',ect
            print *,'eci: ',eci
            print *,'egt: ',egt
            print *,'egi: ',egi
            print *,'egs: ',egs
            print *,'dtt: ',dtt
         endif

! soil moisture
         w1i = WWW(1) * POROS * ZDEPTH(1) * 1000.
         w2i = WWW(2) * POROS * ZDEPTH(2) * 1000.
         w3i = WWW(3) * POROS * ZDEPTH(3) * 1000.
         if (swdown.gt.0.0) then
            suralbedo = ((radn(1,1) * salb(1,1)) + (radn(1,2) * salb(1,2)) + &
                         (radn(2,1) * salb(2,1)) + (radn(2,2) * salb(2,2))) / swdown
         else
            suralbedo = LIS%d%UDEF
         endif
         soilwet = (zdepth(1) + zdepth(2) + zdepth(3)) * &
                   ((-exp(ph2(1)) / phsat) ** (-bee))
         soilwet = (((w1i + w2i + w3i) / poros / 1000.) - soilwet) / &
                   ((zdepth(1) + zdepth(2) + zdepth(3)) - soilwet)

!=== ===============END OF SIMULATION===============================
!===================================================================
!===================================================================

!      PRINT*,'  --------------------------------------'
!      PRINT*,'  State Variables '
!      PRINT*,'  --------------------------------------'
!      WRITE(*,*) SSIB(T)%T1,' T1...Skin temperature (K)'
!      WRITE(*,*)(SSIB(T)%STC(IJ), IJ=1,NSOIL),' STC'
!      WRITE(*,*)(SSIB(T)%SMC(IJ), IJ=1,NSOIL),' SMC'
!      WRITE(*,*)(SSIB(T)%SH2O(IJ), IJ=1,NSOIL),' SH2O'
!      WRITE(*,*) SSIB(T)%CMC,' CMC...Canopy water content (m)'
!      WRITE(*,*) SSIB(T)%SNOWH,' SNOWH...Actual snow depth (m)'
!      WRITE(*,*) SSIB(T)%SNEQV,' SNEQV...Water equiv snow depth (m)'
!      WRITE(*,*) 'CH= ',SSIB(T)%CH,'   CM= ',SSIB(T)%CM
!      PRINT*,'  --------------------------------------'

!=== Collect the output variables into SSIB(T)%RETURN
!      IF (ABS(ERROR).GT.0.0001) THEN
!         write(24,*) ERROR, (WWW(1) * POROS * ZDEPTH(1)),  &
!                            (WWW(2) * POROS * ZDEPTH(2)),  &
!                            (WWW(3) * POROS * ZDEPTH(3)),  &
!         CAPAC(1) , CAPAC(2) , PPL/1000. , ETMASS/1000. , ROFF 
!      ENDIF

!=== Collect state variables

      SSIB(T)%TCINI	   = TC       ! CANOPY TEMPERATURE (K)
      SSIB(T)%TGSINI	   = TGS      ! SOIL SURFACE TEMPERATURE (K)
      SSIB(T)%TDINI	   = TD       ! DEEP SOIL TEMPERATURE (K)
      SSIB(T)%TAINI	   = TA       ! TEMPERATURE AT CANOPY AIR SPACE (K)
      SSIB(T)%TMINI	   = TM       ! TEMPERATURE AT LOWEST MODEL LAYER (K)
      SSIB(T)%HTINI	   = HT
      SSIB(T)%QAINI	   = QA
      SSIB(T)%WWWINI(1)    = WWW(1)   !SOIL MOISTURE
      SSIB(T)%WWWINI(2)    = WWW(2)
      SSIB(T)%WWWINI(3)    = WWW(3)
      SSIB(T)%CAPACINI(1)  = CAPAC(1) !INTERCEPTION AT CANOPY
      SSIB(T)%CAPACINI(2)  = CAPAC(2) !SNOW DEPTH

!=== Collect the ALMA output variables
      SSIB(T)%swnet        = SSIB(T)%swnet    + swcan + swgnd
      SSIB(T)%lwnet        = SSIB(T)%lwnet    + radn(3,2) - zlwup
      SSIB(T)%qle          = SSIB(T)%qle      + (ectw+eciw+egtw+egiw+egsw)
      SSIB(T)%qh           = SSIB(T)%qh       + hflux
      SSIB(T)%qg           = SSIB(T)%qg       + shf
      SSIB(T)%qf           = SSIB(T)%qf       + (smelt * snomel / dtt)
      SSIB(T)%qtau         = SSIB(T)%qtau     + (drag * um)
      SSIB(T)%delsurfheat  = SSIB(T)%delsurfheat + (chf*dtt)
      if (ssib(t)%FORCING(1).lt.tf) then
         SSIB(T)%snowf        = SSIB(T)%snowf + tprec
      else
         SSIB(T)%rainf        = SSIB(T)%rainf + tprec
      endif
      SSIB(T)%evap         = SSIB(T)%evap     + (etmass/dtt)
      SSIB(T)%qs           = SSIB(T)%qs       + max((roff - (q3g * dtt)),0.0) * 1000. / dtt
      SSIB(T)%qsb          = SSIB(T)%qsb      + (q3g * dtt) * 1000. / dtt
      SSIB(T)%qsm          = SSIB(T)%qsm      + smelt * 1000. / dtt
      SSIB(T)%delsoilmoist = SSIB(T)%delsoilmoist + ((((www(1) * zdepth(1)) + &
                                                       (www(2) * zdepth(2)) + &
                                                       (www(3) * zdepth(3))) * &
                                                        poros) - startsm) * 1000.
      SSIB(T)%delswe       = SSIB(T)%delswe       + (capac(2) - startswe) * 1000.
      SSIB(T)%delintercept = SSIB(T)%delintercept + (capac(1) - startint) * 1000.
! These variables are either instantanous or time-averaged,
! depending on value of STATEVAR_AVG flag in lis.crd namelist.
      if (ssibdrv%STATEVAR_AVG.eq.0) then
         SSIB(T)%vegtc           = tc
         SSIB(T)%baresoilt       = tgs
         SSIB(T)%avgsurft        = (tgs * (1.0 - vcover(1))) + (tc * vcover(1))
         SSIB(T)%radteff         = tgeff
         if (swdown.gt.0.0) then
            SSIB(T)%albedo       = suralbedo
         else
            SSIB(T)%albedo       = LIS%d%UDEF
         endif
         SSIB(T)%swe             = capac(2) * 1000.
         if (tc.lt.tf) then
            SSIB(T)%sweveg       = capac(1) * 1000.
         else
            SSIB(T)%sweveg       = 0.0
         endif
         SSIB(T)%soilmoist1      = w1i
         SSIB(T)%soilmoist2      = w2i
         SSIB(T)%soilmoist3      = w3i
         SSIB(T)%soiltemp        = td
         SSIB(T)%soilwet         = soilwet
      else
         SSIB(T)%vegtc           = SSIB(T)%vegtc       + tc
         SSIB(T)%baresoilt       = SSIB(T)%baresoilt   + tgs
         SSIB(T)%avgsurft        = SSIB(T)%avgsurft    + (tgs * (1.0 - vcover(1))) + (tc * vcover(1))
         SSIB(T)%radteff         = SSIB(T)%radteff     + tgeff
         if (swdown.gt.0.0) then
            SSIB(T)%albedo       = SSIB(T)%albedo      + suralbedo
            SSIB(T)%albedocount  = SSIB(T)%albedocount + 1
         endif
         SSIB(T)%swe             = SSIB(T)%swe         + (capac(2) * 1000.)
         if (tc.lt.tf) then
            SSIB(T)%sweveg       = SSIB(T)%sweveg      + (capac(1) * 1000.)
         else
            SSIB(T)%sweveg       = SSIB(T)%sweveg      + 0.0
         endif
         SSIB(T)%soilmoist1      = SSIB(T)%soilmoist1  + w1i
         SSIB(T)%soilmoist2      = SSIB(T)%soilmoist2  + w2i
         SSIB(T)%soilmoist3      = SSIB(T)%soilmoist3  + w3i
         SSIB(T)%soiltemp        = SSIB(T)%soiltemp    + td
         SSIB(T)%soilwet         = SSIB(T)%soilwet     + soilwet
      endif
! Back to typical output variables
      SSIB(T)%ecanop       = SSIB(T)%ecanop    + ((eciw + egiw) / hlat)
      SSIB(T)%tveg         = SSIB(T)%tveg      + ((ectw + egtw) / hlat)
      SSIB(T)%esoil        = SSIB(T)%esoil     + (egsw / hlat)
      SSIB(T)%rootmoist    = SSIB(T)%rootmoist + (w1i + w2i)
      SSIB(T)%canopint     = SSIB(T)%canopint  + (capac(1) * 1000.)
      SSIB(T)%acond        = SSIB(T)%acond     + (1.0 / ra)
      SSIB(T)%snowfrac     = SSIB(T)%snowfrac  + scov2

      enddo                     ! end di_array

      ssib%count = ssib%count + 1

!>>>  END OF SSIB_MAIN <<<<<<<<
!      stop

      return
!EOC
      end subroutine ssib_main

!*** SSIB SUBROUTINES ****************************************************

!=======================================================================
!                                                                       
      SUBROUTINE CONSTS                                                 
!                                                          1 AUGUST 1988
!=======================================================================
!                                                                       
!     INITIALIZATION OF PHYSICAL CONSTANTS                              
!                                                                       
!-----------------------------------------------------------------------
      use ssib_common                                                  
!                                                                       
      PSUR     = 1000.                                                 
      CPAIR    = 1010.                                                  
      RHOAIR   = 1.225                                                  
      STEFAN   = 5.669 * 10E-9                                          
      GRAV     = 9.81                                                   
      VKC      = 0.378                                                  
      PIE      = 3.14159265                                             
      TIMCON   = PIE/86400.                                             
      CLAI     = 4.2 * 1000. * 0.2                                      
      CW       = 4.2 * 1000. * 1000.                                    
      TF       = 273.15
      GASR     = 287.05                                               
!-----------------------------------------------------------------------
!     N.B. : HLAT IS EXPRESSED IN J KG-1                                
!            SNOMEL IS EXPRESSED IN J M-1                               
!-----------------------------------------------------------------------
      HLAT     = ( 3150.19 - 2.378 * TM ) * 1000.                       
      SNOMEL   = 370518.5 * 1000.                                       
!      PSY      = CPAIR / HLAT * PSUR / .622                             

      RETURN                                                            
      END                                                               

!=======================================================================
!                                                                       
      SUBROUTINE INTERC                                                 
!                                                          1 AUGUST 1988
!=======================================================================
!                                                                       
!     CALCULATION OF (1) INTERCEPTION AND DRAINAGE OF RAINFALL AND SNOW 
!                    (2) SPECIFIC HEAT TERMS FIXED FOR TIME STEP        
!                                                                       
!     MODIFICATION 30 DEC 1985 : NON-UNIFORM PRECIPITATION             
!     ------------      CONVECTIVE PPN. IS DESCRIBED BY AREA-INTENSITY  
!                       RELATIONSHIP :-                                 
!                                                                       
!                                        F(X) = A*EXP(-B*X)+C           
!                                                                       
!                       THROUGHFALL, INTERCEPTION AND INFILTRATION      
!                       EXCESS ARE FUNCTIONAL ON THIS RELATIONSHIP      
!                       AND PROPORTION OF LARGE-SCALE PPN.              
!---------------------------------------------------------------------- 
      use ssib_common                                                  
!                                                                       
      REAL CAPACP(2), SNOWP(2), PCOEFS(2,2)                        
      DATA PCOEFS(1,1)/ 20. /, PCOEFS(1,2)/ .206E-8 /,&                  
          PCOEFS(2,1)/ 0.0001 /, PCOEFS(2,2)/ 0.9999 /, BP /20. /      
!                                                                       
      AP = PCOEFS(2,1)                                                  
      CP = PCOEFS(2,2)                                                  
      TOTALP = PPC + PPL                                                
      IF(TOTALP.LT.1.E-8)GO TO 6000                                     
      AP = PPC/TOTALP * PCOEFS(1,1) + PPL/TOTALP * PCOEFS(2,1)          
      CP = PPC/TOTALP * PCOEFS(1,2) + PPL/TOTALP * PCOEFS(2,2)          
 6000 CONTINUE                                                          
!                                                                       
      ROFF = 0.                                                         
      THRU = 0.                                                         
      FPI  = 0.                                                         
!                                                                       
!---------------------------------------------------------------------- 
!     THERMAL CONDUCTIVITY OF THE SOIL, TAKING INTO ACCOUNT POROSITY    
!---------------------------------------------------------------------- 
!                                                                       
      THETA=WWW(1)*POROS                                                
      CHISL=( 9.8E-4+1.2E-3*THETA )/( 1.1-0.4*THETA )                   
      CHISL=CHISL*4.186E2                                               
!                                                                       
!---------------------------------------------------------------------- 
!     THERMAL DIFFUSIVITY AND HEAT CAPACITYOF THE SOIL                  
!---------------------------------------------------------------------- 
!                                                                       
      DIFSL=5.E-7                                                       
!                                                                       
      ROCS =CHISL/DIFSL                                                 
      D1   =SQRT(DIFSL*86400.0)                                         
      CSOIL=ROCS*D1/SQRT(PIE)/2.0 
      THALAS=0.
      OCEANS=0.
      POLAR=0.                                 
      CSOIL=CSOIL*(1.0-THALAS)+10.E10*OCEANS+POLAR*3.6*4.2E4            
!                                                                       
!                                                                       
      P0 = TOTALP * 0.001                                               
!                                                                       
!---------------------------------------------------------------------- 
!     INPUT PRECIPITATION IS GIVEN IN MM, CONVERTED TO M TO GIVE P0.    
!---------------------------------------------------------------------- 
!                                                                       
      DO 1000 IVEG = 1, 2                                               
!                                                                       
      SPWET1 = AMIN1 ( 0.05, CAPAC(IVEG))*CW                            
!                                                                       
      TS = TC                                                           
      SPECHT = ZLT(1) * CLAI                                            
      IF ( IVEG .EQ. 1 ) GO TO 1100                                     
      TS = TGS                                                          
      SPECHT = CSOIL                                                    
1100  CONTINUE                                                          
!                                                                       
      XSC = AMAX1(0., CAPAC(IVEG) - SATCAP(IVEG) )                      
      IF(IVEG.EQ.2 .AND. TS.LE.TF )GO TO 1170                           
      CAPAC(IVEG) = CAPAC(IVEG) - XSC                                   
      ROFF = ROFF + XSC                                                 
      RNOFFS = XSC*1000. + RNOFFS
1170  CONTINUE                                                          
      CAPACP(IVEG) = 0.                                                 
      SNOWP(IVEG) = 0.                                                  
!                                                                       
      IF( TS .GT. TF ) CAPACP(IVEG) = CAPAC(IVEG)                       
      IF( TS .LE. TF ) SNOWP(IVEG) = CAPAC(IVEG)                        
      CAPAC(IVEG) = CAPACP(IVEG)                                        
      SNOWW(IVEG) = SNOWP(IVEG)                                         
      ZLOAD = CAPAC(IVEG) + SNOWW(IVEG)                                 
!                                                                       
      FPI = ( 1.-EXP( - EXTK(IVEG,3,1) * ZLT(IVEG)/VCOVER(IVEG) ) )&     
           * VCOVER(IVEG)                                              
      TTI = P0 * ( 1.-FPI )                                            
!                                                                       
!---------------------------------------------------------------------- 
!    PROPORTIONAL SATURATED AREA (XS) AND LEAF DRAINAGE(TEX)            
!---------------------------------------------------------------------- 
!                                                                       
      XS = 1.                                                           
      IF ( P0 .LT. 1.E-9 ) GO TO 1150                                   
      ARG =  ( SATCAP(IVEG)-ZLOAD )/( P0*FPI*AP ) -CP/AP                
      IF ( ARG .LT. 1.E-9 ) GO TO 1150                                  
      XS = -1./BP * ALOG( ARG )                                         
      XS = AMIN1( XS, 1. )                                              
      XS = AMAX1( XS, 0. )                                              
1150  TEX = P0*FPI * ( AP/BP*( 1.- EXP( -BP*XS )) + CP*XS ) - &          
           ( SATCAP(IVEG) - ZLOAD ) * XS                               
      TEX = AMAX1( TEX, 0. )                                            
!                                                                       
!---------------------------------------------------------------------- 
!    TOTAL THROUGHFALL (THRU) AND STORE AUGMENTATION                    
!---------------------------------------------------------------------- 
!                                                                       
      THRU = TTI + TEX                                                  
      IF(IVEG.EQ.2.AND.TGS.LE.TF)THRU = 0.                              
!                                                                       
      PINF = P0 - THRU                                                  
      IF( TM .GT. TF ) CAPAC(IVEG) = CAPAC(IVEG) + PINF                 
      IF( TM .LE. TF ) SNOWW(IVEG) = SNOWW(IVEG) + PINF                 
!                                                                       
      IF( IVEG .EQ. 1 ) GO TO 1300                                      
      IF( TM .GT. TF ) GO TO 1200                                      
      SNOWW(IVEG) = SNOWP(IVEG) + P0                                    
      THRU = 0.                                                         
      GO TO 1300                                                        
!                                                                       
!---------------------------------------------------------------------- 
!    INSTANTANEOUS OVERLAND FLOW CONTRIBUTION ( ROFF )                  
!---------------------------------------------------------------------- 
!                                                                       
1200  EQUDEP = SATCO * DTT                                              
!                                                                       
      XS = 1.                                                           
      IF ( THRU .LT. 1.E-9 ) GO TO 1250                                 
      ARG = EQUDEP / ( THRU * AP ) -CP/AP                               
      IF ( ARG .LT. 1.E-9 ) GO TO 1250                                  
      XS = -1./BP * ALOG( ARG )                                         
      XS = AMIN1( XS, 1. )                                              
      XS = AMAX1( XS, 0. )                                              
1250  ROFFO = THRU * ( AP/BP * ( 1.-EXP( -BP*XS )) + CP*XS )&            
            -EQUDEP*XS                                                 
      ROFFO = AMAX1 ( ROFFO, 0. )                                       
      ROFF = ROFF + ROFFO                                               
      RNOFFS = RNOFFS + roffo*1000.
      filtr =  filtr + (THRU - ROFFO)   
      WWW(1) = WWW(1) + (THRU - ROFFO) / ( POROS*ZDEPTH(1) )            
1300  CONTINUE                                                          
!                                                                       
!---------------------------------------------------------------------- 
!    TEMPERATURE CHANGE DUE TO ADDITION OF PRECIPITATION                
!---------------------------------------------------------------------- 
!                                                                       
      DIFF = ( CAPAC(IVEG)+SNOWW(IVEG) - CAPACP(IVEG)-SNOWP(IVEG) )*CW 
      CCP = SPECHT + SPWET1                                             
      CCT = SPECHT + SPWET1 + DIFF                                      
!                                                                       
      TSD = ( TS * CCP + TM * DIFF ) / CCT                              
!                                                                       
      FREEZE = 0.                                                       
      IF ( TS .GT. TF .AND. TM .GT. TF ) GO TO 2000                     
      IF ( TS .LE. TF .AND. TM .LE. TF ) GO TO 2000                     
!                                                                       
      TTA = TS                                                          
      TTB = TM                                                          
      CCA = CCP                                                         
      CCB = DIFF                                                        
      IF ( TSD .GT. TF ) GO TO 2100                                     
!                                                                       
!---------------------------------------------------------------------- 
!    FREEZING OF WATER ON CANOPY OR GROUND                              
!---------------------------------------------------------------------- 
!                                                                       
      CCC = CAPACP(IVEG) * SNOMEL                                       
      IF ( TS .LT. TM ) CCC = DIFF * SNOMEL / CW                        
      TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT                       
!                                                                       
      FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )                 
      FREEZE = (AMIN1 ( CCC, FREEZE )) / SNOMEL                         
      IF(TSD .GT. TF)TSD = TF - 0.1                                     
!                                                                       
      GO TO 2000                                                        
!                                                                       
2100  CONTINUE                                                          
!                                                                       
!---------------------------------------------------------------------- 
!    MELTING OF SNOW ON CANOPY OR GROUND                                
!---------------------------------------------------------------------- 
!                                                                       
      CCC = - SNOWW(IVEG) * SNOMEL                                      
      IF ( TS .GT. TM ) CCC = - DIFF * SNOMEL / CW                      
!                                                                       
      TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT                       
!                                                                       
      FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )                 
      FREEZE = (AMAX1( CCC, FREEZE )) / SNOMEL                          
      IF(TSD .LE. TF)TSD = TF - 0.1                                     
!                                                                       
2000  CONTINUE
      SMELT = FREEZE
      SNOWW(IVEG) = SNOWW(IVEG) + FREEZE                                
      CAPAC(IVEG) = CAPAC(IVEG) - FREEZE                                
!                                                                       
      IF( IVEG .EQ. 1 ) TC = TSD                                        
      IF( IVEG .EQ. 2 ) TGS = TSD                                       
      IF( SNOWW(IVEG) .LT. 0.0000001 ) GO TO 3000                       
      ZMELT = 0.                                                        
!     modified to force water into soil. Xue Feb. 1994
      ZMELT = CAPAC(IVEG)                             
!     IF ( TD .GT. TF ) ZMELT = CAPAC(IVEG)                             
!     IF ( TD .LE. TF ) ROFF = ROFF + CAPAC(IVEG)                       
      CAPAC(IVEG) = 0.                                                  
      WWW(1) = WWW(1) + ZMELT / ( POROS * ZDEPTH(1) )                   
      filtr = filtr + ZMELT                    
!                                                                       
3000  CONTINUE                                                          
!                                                                       
      CAPAC(IVEG) = CAPAC(IVEG) + SNOWW(IVEG)                           
      SNOWW(IVEG) = 0.                                                  
!
!     **** LOAD PILPS PARAMETER
!
!      if (freeze.lt.0) snm=snm-freeze
      freeze=0.0
!                                                                       
      P0 = THRU                                                         
!                                                                       
1000  CONTINUE                                                          
!                                                                       
!---------------------------------------------------------------------- 
!    CALCULATION OF CANOPY AND GROUND HEAT CAPACITIES.                  
!    N.B. THIS SPECIFICATION DOES NOT NECESSARILY CONSERVE ENERGY WHEN  
!    DEALING WITH VERY LATGE SNOWPACKS.                                 
!---------------------------------------------------------------------- 
!                                                                       
      CCX = ZLT(1) * CLAI + CAPAC(1) * CW                               
      SPWET = AMIN1 ( 0.05, CAPAC(2))*CW                                
      CG = (CSOIL + SPWET)                                              
!                                                                       
      RETURN                                                            
      END                                                               
!====================================================================   
!                                                                       
      SUBROUTINE NEWTON(A1,Y,FINC,NOX,NONPOS,IWOLK,L)
!                                                                       
!====================================================================== 
! ** VERSION ACQUIRED FROM EROS 2/19/86.                                
!                                                                       
!=======================================================================
!                                                                       
! ** THE NEWTON RAPHSON ITERATIVE ROUTINE WILL BE USED TO GENERATE NEW  
! ** VALUES OF A1 IF DABSOLUTE VALUE OF Y IS GREATER THAN ERTOL;        
! ** A1 IS ESTIMATE, Y IS RESULTANT ERROR                               
! ** NEX IS EXIT CONDITION  (0=NO EXIT) OR (1 WHEN DABS(Y) LT ERTOL)    
! ** ERTOL IS THE DABSOLUTE VALUE OF Y NECESSARY TO OBTAIN AN EXIT      
! ** FINC IS INITIAL INCREMENT SIZE FOR SECOND ESTIMATE OF A1           
! ** NONPOS=0 IF QUANTITY TO BE MINIMIZED CAN BE LESS THAN ZERO;        
! ** NONPOS=1 IF QUANTITY CAN ONLY BE POSITIVE                          
! ** L IDENTIFIES WHICH QUANTITY IS BEING CALCULATED.                   
!                                                                       
! ** CONTROL VALUES: FINC,ERTOL,NOX,NONPOS,L:MUST BE SET BY USER        
!-----------------------------------------------------------------------
!                                                                       
       use ssib_common
       REAL   IWALK(3), NEX(3)                                     
       DATA CONS/1.0/                                                   
!                                                                       
       ERTOL = 0.05 * FINC                                              
       IWALK(L) = IWOLK                                                 
       NEX(L)=NOX                                                       
!                                              
       IF ( ITER(L) .GE. 490 ) GO TO 160                                
       IF (ERTOL .LT. 0.00000001) ERTOL=0.000001                        
       IF (ABS(Y) .LE. ERTOL) GO TO 150                                 
       IF((ABS(Y-Y1(L))).LE.0.01*ERTOL .AND. IWALK(L).EQ.0 ) GO TO 8    
!                                                                     
       IF(ABS(Y1(L)).GT.ERTOL) GO TO 1                                  
       A2(L)=A1                                                         
       A1=A1-Y                                                          
       NEX(L)=0                                                         
       Y1(L)=Y                                                          
       ITER(L)=1                                                        
       IF (IWALK(L) .EQ. 3) GO TO 101                                   
       IWALK(L)=0                                                       
       GO TO 101                                                        
   1   ITER(L)=ITER(L)+1                                                
       IF(ITER(L) .EQ. 10) IWALK(L)=1                                   
       IF(IWALK(L) .NE. 0) GO TO 2                                      
       IF(ABS(Y) .GT. ERTOL) GO TO 3                                    
       NEX(L)=1                                                         
       GO TO 150                                                        
   3   A=A1-Y*(A1-A2(L))/(Y-Y1(L))                                      
       IF(ABS(A-A1).GT.(10.0*FINC))&                                    
                 A=A1+10.0*FINC*SIGN(CONS,(A-A1))                      
       A2(L)=A1                                                         
       A1=A                                                             
       Y1(L)=Y                                                          
       GO TO 101                                                        
   2   IF(IWALK(L).EQ.2)GO TO 4                                         
       IF(IWALK(L).EQ.3) GO TO 6                                        
       IF(SIGN(CONS,Y).EQ.SIGN(CONS,Y1(L))) GO TO  3                    
       ZINC(L)=(A1-A2(L))/4.0                                           
       A1=A2(L)+ZINC(L)                                                 
       IWALK(L)=2                                                       
       NEX(L)=0                                                         
       GO TO 101                                                        
   4   IF(SIGN(CONS,Y) .EQ.SIGN(CONS,Y1(L))) GO TO 5                    
       ZINC(L)=-ZINC(L)/4.0                                             
       A2(L)=A1                                                         
       A1=A1+ZINC(L)                                                    
       NEX(L)=0                                                         
       Y1(L)=Y                                                          
       GO TO 101                                                        
   5   A2(L)=A1                                                         
       A1=A1+ZINC(L)                                                    
       Y1(L)=Y                                                          
       NEX(L)=0                                                         
       GO TO 101                                                        
   6   IF(SIGN(CONS,Y).EQ.SIGN(CONS,Y1(L))) GO TO 7                     
       IWALK(L)=1                                                       
       GO TO 2                                                          
   7   A2(L) = A1                                                       
       A1 = A1+FINC                                                     
       Y1(L)=Y                                                          
       NEX(L) = 0                                                       
       GO TO 101                                                        
   8   A1 = A1 + FINC*2.0                                               
       NEX(L)=0                                                         
       GO TO 101                                                        
160    CONTINUE                                                         
       WRITE(6,900) Y, A1                                               
!      STOP                                                             
 900   FORMAT ( 3X,' FAILURE TO CONVERGE AFTER 490 ITERATIONS',&         
      /, 3X,' Y = ',2G12.5,2X,I14)                                     
!                                                                       
 150   NEX(L) = 1                                                       
       ZINC(L)=0.0                                                      
       ITER(L) = 0                                                      
       IWALK(L)=0                                                       
       Y1(L)=0.0                                                        
       Y=0.0                                                            
       A2(L)=0.0                                                        
 101   CONTINUE                                                         
       IF(NONPOS.EQ.1.AND.A1.LT.0.0) A1=A2(L)/2.0                       
       NOX = NEX(L)                                                     
       IWOLK = IWALK(L)                                                 
!                                                                       
       RETURN                                                           
       END                                                              

!=======================================================================
!                                                                       
      SUBROUTINE RADAB
!                                                          1 AUGUST 1988
!=======================================================================
!                                                                       
!     CALCULATION OF ALBEDOS VIA TWO STREAM APPROXIMATION( DIRECT       
!     AND DIFFUSE ) AND PARTITION OF RADIANT ENERGY                     
!                                                                       
!-----------------------------------------------------------------------
!     DIMENSION XQQ(2,2,2),sr(2)
      use ssib_common                                                  
      REAL sr(2)
      data sr/0.85,0.65/
      REAL TRANC1(2), TRANC2(2), TRANC3(2)                         

      F = SUNANG                                                        

!---------------------------------------------------------------------- 
!     CALCULATION OF MAXIMUM WATER STORAGE VALUES.                      
!---------------------------------------------------------------------- 
      FMELT = 1.      
      nmm=1                                                  
      IF ( ABS(TF-TGS) .LT. 0.5 ) FMELT = 0.6                           
      SATCAP(1) =  ZLT(1) * 0.0001                                      
      SATCAP(2) =  ZLT(2) * 0.0001                                      
      DEPCOV = AMAX1( 0., (CAPAC(2)*5.-Z1) )                            
      DEPCOV = AMIN1( DEPCOV, (Z2-Z1)*0.95 )                            
      SATCAP(1) = SATCAP(1) * ( 1. - DEPCOV / ( Z2 - Z1 ) )             
!                                                                       
!---------------------------------------------------------------------- 
!        
      DO 1000 IWAVE = 1, 2                                              
!                                                                       
      DO 2000 IVDUM = 1, 2                                              
!                                                                       
      IF ( IVDUM .EQ. 1 ) IVEG = 2                                      
      IF ( IVDUM .EQ. 2 ) IVEG = 1                                      
!---------------------------------------------------------------------- 
!     MODIFICATION FOR EFFECT OF SNOW ON UPPER STOREY ALBEDO            
!         SNOW REFLECTANCE   = 0.80, 0.40 . MULTIPLY BY 0.6 IF MELTING  
!         SNOW TRANSMITTANCE = 0.20, 0.54                               
!         SNOW REFLECTANCE   = 0.85, 0.65 . MULTIPLY BY 0.6 IF MELTING  
!                                                                       
!---------------------------------------------------------------------- 
      SCOV = 0.                                                         
      IF( IVEG .EQ. 2 ) GO TO 100                                       
      IF( TC .LE. TF ) SCOV =  AMIN1( 0.5, CAPAC(1) / SATCAP(1) )       
100   CONTINUE                                                          
      REFF1 = ( 1. - SCOV ) * REF(IVEG,IWAVE,1) + SCOV * ( 1.2 - &       
             IWAVE * 0.4 ) * FMELT                                     
      REFF2 = ( 1. - SCOV ) * REF(IVEG,IWAVE,2) + SCOV * ( 1.2 - &       
             IWAVE * 0.4 ) * FMELT                                     
      TRAN1 = TRAN(IVEG,IWAVE,1) * ( 1. - SCOV )             &  	 
     	     + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT )  &	 	
     	     * TRAN(IVEG,IWAVE,1)				       
      TRAN2 = TRAN(IVEG,IWAVE,2) * ( 1. - SCOV )                  &	 
     	     + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT ) * 0.9 &	
     	     * TRAN(IVEG,IWAVE,2)				       
                                                                        
!---------------------------------------------------------------------- 
!                                                                       
      SCAT = GREEN(IVEG)*( TRAN1 + REFF1 ) +( 1. - GREEN(IVEG) ) * &     
            ( TRAN2 + REFF2)                                           
      CHIV = CHIL(IVEG)                                                 
!                            
      IF ( ABS(CHIV) .LE. 0.01 ) CHIV = 0.01                            
      AA = 0.5 - 0.633 * CHIV - 0.33 * CHIV * CHIV                      
      BB = 0.877 * ( 1. - 2. * AA )                                     
!                                                                       
      PROJ = AA + BB * F                                                
      EXTKB = ( AA + BB * F ) / F                                       
      ZMEW = 1. / BB * ( 1. - AA / BB * ALOG ( ( AA + BB ) / AA ) )     
      ACSS = SCAT / 2. * PROJ / ( PROJ + F * BB )                       
      ACSS = ACSS * ( 1. - F * AA / ( PROJ + F * BB ) * ALOG ( ( PROJ &  
            +   F * BB + F * AA ) / ( F * AA ) ) )                     
!                                                                       
      EXTK( IVEG, IWAVE, 1 ) = PROJ / F * SQRT( 1.-SCAT )               
      EXTK( IVEG, IWAVE, 2 ) = 1. / ZMEW * SQRT( 1.-SCAT )              
      EXTK( IVEG, 3, 1 ) = AA + BB                                      
      EXTK( IVEG, 3, 2 ) = 1./ZMEW                                      
!                                                                       
      UPSCAT = GREEN(IVEG) * TRAN1 + ( 1. - GREEN(IVEG) ) * TRAN2       
      UPSCAT = 0.5 * ( SCAT + ( SCAT - 2. * UPSCAT ) *  &                
              (( 1. - CHIV ) / 2. ) ** 2 )                             
!                                                                       
      BETAO = ( 1. + ZMEW * EXTKB ) / ( SCAT * ZMEW * EXTKB ) * ACSS    
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     DICKINSON'S VALUES                                                
!                                                                       
      BE = 1. - SCAT + UPSCAT                                           
      CE = UPSCAT                                                       
      BOT = ( ZMEW * EXTKB ) ** 2 + ( CE**2 - BE**2 )                   
      IF ( ABS(BOT) .GT. 1.E-10) GO TO 200                              
      SCAT = SCAT* 0.98                                                 
      BE = 1. - SCAT + UPSCAT                                           
      BOT = ( ZMEW * EXTKB ) ** 2 + ( CE**2 - BE**2 )                   
200   CONTINUE                                                          
      DE = SCAT * ZMEW * EXTKB * BETAO                                  
      FE = SCAT * ZMEW * EXTKB * ( 1. - BETAO )                         
!---------------------------------------------------------------------- 
!                                                                       
      CCE = DE * BE - ZMEW * DE * EXTKB + CE * FE                       
      FFE = BE * FE + ZMEW * FE * EXTKB + CE * DE                       
!                                                                       
      TORE = -CCE / BOT                                                 
      SIGE = -FFE / BOT                                                 
!                                                                       
      PSI = SQRT(BE**2 - CE**2)/ZMEW                                    
!                                                                       
!---------------------------------------------------------------------- 
!     REDUCTION IN EXPOSED HEIGHT OF UPPER STOREY AS SNOW ACCUMULATES   
!                                                                       
      SDEP = CAPAC(2) * 5.                                              
      FAC = ( SDEP - Z1 ) / ( Z2 - Z1 )                                 
      FAC = AMAX1( 0., FAC )                                            
      FAC = AMIN1( 0.99, FAC )                                          
!                                                                       
      ZAT = ZLT(IVEG) / VCOVER(IVEG)                                    
      IF ( IVEG .EQ. 1 ) ZAT = ZAT * (1.-FAC)                           
!                                                                       
      POWER1 = AMIN1( PSI*ZAT, 50. )                                    
      POWER2 = AMIN1( EXTKB*ZAT, 50. )                                  
      EPSI = EXP( - POWER1 )                                            
      EK = EXP ( - POWER2 )                                             
!                                                                       
      ROSB = SOREF(IWAVE)                                               
      ROSD = SOREF(IWAVE)                                               
      IF ( IVEG .EQ. 2 ) GO TO 300                                      
      ROSB = ALBEDO(2,IWAVE,1)                                          
      ROSD = ALBEDO(2,IWAVE,2)                                          
300   CONTINUE                                                          
!                                                                       
      GE = ROSB / ROSD                                                  
!                                                                       
!-----------------------------------------------------------------------
!     CALCULATION OF DIFFUSE ALBEDOS                                    
!-----------------------------------------------------------------------
!                                                                       
      F1 = BE - CE / ROSD                                               
      ZP = ZMEW * PSI                                                   
!                                                                       
      DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI - &                         
           ( BE - ZP ) * ( F1 + ZP ) * EPSI                            
      ALPHA = CE * ( F1 - ZP ) / EPSI / DEN                             
      BETA = -CE * ( F1 + ZP ) * EPSI / DEN                             
      F1 = BE - CE * ROSD                                               
      DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI                     
!                                                                       
      GAMMA = ( F1 + ZP ) / EPSI / DEN                                  
      DELTA = - ( F1 - ZP ) * EPSI / DEN                                
!                                                                       
      ALBEDO(IVEG,IWAVE,2) =  ALPHA + BETA
!     XQQ(IVEG,IWAVE,2) = ALBEDO(IVEG, IWAVE, 2)                        
!                                                                       
      IF ( IVEG .EQ. 1 ) GO TO 400                                      
      SCOV2 = 0.                                                        
      IF ( TGS .LE. TF ) SCOV2 = AMIN1( 1., CAPAC(2) / 0.004 )          
      ALBEDO(2,IWAVE,2) = &                                              
      ROSD * ( 1. - VCOVER(2) ) + ALBEDO(2,IWAVE,2) * VCOVER(2)        
      ALBEDO(2,IWAVE,2) =                          &			 
      ( 1. - SCOV2 ) * ALBEDO(2,IWAVE,2) + SCOV2 * &	 		
      ( 1.2-IWAVE*0.4 ) *			   &	 		
      FMELT							       
400   CONTINUE                                                          
!                                                                       
      TRANC2(IWAVE) = GAMMA * EPSI + DELTA / EPSI                       
!                                                                       
!-----------------------------------------------------------------------
!     CALCULATION OF DIRECT ALBEDOS                                     
!-----------------------------------------------------------------------
!                                                                       
      F1 = BE - CE / ROSD                                               
      ZMK = ZMEW * EXTKB                                                
!                                                                       
      DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI - &                         
           ( BE - ZP ) * ( F1 + ZP ) * EPSI                            
      ALPHA = ( DE - TORE * ( BE + ZMK ) ) * ( F1 - ZP ) / EPSI - &      
             ( BE - ZP ) * ( DE - CE*GE - TORE * ( F1 + ZMK ) ) * EK   
      ALPHA = ALPHA / DEN                                               
      BETA = ( BE + ZP ) * (DE - CE*GE - TORE * ( F1 + ZMK ))* EK - &    
            ( DE - TORE * ( BE + ZMK ) ) * ( F1 + ZP ) * EPSI          
      BETA = BETA / DEN                                                 
      F1 = BE - CE * ROSD                                               
      DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI                     
      GAMMA = - SIGE * ( F1 + ZP ) / EPSI - &                            
             ( FE + CE * GE * ROSD + SIGE * ( ZMK - F1 ) ) * EK        
      GAMMA = GAMMA / DEN                                               
      DELTA = ( CE * GE * ROSD + FE + SIGE * ( ZMK - F1 ) ) * EK &       
             + SIGE * ( F1 - ZP ) * EPSI                               
      DELTA = DELTA / DEN                                               
!                                                                       
      ALBEDO(IVEG,IWAVE,1) = TORE + ALPHA + BETA                        
!     XQQ(IVEG,IWAVE,1) = ALBEDO(IVEG, IWAVE, 1)                        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      IF( IVEG .EQ. 1 ) GO TO 500                                       
      ALBEDO(2,IWAVE,1) = ROSB * ( 1. - VCOVER(2) ) &                    
                         + ALBEDO(2,IWAVE,1) * VCOVER(2)               
      ALBEDO(2,IWAVE,1) = ( 1. - SCOV2 ) * ALBEDO(2,IWAVE,1) + &         
                         SCOV2 * ( 1.2-IWAVE*0.4 ) * FMELT             
!                                                                       
500   CONTINUE                                                          
!                                                                       
      TRANC1(IWAVE) = EK                                                
      TRANC3(IWAVE) = SIGE * EK + GAMMA * EPSI + DELTA / EPSI           
!                                                                       
2000  CONTINUE                                                          
!                                                                       
!---------------------------------------------------------------------- 
!     CALCULATION OF TERMS WHICH MULTIPLY INCOMING SHORT WAVE FLUXES    
!     TO GIVE ABSORPTION OF RADIATION BY CANOPY AND GROUND              
!---------------------------------------------------------------------- 
!                                                                       
      RADFAC(2,IWAVE,1) = ( 1.-VCOVER(1) ) * ( 1.-ALBEDO(2,IWAVE,1) ) &  
     	    + VCOVER(1) * ( TRANC1(IWAVE) * ( 1.-ALBEDO(2,IWAVE,1) )  & 
     	    + TRANC3(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )	       
!                                                                       
      RADFAC(2,IWAVE,2) = ( 1.-VCOVER(1) ) * ( 1.-ALBEDO(2,IWAVE,2) ) &  
            + VCOVER(1) *  TRANC2(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) )    
!                                                                       
      RADFAC(1,IWAVE,1) = VCOVER(1) * ( ( 1.-ALBEDO(1,IWAVE,1) ) &	 
     	    - TRANC1(IWAVE) * ( 1.-ALBEDO(2,IWAVE,1) )  	 &	
     	    - TRANC3(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )	       
!                                                                       
      RADFAC(1,IWAVE,2) = VCOVER(1) * ( ( 1.-ALBEDO(1,IWAVE,2) ) &       
            - TRANC2(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )               
!                                                                       
!     XQQ(1,IWAVE,1) = RADFAC(1,IWAVE,1)                                
!     XQQ(1,IWAVE,2) = RADFAC(1,IWAVE,2)                                
!     XQQ(2,IWAVE,1) = RADFAC(2,IWAVE,1)                                
!     XQQ(2,IWAVE,2) = RADFAC(2,IWAVE,2)                                
!                                                                       
!                                                                       
!---------------------------------------------------------------------- 
!     CALCULATION OF TOTAL SURFACE ALBEDOS ( SALB )                     
!---------------------------------------------------------------------- 
!                                                                       
      DO 3000 IRAD = 1, 2                                               
      SALB(IWAVE,IRAD) = ( 1.-VCOVER(1) ) * ALBEDO(2,IWAVE,IRAD) + &     
                        VCOVER(1) * ALBEDO(1,IWAVE,IRAD)               
3000  CONTINUE                                                          
!                                                                       
!---------------------------------------------------------------------- 
!     SAVING OF EXTINCTION COEFFICIENTS ( PAR ) FOR STOMAT CALCULATION  
!---------------------------------------------------------------------- 
      IF ( IWAVE .EQ. 2 ) GO TO 600                                     
      RADSAV(1) = 1. - VCOVER(1) + &                                     
                VCOVER(1) * ( TRANC1(IWAVE) + TRANC3(IWAVE) )         
      RADSAV(2) = 1. - VCOVER(1) + VCOVER(1) * TRANC2(IWAVE)            
!     XQQ(1,1,1) = RADSAV(1)                                            
!     XQQ(1,2,1) = RADSAV(2)                                            
600   CONTINUE                                                          
!                                                                       
1000  CONTINUE                                                          
!
!     albedo adjustment ==============================================
!
      if (xadj.eq.0.) go to 730
      xx = radfac(1,1,2) + radsav(2)
      xy = radfac(1,1,1) + radsav(1)
      ssum = salb(1,1)*frac(1,1) + salb(1,2)*frac(1,2)+&
            salb(2,1)*frac(2,1) + salb(2,2)*frac(2,2)
!     for diffuse albedo
      do 650 iwave = 1, 2
      salb(iwave,2) = salb(iwave,2) + xadj * salb(iwave,2) / ssum
      x0 = 1. - salb(iwave,2)
      x1 = radfac(1,iwave,2) + radfac(2,iwave,2)
      x2 = radfac(1,iwave,2) / x1
      x3 = radfac(2,iwave,2) / x1
      radfac(1,iwave,2) = x0 * x2
      radfac(2,iwave,2) = x0 * x3
      if (salb(iwave,2).gt.1..or.radfac(1,iwave,2).gt.1..or. &
     	 radfac(2,iwave,2).gt.1..or.salb(iwave,2).lt.0..or.  &
     	 radfac(1,iwave,2).lt.0..or.radfac(2,iwave,2).lt.0.) then
          stop 999
      end if
 650  continue
 640  format(1x,'unrealistic value, dif',2i12,4e11.4)
!     for direct albedo
      do 750 iwave = 1, 2
      salb(iwave,1) = salb(iwave,1) + xadj * salb(iwave,1) / ssum
      x0 = 1. - salb(iwave,1)
      x1 = radfac(1,iwave,1) + radfac(2,iwave,1)
      x2 = radfac(1,iwave,1) / x1
      x3 = radfac(2,iwave,1) / x1
      radfac(1,iwave,1) = x0 * x2
      radfac(2,iwave,1) = x0 * x3
      radsav(1) =  xy - radfac(1,1,1)
      radsav(2) =  xx - radfac(1,1,2)
      if (salb(iwave,1).gt.1..or.radfac(1,iwave,1).gt.1..or. &
     	 radfac(2,iwave,1).gt.1..or.salb(iwave,1).lt.0..or.  &
     	 radfac(1,iwave,1).lt.0..or.radfac(2,iwave,1).lt.0.) then
          stop 999
      end if
 750  continue
 740  format(1x,'unrealistic value',2i12,4e11.4)
 730  continue
!**************** end adjustment *******************************
!      sibsu = radn(1,1)*salb(1,1) + radn(1,2)*salb(1,2)&
!                          + radn(2,1)*salb(2,1) + radn(2,2)*salb(2,2)
!      if ((swdown.gt.0.1).and.(sibsu.gt.0.1)) then
!         bedo = sibsu / swdown
!         if (bedo.gt.1.) then
!            sibsu =  0.
!            bedo = 999.
!         endif
!      else
!         sibsu = 0.0
!         bedo = 999.
!      endif

!---------------------------------------------------------------------- 
!                                                                       
!     CALCULATION OF LONG-WAVE FLUX TERMS FROM CANOPY AND GROUND        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      TC4 = TC * TC * TC * TC                                           
      TG4 = TGS * TGS * TGS * TGS                                       
!                                                                       
      ZKAT = EXTK(1,3,2) * ZLT(1) / VCOVER(1)                           
      ZKAT = AMIN1( 50. , ZKAT )                                        
      ZKAT = AMAX1( 1.E-5, ZKAT )                                       
      THERMK = EXP(-ZKAT)                                               
!                                                                       
      FAC1 =  VCOVER(1) * ( 1.-THERMK )                                 
      FAC2 =  1.                                                        
      CLOSS =  2. * FAC1 * STEFAN * TC4                                 
      CLOSS =  CLOSS - FAC2 * FAC1 * STEFAN * TG4                       
      GLOSS =  FAC2 * STEFAN * TG4                                      
      GLOSS =  GLOSS - FAC1 * FAC2 * STEFAN * TC4                       
!                                                                       
      ZLWUP =  FAC1 * STEFAN * TC4 + (1. - FAC1 ) * FAC2 * STEFAN * TG4 
      TGEFF = SQRT( SQRT ( ( ZLWUP / STEFAN ) ) )                       
!                                                                       
      RADSAV(3) = EXTK(1,1,1)                                           
      RADSAV(4) = EXTK(1,1,2)                                           
      RADSAV(5) = EXTK(2,1,1)                                           
      RADSAV(6) = EXTK(2,1,2)                                           
      RADSAV(7) = THERMK                                                
      RADSAV(8) = EXTK(1,3,1)                                           
      RADSAV(9) = EXTK(2,3,1)                                           
      RADSAV(10)= CLOSS                                                 
      RADSAV(11)= GLOSS                                                 
      RADSAV(12)= TGEFF                                                 
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      CALL LONGRN( TRANC1, TRANC2, TRANC3)                              
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      CALL RADUSE                                                       
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      RETURN                                                            
      END                                                               
!
!=======================================================================
!                                                                       
      SUBROUTINE RADUSE                                                 
!                                                          1 AUGUST 1988
!=======================================================================
!                                                                       
!     CALCULATION OF ABSORPTION OF RADIATION BY SURFACE                 
!                                                                       
!-----------------------------------------------------------------------
      use ssib_common                                                  
!                                                                       
      P1F         = RADSAV(1)                                           
      P2F         = RADSAV(2)                                           
      EXTK(1,1,1) = RADSAV(3)                                           
      EXTK(1,1,2) = RADSAV(4)                                           
      EXTK(2,1,1) = RADSAV(5)                                           
      EXTK(2,1,2) = RADSAV(6)                                           
      THERMK      = RADSAV(7)                                           
      EXTK(1,3,1) = RADSAV(8)                                           
      EXTK(2,3,1) = RADSAV(9)                                           
      CLOSS       = RADSAV(10)                                          
      GLOSS       = RADSAV(11)                                          
      TGEFF       = RADSAV(12)                                          
!                                                                       
!---------------------------------------------------------------------- 
!     SUMMATION OF SHORT-WAVE RADIATION ABSORBED BY CANOPY AND GROUND   
!---------------------------------------------------------------------- 
!                                                                       
      RADT(1) = 0.                                                      
      RADT(2) = 0.                                                      
!                                                                       
      DO 1000 IVEG  = 1, 2                                              
      DO 1000 IWAVE = 1, 2                                              
      DO 1000 IRAD  = 1, 2                                              
!                                                                       
      RADT(IVEG) = RADT(IVEG)+RADFAC(IVEG,IWAVE,IRAD)*RADN(IWAVE,IRAD)  
!                                                                       
1000  CONTINUE                                                          
!                                                                       
      swcan=radt(1)
      swgnd=radt(2)
!                                                                       
      RADT(1) = RADT(1) + RADN(3,2)*VCOVER(1)*(1.- THERMK)&              
             - CLOSS                                                   
      RADT(2) = RADT(2) + RADN(3,2)*( 1.-VCOVER(1)*(1-THERMK) ) &        
             - GLOSS                                                   
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      PAR(1) = RADN(1,1) + RADN(1,2) + 0.001                            
      PD(1) = ( RADN(1,1) + 0.001 ) / PAR(1)                            
      P1 = P1F * RADN(1,1) + 0.001                                      
      P2 = P2F * RADN(1,2)                                              
      PAR(2) = P1 + P2                                                  
      PD(2) = P1 / PAR(2)                                               
!                                                                       
      RETURN                                                            
      END                                                               
!                                                                       
!=======================================================================
!                                                                       
      SUBROUTINE LONGRN( TRANC1, TRANC2, TRANC3)                        
!                                                          1 AUGUST 1988
!=======================================================================
!                                                                       
!     CALCULATION OF DOWNWARD LONGWAVE.                                 
!                                                                       
!-----------------------------------------------------------------------
      use ssib_common                                                  
!                                                                       
      DIMENSION TRANC1(2), TRANC2(2), TRANC3(2)                         
!                                                                       
      IF(ILW .EQ. 1)GO TO 101                                           
      IF(ILW .EQ. 2)GO TO 102                                           
      IF(ILW .EQ. 3)GO TO 103                                           
101   CONTINUE                                                          
!---------------------------------------------------------------------- 
!     DOWNWARD LONG-WAVE ASSUMED TO BE PROVIDED AS RADN(3,2)            
!---------------------------------------------------------------------- 
      GO TO 200                                                         
!                                                                       
102   CONTINUE                                                          
!---------------------------------------------------------------------- 
!     DOWNWARD LONG-WAVE FROM BRUNT'S EQUATION, MONTEITH(1973), P37.    
!---------------------------------------------------------------------- 
      ESKY = 0.53 + 0.06*SQRT(EM)                                       
      RADN(3,2)  =  ESKY*(1.+0.2*(CLOUD*CLOUD))*STEFAN*TM**4            
      GO TO 200                                                         
!                                                                       
103   CONTINUE                                                          
!---------------------------------------------------------------------- 
!     DOWNWARD LONG-WAVE FLUX CALCULATED AS RESIDUAL FROM MEASURED      
!     NET RADIATION AND OUTGOING LONGWAVE RADIATION.                    
!                                                                       
!     CALCULATION OF ABSORBED FRACTIONS OF RADIATION ( EXPENDABLE )     
!---------------------------------------------------------------------- 
      DO 2000 IWAVE = 1, 2                                              
!                                                                       
      RAB(2,IWAVE,1) =  ( 1. - VCOVER(1) ) *  &                          
       ( RADN(IWAVE,1) * ( 1. - ALBEDO(2,IWAVE,1) ) )                  
      RAB(2,IWAVE,2) =  ( 1. - VCOVER(1) ) *  &                          
         RADN(IWAVE,2) * ( 1. - ALBEDO(2,IWAVE,2) )                    
!                                                                       
      RAB(2,IWAVE,1) = RAB(2,IWAVE,1) + VCOVER(1) *                     &
       ( RADN(IWAVE,1) * ( TRANC1(IWAVE) * ( 1. - ALBEDO(2,IWAVE,1) ) + &
     	 TRANC3(IWAVE) * ( 1. - ALBEDO(2,IWAVE,2) ) ) ) 	     	
      RAB(2,IWAVE,2) = RAB(2,IWAVE,2) + VCOVER(1) * &                   
         RADN(IWAVE,2) * TRANC2(IWAVE) * ( 1. - ALBEDO(2,IWAVE,2) )    
!                                                                     
      RAB(1,IWAVE,1) =  VCOVER(1) *                     &		
     	 RADN(IWAVE,1) * ( ( 1. - ALBEDO(1,IWAVE,1) ) - &	
     	 TRANC1(IWAVE) * ( 1. - ALBEDO(2,IWAVE,1) ) -	&	 
     	 TRANC3(IWAVE) * ( 1. - ALBEDO(2,IWAVE,2) ) )		 
      RAB(1,IWAVE,2) =  VCOVER(1) *                      &	    
     	 RADN(IWAVE,2) * ( ( 1. - ALBEDO(1,IWAVE,2) ) -  &	    
     	 TRANC2(IWAVE) * ( 1. - ALBEDO(2,IWAVE,2) ) )		    
2000  CONTINUE                                                        
!                                                                      
      SWAB = RAB(1,1,1) + RAB(1,1,2) + RAB(1,2,1) + RAB(1,2,2) + &       
            RAB(2,1,1) + RAB(2,1,2) + RAB(2,2,1) + RAB(2,2,2)         
      SWUP = SWDOWN - SWAB                                              
      RADN(3,2) = RNETM - SWAB + ZLWUP                                  
!                                                                       
200   CONTINUE                                                          
!                                                                       
      RETURN                                                           
      END                                                               

!=======================================================================
!                                                                       
      SUBROUTINE ROOT1
!                                                          1 DEC 1988   
!=======================================================================
!                                                                       
!    CALCULATION OF SOIL MOISTURE POTENTIALS IN ROOT ZONE OF EACH       
!    VEGETATION LAYER AND SUMMED SOIL+ROOT RESISTANCE                   
!                                                                       
!-----------------------------------------------------------------------
      use ssib_common                                                   
!                                                                       
      DO 1000 IL = 1, 3                                                 
      PHSOIL(IL) = PHSAT * AMAX1( 0.05, WWW(IL) ) ** ( - BEE )          
 1000 CONTINUE                                                          
!                                                                       
!-----------------------------------------------------------------------
!     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE USED FOR SOURCE      
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!      PHROOT(1) = PHSOIL(1)-0.01                                        
!                                                                       
!      DO 1200 I = 2 ,3                                                  
! 1200 PHROOT(1) = AMAX1( PHROOT(1), PHSOIL(I) )                         
!      PHROOT(2) = PHROOT(1)                                             
!                                                                       
!                                                                       
      RETURN                                                            
      END  

!=======================================================================
!                                                                       
      SUBROUTINE STOMAT1                                                 
!                                                         19 DECEMB 1988
!=======================================================================
!                                                                       
!     CALCULATION OF PAR-LIMITED STOMATAL RESISTANCE                    
!                                                                       
!-----------------------------------------------------------------------
      use ssib_common                                                   
!                                                                       
      DO 1000 IVEG = 1, 2                                               
!                                                                       
      AT = ZLT(IVEG) / VCOVER(IVEG)                                     
!                                                                       
      IF (SUNANG .LE. 0.02) THEN                                        
         XABC = RSTPAR(IVEG,1) / RSTPAR(IVEG,2) + RSTPAR(IVEG,3)        
         RST(IVEG) = 0.5 / XABC * AT                                     
         IF (RST(IVEG) .LT. 0.) RST(IVEG) = 0.00001                     
         GO TO 1010                                                     
      END IF                                                            
!                                                                       
      GAMMA = ( RSTPAR(IVEG,1) + RSTPAR(IVEG,2) *   &  
              RSTPAR(IVEG,3) ) / RSTPAR(IVEG,3)                                          
!                                                                       
      POWER1 = AMIN1( 50., AT * EXTK(IVEG,1,1) )                        
      POWER2 = AMIN1( 50., AT * EXTK(IVEG,1,2) )                        
!                                                                       
!-----------------------------------------------------------------------
!     ROSS INCLINATION FUNCTION                                         
!-----------------------------------------------------------------------
!                                                                       
      AA = 0.5 - 0.633 * CHIL(IVEG)- 0.33 * CHIL(IVEG)* CHIL(IVEG)      
      BB = 0.877 * ( 1. - 2. * AA )                                     
!                                                                       
!-----------------------------------------------------------------------
!     COMBINED ESTIMATE OF K-PAR USING WEIGHTS FOR DIFFERENT COMPONENTS 
!-----------------------------------------------------------------------
!                                                                       
      ZAT = ALOG( ( EXP(-POWER1) + 1. )/2. ) * PD(IVEG) &                
           / ( POWER1/AT )                                             
      ZAT = ZAT + ALOG( ( EXP(-POWER2) + 1. )/2. )   &                   
      * ( 1. - PD(IVEG) ) / ( POWER2/AT )                              
!                                                                       
      POW1 = AMIN1( 50., (POWER1*ZAT/AT) )                              
      POW2 = AMIN1( 50., (POWER2*ZAT/AT) )                              
!                                                                       
      ZK = 1. / ZAT * ALOG( PD(IVEG) * EXP ( POW1 ) &                    
           + ( 1. - PD(IVEG) ) * EXP ( POW2 ) )                        
!                                                                       
!                                                                       
      POW = AMIN1( 50., ZK*AT )                                         
      EKAT = EXP ( POW )                                                
!                                                                       
      AVFLUX = PAR(IVEG) * ( PD(IVEG) / SUNANG * ( AA + BB * SUNANG ) &  
     	   + ( 1. - PD(IVEG) )*( BB / 3. + AA * 1.5		      & 
     	   + BB / 4. * PIE ))					       
!                                                                       
      RHO4 = GAMMA / AVFLUX                                             
!                                                                       
      RST(IVEG) = RSTPAR(IVEG,2)/GAMMA * ALOG(( RHO4 * EKAT + 1. ) / &   
                   ( RHO4 + 1. ) )                                     
      RST(IVEG) = RST(IVEG) - ALOG (( RHO4 + 1. / EKAT ) / &             
                   ( RHO4 + 1. ) )                                     
      RST(IVEG) = RST(IVEG) / ( ZK * RSTPAR(IVEG,3) )                   
!                                                                       
!---------------------------------------------------------------------- 
!     MODIFICATIONS FOR GREEN FRACTION : RST UPRIGHT                    
!---------------------------------------------------------------------- 
!                                                                       
1010  RST(IVEG) = 1. / ( RST(IVEG) * GREEN(IVEG) + 0.0000001)           
1000  CONTINUE                                                          
!                                                                       
!gustavo      rst(1) = rst(1) * ctlpa
      rst(1) = rst(1) * 1
      RETURN                                                            
      END                                                               

!====================================================================== 
!                                                                       
      SUBROUTINE STRES1( IFIRST ,RSTM)                                  
!                                                                       
!====================================================================== 
!                                                                       
!     CALCULATION OF ADJUSTMENT TO LIGHT DEPENDENT STOMATAL RESISTANCE  
!     BY TEMPERATURE, HUMIDITY AND STRESS FACTORS                       
!     SIMPLIFIED SEE XUE ET AL(1991)                                    
!                                                                       
!         RSTFAC(IVEG,1) = FD                                           
!         RSTFAC(IVEG,2) = FP                                           
!         RSTFAC(IVEG,3) = FT                                           
!         RSTFAC(IVEG,4) = FTPD                                         
!                                                                       
!---------------------------------------------------------------------- 
      use ssib_common                                                   
!     REAL RSTM(2), FT(2), DEP(3)                                  
      REAL RSTM(2), DEP(3)                                  
!---------------------------------------------------------------------- 
!     HUMIDITY, TEMPERATURE AND TRANSPIRATION FACTORS                   
!---------------------------------------------------------------------- 
!                                                                       
      DO 1000 IVEG = 1, 2                                               
!                                                                       
      TV = TC                                                           
      ETV = ETC                                                         
      RAIR = RB * 2.                                                    
      IF ( IVEG .EQ. 1 ) GO TO 100                                      
      TV = TGS                                                          
      ETV = ETGS                                                        
      RAIR = RD                                                         
100   CONTINUE                                                          
!                                                                       
      TV = AMIN1 ( ( TU(IVEG) - 0.1 ), TV )                             
      TV = AMAX1 ( ( TLL(IVEG) + 0.1 ), TV )                            
!                                                                       
      IF( IFIRST .EQ. 0 ) GO TO 200                                     
      RSTM(IVEG) = RST(IVEG)                                            
      D2 = ( TU(IVEG) - TOPT(IVEG) ) / ( TOPT(IVEG) - TLL(IVEG) )       
      D1 = 1. /(( TOPT(IVEG) - TLL(IVEG) )*    &                         
             EXP( ALOG( TU(IVEG) - TOPT(IVEG))*D2))                    
      RSTFAC(IVEG,3) = D1*( TV-TLL(IVEG)) * EXP(ALOG(TU(IVEG)-TV)*D2)   
!                                                                       
      IF (RSTFAC(IVEG,3).LT.0.) RSTFAC(IVEG,3) = 0.                     
      IF (RSTFAC(IVEG,3).GT.1.) RSTFAC(IVEG,3) = 1.                     
!                                                                       
!                                                                       
!---------------------------------------------------------------------- 
!      SIMPLIFIED CALCULATION OF LEAF WATER POTENTIAL FACTOR , FP       
!---------------------------------------------------------------------- 
!                                                                       
! Gustavo, what does this code do?  nroot and rootp never defined!
!      if (nroot.eq.1) then
      XROT = ROOTD(1)                                                   
      DO 7400 I = 1, 3                                                  
 7400 DEP(I) = 0.                                                       
      DO 7500 I = 1, 3                                                  
      DEP(I) = MIN(ZDEPTH(I), XROT)                                     
      XROT = XROT - ZDEPTH(I)                                           
      IF (XROT.LE.0.) GO TO 7410                                        
 7500 CONTINUE                                                          
 7410 CONTINUE                                                          
      XDR = (PHSOIL(1) * DEP(1) + PHSOIL(2) * DEP(2)  &                  
           +PHSOIL(3) * DEP(3)) /ROOTD(1)                              
! Gustavo, what does this code do?  nroot and rootp never defined!
!      else
!      XDR = PHSOIL(1) * rootp(1) + PHSOIL(2) * rootp(2)  &                  
!           +PHSOIL(3) * rootp(3)                             
!      end if
      XDR = - XDR                                                       
      IF (XDR .LE. 0.001) XDR = 0.001                                   
      XDR = ALOG (XDR)                                                  
      RSTFAC(IVEG,2) = 1. - EXP(- PH1(IVEG) * (PH2(IVEG) - XDR))       
      IF (RSTFAC(IVEG,2).GT.1.) RSTFAC(IVEG,2) = 1.                     
      IF (RSTFAC(IVEG,2).LT.0.) RSTFAC(IVEG,2) = 0.                     
!                                                                       
200   RST(IVEG) = RSTM(IVEG)                                            
!                                                                       
      EPOT = ETV - EA                                                   
      EPOT = AMAX1(0.0001,(ETV-EA)) 
!
!
!               ***** PJS mod 10/9/92 ***** 
! ***** based on Verma FIFE-87 function for C4 grasses *****
!                                    
!     RSTFAC(IVEG,1) = 1. - DROP * DEFAC(IVEG)
!
      rstfac(iveg,1) = 1./ ( 1 + defac(iveg)*drop )                          
!
      IF (RSTFAC(IVEG,1).LT.0.) RSTFAC(IVEG,1) = 0.                     
      IF (RSTFAC(IVEG,1).GT.1.) RSTFAC(IVEG,1) = 1.                     
!                                                                       
!                                                                       
!---------------------------------------------------------------------- 
!     VALUE OF FP FOUND                                                 
!---------------------------------------------------------------------- 
!                                                                       
300   FTPD = RSTFAC(IVEG,1) * RSTFAC(IVEG,2) * RSTFAC(IVEG,3)           
      RSTFAC(IVEG,4) = AMAX1( FTPD, 0.00001 )                           
!---------------------------------------------------------------------- 
!                                                                       
      RST(IVEG) = RST(IVEG) / RSTFAC(IVEG,4) / VCOVER(IVEG)             
!                                                                       
      RST(IVEG) = AMIN1( RST(IVEG), 100000. )                           
!                                                                       
1000  CONTINUE                                                          
!                                                                       
      RETURN                                                            
      END                                                               
                                                                        
!=======================================================================
      subroutine zenithh(zlat,zlon,dtt,year,month,day,sec,zenang)

      implicit none
!***********************************************************************
!***  compute cosine solar zenith angle averaged over timestep
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! nymd,nhms    (-)       date/time
! lonco        (deg)     longitude (-180 west to 180 east)
! latco        (deg)     latitude (-90 south to 90 north)
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! cosz         (radians) cosine of solar zenith angle
!-----------------------------------------------------------------------
      integer istrip, nymd, nhms, icnt, numstep
      integer n, nnn, sec1, sec2, imadd
      integer year, month, day, sec
      real zlat, zlon, zenang, dtt
      real cosz1, cosz2, sum

      icnt = 0
      sum = 0.

      numstep = nint(dtt / 120.0)
      do nnn = 1,numstep
         imadd = nnn - 1
         sec1 = sec + (120 * imadd)
         sec2 = sec + (120 * imadd) + 60
         call astro_ssib(zlat,zlon,year,month,day,sec1,cosz1)
         call astro_ssib(zlat,zlon,year,month,day,sec2,cosz2)

         sum = sum + max(cosz1,1e-10) + max(cosz2,1e-10)
         icnt = icnt + 2
      enddo

      if (icnt.gt.0) then
         zenang = sum / icnt
      else
         zenang = -999.
      endif

      return
      end

!***********************************************************************
      subroutine astro_ssib(zlat,zlon,year,month,day,sec,zenang)
!***********************************************************************
!***  compute cosine solar zenith angle at individual unique time
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! zlat        (deg)     latitude (-90 south to 90 north)
! zlon        (deg)     longitude (-180 west to 180 east)
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! zenang       (radians) cosine of solar zenith angle
!-----------------------------------------------------------------------
      integer istrip, nymd, nhms, i, km, k, kp, iday, idayp1
      integer year, month, day, sec, n, dayscy
      real pi, zero, one, two, six, dg2rd, eqnx, ob
      real daylen, fac, thm, thp, thnow, zs, zc, sj, cj
      real zenang, zlat, zlon
      parameter(pi = 3.1415926535898)
      parameter(eqnx = 80.9028, dg2rd = pi/180., daylen = 86400.)
      parameter(dayscy = 365*4+1, ob = 23.45*dg2rd)
      parameter(zero = 0., one = 1., two = 2., six = 6.)
      real th(dayscy), t0, t1, t2, t3, t4, hc, mndy(12,4)
      data mndy /0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305,     &
                 335, 366, 397, 34*0/
      logical first
      data first /.true./

! Compute day-angles for 4-year cycle
      if (first) then
         print *,'doing SSiB astro for the first time'
         do i = 15,48
            mndy(i,1) = mndy(i-12,1) + 365
         enddo

         km = int(eqnx) + 1
         fac = km - eqnx
         t0 = zero
         t1 = fun_astro(t0) * fac
         t2 = fun_astro(zero + (t1 / two)) * fac
         t3 = fun_astro(zero + (t2 / two)) * fac
         t4 = fun_astro(zero + t3) * fac
         th(km) = (t1 + two * (t2 + t3) + t4) / six

         do k = 2,dayscy
            t1 = fun_astro(th(km))
            t2 = fun_astro(th(km) + (t1 / two))
            t3 = fun_astro(th(km) + (t2 / two))
            t4 = fun_astro(th(km) + t3)
            kp = mod(km,dayscy) + 1
            th(kp) = th(km) + (t1 + two * (t2 + t3) + t4) / six
            km = kp
         enddo

         first = .false.
      endif

! Compute earth-sun distance to current second
      iday = day + mndy(month,mod(year, 4)+1)
      idayp1 = mod(iday, dayscy) + 1
      thm = mod(th(iday), two*pi)
      thp = mod(th(idayp1), two*pi)

      if (thp.lt.thm) thp = thp + two * pi
      fac = float(sec) / daylen
      thnow = thm * (one - fac) + thp * fac

      zs = sin(thnow) * sin(ob)
      zc = sqrt(one - zs * zs)

! Compute cosine of the zenith angle
      fac = fac * two * pi + pi

      sj = sin(zlat * dg2rd)
      cj = sqrt(one - sj * sj)
      hc = cos(fac + zlon * dg2rd)
      zenang = sj * zs + cj * zc * hc
      if (zenang.lt.zero) zenang = zero

      return
      end

!=============================================                                                                           
!=== Begins functions area where all the old F77 functions defined inside the code
! are taken out
! Gustavo 11/25/03
!=============================================                                                                           
     FUNCTION SSIB_E(x)
       real x,e
       SSIB_E = EXP( 21.18123 - 5418. / X ) / .622      
     END FUNCTION                   
     FUNCTION SSIB_GE(x)
       real x,ge
       SSIB_GE = EXP( 21.18123 - 5418. / X ) * 5418. / (X*X) / .622  
     END FUNCTION                   
     FUNCTION SSIB_FS(x)
       REAL X
       SSIB_FS = 66.85 * X                                              
     END FUNCTION                   
     FUNCTION SSIB_FT(x)
       REAL X
       SSIB_FT = 0.904 * X                                          
     END FUNCTION                   
     FUNCTION SSIB_FV(x)
       REAL X
       SSIB_FV = 0.315 * X                                              
     END FUNCTION                   
!=======================================================================
!                                                                       
      SUBROUTINE CROPS(XLAT,MONTH,DAY,NSX,XCOVER)
!,CHIL,ZLT,GREEN,XCOVER,RSTPAR,TOPT,TL,TU,DEFAC,PH2,PH1)
!
!=======================================================================
!
!     A NEW CROP VERSION BY XUE.                            AUG., 1998
!
!     XLAT IS FROM -90 TO 90 DEGREES  FROM S. TO N.
!
      use ssib_common
      DIMENSION XCOVER(2)
!      DIMENSION GREEN (2),XCOVER(2),&
!     	       CHIL  (2),ZLT   (2),&
!     	       RSTPAR(2,3), TOPT(2),&
!     	       TL(2), TU(2), DEFAC(2),&
!     	       PH1(2), PH2(2)
!
      DIMENSION PHENST(9),WLAI(9),WGRN(9)
!
!C-----------------------------------------------------------------
!**               E   J    H    SD   R   HRV  CUT  PRE-E   E
!     SAVE WLAI,WGRN,IHEAD,IEND,DEND,IWHEAT,SYR
      DATA WLAI/1.0, 2.0, 6.0, 4.0, 3.0, 1.0, 0.01, 0.01, 1.0/
      DATA WGRN/0.6, 0.9, 0.8, 0.5, 0.2, 0.1, 0.01, 0.01, 0.6/
      DATA IHEAD,IEND,DEND,IWHEAT/3,9,244.,12/,SYR/365.25E0/
      NSX = 0
      IF (XLAT.LT.0.) THEN
      RDAY= DAY+184
      IF (RDAY.GT.365) RDAY=RDAY-365
      ELSE
      RDAY= DAY
      END IF
      JULDAY=INT(RDAY+0.2)
      PHI=XLAT
!     PHI=90.0    -180.0 E0/PIE * XLAT
      APHI = ABS(PHI)
      IF (APHI.GT.55.) PHI=SIGN(55.,PHI)
      IF (APHI.LT.20.) PHI=SIGN(20.,PHI)
!
      FLIP =   0.0
!     IF(PHI.LT.0.0 E0 )FLIP = 182.5
!
! ** DETERMINE WHEAT PHENOLOGY FOR LATITUDE AND JULIAN DAY
       PHENST(2) = 4.50    *ABS(PHI) - 64.0     + FLIP
       PHENST(3) = 4.74    *ABS(PHI) - 46.2     + FLIP
       PHENST(4) = 4.86    *ABS(PHI) - 30.8     + FLIP
       PHENST(5) = 4.55    *ABS(PHI) -  3.0     + FLIP
       PHENST(6) = 4.35    *ABS(PHI) + 11.5     + FLIP
       PHENST(7) = PHENST(6) + 3.0
       DEMG      = ABS( 5.21    *ABS(PHI) - 0.3    )
       PHENST(1) = PHENST(2) - DEMG
       PHENST(9) = PHENST(1)
       PHENST(8) = PHENST(9) - 5.0
!
       DO 10 NS = 1,9
       IF(PHENST(NS) .LT. 0.0E0)PHENST(NS) = PHENST(NS) + 365.
       IF(PHENST(NS) .GT. 365. )PHENST(NS) = PHENST(NS) - 365.
   10  CONTINUE
!
       ROOTGC = 1.0
       CHILW  =-0.02
       TLAI   = 0.5
       GRLF   = 0.6
!
! ** FIND GROWTH STAGE GIVEN LATITUDE AND DAY
       DO 50 NS = 1,8
       TOP = PHENST(NS+1)
       BOT = PHENST(NS)
       DIFF1 = TOP-BOT
       DIFF2 = RDAY-BOT
       IF(RDAY.GE. BOT .AND. RDAY .LE. TOP ) GO TO 40
       IF(BOT .LT. TOP ) GO TO 50
!
! ** PHENOLOGY STAGES OVERLAP THE END OF YEAR?
       ICOND = 0
       IF(RDAY .GE. BOT   .AND. RDAY .LE. 365.) ICOND = 1
       IF(RDAY .GE. 0.0   .AND. RDAY .LE. TOP ) ICOND = 2
!
       IF(ICOND .EQ. 0)GO TO 50
       IF(ICOND .EQ. 2)GO TO 35
           DIFF1 = 365.    - BOT + TOP
           DIFF2 = RDAY     - BOT
           GO TO 40
!
   35  CONTINUE
           DIFF1 = 365.   - BOT + TOP
           DIFF2 = 365.   - BOT + RDAY
!
! ** DATE FOUND IN PHENOLOGY STAGE
   40  CONTINUE
       IF ((RDAY.GT.PHENST(IHEAD)).AND.(RDAY.LE.DEND)) THEN      
           TLAI=WLAI(IHEAD)
           GRLF=WGRN(IHEAD)
           GO TO 77
       END IF
       IF ((RDAY.GT.DEND).AND.(RDAY.LE.PHENST(IEND))) THEN
          DIFF1=PHENST(IEND)-DEND
          DIFF2=RDAY-DEND
          PERC =  DIFF2/DIFF1
          TLAI =  PERC*(WLAI(IEND)-WLAI(IHEAD)) + WLAI(IHEAD)
          GRLF =  PERC*(WGRN(IEND)-WGRN(IHEAD)) + WGRN(IHEAD)
          GO TO 77
       END IF
       PERC =  DIFF2/DIFF1
       TLAI =  PERC*(WLAI(NS+1)-WLAI(NS)) + WLAI(NS)
       GRLF =  PERC*(WGRN(NS+1)-WGRN(NS)) + WGRN(NS)
   77  CONTINUE
       GO TO  95
   50  CONTINUE
   95  CONTINUE
!xx    XCOVER(IWHEAT,MONTH,1)=0.90*(1.0 - EXP(-TLAI))
       XCOVER(1)=0.90*(1.0 - EXP(-TLAI))
!xx
       ZLTGMX = WLAI(IHEAD)
       ROOTGC = 2910.0    * (0.5    +0.5    *TLAI/ZLTGMX * GRLF)
       IF (NS.NE.1.AND.NS.NE.2) CHILW=-0.2
!
!xx    ZLT   (IWHEAT,MONTH,1) = TLAI
!xx    GREEN (IWHEAT,MONTH,1) = GRLF
!xx    CHIL  (IWHEAT,1) = CHILW
       ZLT   (1) = TLAI
       GREEN (1) = GRLF
       CHIL  (1) = CHILW
!xx
!      ROOTL (IWHEAT,MONTH,1) = ROOTGC
!xx    TOPT  (1) = YOPT(2)
!xx    TL    (1) = YLL (2)
!xx    TU    (1) = YU  (2)
!xx    DEFAC (1) = YEFAC(2)
!xx    PH1   (1) = YH1(2)
!xx    PH2   (1) = YH2(2)
!xx    DO 155 NN = 1, 3
!x155  RSTPAR(1,NN) = YSTPAR(2,NN)
!xx
       NSX = NS
       IF (NSX.EQ.9) NSX = 1
       IF (NSX.GT.6) NSX = 6
!
       RETURN
       END

!======================================================================
!                                                                       
      SUBROUTINE TEMRS1(prin)
!                                                          11 AUG 2000 
!=======================================================================
!     A MODIFIED SIMPLIFIED VERSION (XUE ET AL. 1991)
!     FLUX COUPLING
!     CORE ROUTINE: CALCULATION OF CANOPY AND GROUND TEMPERATURE        
!     INCREMENTS OVER TIME STEP, FLUXES DERIVED.                        
!-----------------------------------------------------------------------
!                                                                       
!     SUBROUTINES IN THIS SUBROUTINE:
!                                 STRES1                                
!                                 STOMAT1
!                                 NEWTON
      use ssib_common
      logical prin
!      DIMENSION ZINC(3), A2(3), Y1(3)  , ITEX(3)
!      DIMENSION RSTM(2)                                                 
!     COMMON/TONNEW/ ZINC(3), A2(3), Y1(3)                              
!     COMMON/NEWT/ ITEX(3)                                              
!                                                                       
!---------------------------------------------------------------------- 
!     E(X) IS VAPOUR PRESSURE IN MBARS AS A FUNCTION OF TEMPERATURE     
!     GE(X) IS D E(X) / D ( TEMP )                                      
!---------------------------------------------------------------------- 
!                                                                       
!      DIMENSION WWW(3), CAPAC(2), SATCAP(2), ZDEPTH(3)
!      DIMENSION VCOVER(2), ZLT(2), RADT(2),ALBEDO(2,3,2)
!      DIMENSION TOPT(2), TL(2), TU(2), DEFAC(2)
!      DIMENSION PH1(2), PH2(2), RST(2), RSTFAC(2,4)
!      DIMENSION ROOTD(2), ROOTP(3), PHSOIL(3)
!
!      E(X) = EXP( 21.18123 - 5418. / X ) / .622                         
!      GE(X) = EXP( 21.18123 - 5418. / X ) * 5418. / (X*X) / .622 

      ETC   = SSIB_E(TC)                                                     
      ETGS  = SSIB_E(TGS)                                                    
      GETC  = SSIB_GE(TC)                                                    
      GETGS = SSIB_GE(TGS)                                                   

      RCP      = RHOAIR * CPAIR                                              
      PSY      = CPAIR / HLAT * PR / 100. / .622
      BPS      = (pr / 1000.) ** (GASR / CPAIR)

      WC = AMIN1( 1., CAPAC(1)/SATCAP(1) )                              
      WG = AMIN1( 1., CAPAC(2)/SATCAP(2) )                              
!                                                                       
!---------------------------------------------------------------------- 
!      RSOIL FUNCTION FROM FIT TO CAMILLO AND GURNEY (1984) DATA.       
!      WETNESS OF UPPER 0.5 CM OF SOIL CALCULATED FROM APPROXIMATION    
!      TO MILLY FLOW EQUATION WITH REDUCED (1/50 ) CONDUCTIVITY IN      
!      TOP LAYER.                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     WT = WWW(1) + 0.75 * ZDEPTH(1) / ( ZDEPTH(1) + ZDEPTH(2) )        
!    &     * (WWW(1) - (WWW(2)**2)/WWW(1) ) / 2. * 50.                  
!     FAC = AMIN1( WT, 0.99 )                                           
!     FAC = AMAX1( FAC, WWW(1) * 0.1 )                            
!
!------------------------------------------------------------
!     Y.K. Xue changed Jan.18,1994
!------------------------------------------------------------
!                                  
      FAC = AMIN1( www(1), 0.99 )                                           
      FAC = AMAX1( FAC, 0.02 )                            
      RSOIL =  101840. * (1. - fac ** 0.0027) 
! 
!------------------------------------------------------------  
!                                                                       
      FAC = AMIN1( WWW(1), 1. )
      FAC = AMAX1( FAC, 0.02  )  
      PSIT = PHSAT * FAC ** (- BEE )                                    
      ARGG = AMAX1(-10.,(PSIT*GRAV/461.5/TGS))                             
      HR = EXP(ARGG)                                                    
      pilphr=hr
!                                                                       
!---------------------------------------------------------------------- 
!     ALTERATION OF AERODYNAMIC TRANSFER PROPERTIES IN CASE OF SNOW     
!     ACCUMULATION.                                                     
!---------------------------------------------------------------------- 
!                                                                       
      RESD = XDD                                                          
      RESZ0 = Z0                                                        
      RESRDC = RDC                                                      
      RESRBC = RBC                                                      
      RESV2 = VCOVER(2)                                                 
!                                                                       
      IF ( TGS .GT. TF ) GO TO 100                                      
!                                                                       
      SDEP = CAPAC(2) * 5.                                              
      SDEP = AMIN1( SDEP, (Z2*0.95) )                                   
      XDD = Z2 - ( Z2-XDD ) / Z2 * ( Z2 - SDEP )                            
      Z0 = Z0 / ( Z2-RESD ) * ( Z2-XDD )                                  
      RDC = RDC * ( Z2-SDEP ) / Z2                                      
      RBC = RBC * Z2 / ( Z2-SDEP )                                      
      VCOVER(2) = 1.                                                    
      WG = AMIN1( 1., CAPAC(2) / 0.004 )                                
      RST(2) = RSOIL                                                    
!                                                                       
100   CONTINUE                                                          
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!      CALCULATION OF EA, TA, RA, RB, RD AND SOIL MOISTURE STRESS       
!      FOR THE BEGINNING OF THE TIME STEP                               
!                                                                       
!---------------------------------------------------------------------- 
      IFIRST = 1                                                        
      ICOUNT = 0                                                        
      IONCE = 1                                                         
!                                                                       
!     TA = TGS                                                          
      TGEN = TGS
      TCEN = TC
      FC = 1.
      FG = 1.
!xx   TA = TM                                                           
      TRIB = TA
!xx   TRIB = TM                                                         
      EA = EM                                                           
      HT = 0.                                                           
      IONCE = 0                                                         
!                                                                       
1000  CONTINUE                                                          
      ICOUNT = ICOUNT + 1                                               
      
      CALL RASIT52(TRIB)

!     *******   IF ( IFIRST .EQ. 1 ) CALL RBRD1 ******
      IF ( IFIRST .EQ. 1 ) THEN
          TCTA = TC - TA
          RB  = 1.0/(SQRT(U2)/RBC+ZLT(1)*.004)

          X1 = TEMDIF
!
          TGTA = TGS- TA
          TEMDIF = ( TGTA + SQRT(TGTA*TGTA) ) / 2. + 0.1
          FIH = SQRT( 1. + 9.*GRAV*TEMDIF*Z2/TGS/( U2*U2) )
          if (prin) then
             print *,'rdc: ',rdc
             print *,'u2: ',u2
             print *,'fih: ',fih
          endif
          RD  = RDC / U2 / FIH
      END IF
!     ******    END RBRD    ********** 
      D1 = 1./RA + 1./RB + 1./RD                                        
!xx   TA = ( TGS/RD + TC/RB + TM/RA ) / D1                              
      TA = ( TGS/RD + TC/RB + TM/RA *bps) / D1                              
!xx
      HT = ( TA - TM ) * RCP / RA                                       
      RCC = RST(1)*FC + 2. * RB                                         
      COC = (1.-WC)/RCC + WC/(2.*RB)                                    
      RG = RST(2)*FG
      RSURF = RSOIL*FG
      COG1 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)*HR &    
            + VCOVER(2)/(RSURF+RD+44.)*HR
      COG2 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD) &       
            + VCOVER(2)/(RSURF+RD+44.)
      COG1 = COG1 + WG/RD * VCOVER(2) 
      COG2 = COG2 + WG/RD * VCOVER(2)
      D2 = 1./RA + COC + COG2 
      TOP = COC * ETC + COG1 * ETGS + EM / RA                           
      EA = TOP / D2
      if (prin) then
         print *,'d2: ',d2
         print *,'top: ',top
         print *,'coc: ',coc
         print *,'cog1: ',cog1
         print *,'cog2: ',cog2
         print *,'etc: ',etc
         print *,'etgs: ',etgs
         print *,'em: ',em
         print *,'ra: ',ra
         print *,'wg: ',wg
         print *,'rd: ',rd
         print *,'rg: ',rg
         print *,'rsurf: ',rsurf
         print *,'hr: ',hr
         print *,'vcover_2: ',vcover(2)
         print *,'fg: ',fg
         print *,'rsoil: ',rsoil
         print *,'rst_2: ',rst(2)
         print *,'------------'
      endif
      DROP = AMAX1( 0., (E(TA)-EA) )
!---------------------------------------------------------------------- 
!
        CALL STRES1 (IFIRST, RSTM)

!---------------------------------------------------------------------- 
!                                                                       
      IFIRST = 0                                                        
      ERIB = EA                                                         
      TRIB = TA                                                         
!CC                                                                     
      IF ( ICOUNT .LE. 4 ) GO TO 1000                                   
!                                                                       
!---------------------------------------------------------------------- 
!
      TC3 = TC * TC * TC
      if (prin) then
         print *,'tc: ',tc
      endif
      TG3 = TGS * TGS * TGS
      FAC1 = ( 1. - ALBEDO(1,3,2) ) * ( 1.-THERMK ) * VCOVER(1)
      FAC2 =   1. - ALBEDO(2,3,2)
      RNCDTC = - 2. * 4. * FAC1 * STEFAN * TC3
      RNCDTG = 4. * FAC1 * FAC2 * STEFAN * TG3
      RNGDTG = - 4. * FAC2 * STEFAN * TG3
      RNGDTC = 4. * FAC1 * FAC2 * STEFAN * TC3

      if (prin) then
         print *,'fac1: ',fac1
         print *,'rcdtc: ',rncdtc
         print *,'rcdtg: ',rncdtg
         print *,'rgdtc: ',rngdtc
         print *,'rgdtg: ',rngdtg
      endif

!---------------------------------------------------------------------- 
!                                                                       
!     DEW CALCULATION : DEW CONDITION IS SET AT BEGINNING OF TIME STEP. 
!     IF SURFACE CHANGES STATE DURING TIME STEP, LATENT HEAT FLUX IS    
!     SET TO ZERO.                                                      
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      IF ( EA .GT. ETC ) FC = 0.                                        
      IF ( EA .GT. ETGS) FG = 0.                                        

      if (prin) then
         print *,'ea: ',ea
         print *,'etc: ',etc
         print *,'etgs: ',etgs
      endif

!---------------------------------------------------------------------- 
!                                                                       
!     WET FRACTION EXHAUSTION TEST : IF CAPAC(X) IS EXHAUSTED IN        
!     A TIME STEP, INTERCEPTION LOSS IS LIMITED TO CAPAC(X).            
!                                                                       
!---------------------------------------------------------------------- 
!     START OF NON-NEUTRAL RESISTANCE CALCULATION LOOP                  
!---------------------------------------------------------------------- 
!                                                                       
      I = 0                                                             
!                                                                      
!    ----- INITIALIZE NEWTON-RAPHSON ITERATIVE ROUTINE FOR RASIT 3,5,8  
                    NOX = 0                                             
                 NONPOS = 1                                             
                  IWALK = 0                                             
                     LX = 2                                             
                   FINC = 1.                                            
!                   ITEX(LX) = 0.                                        
                   ZINC(LX) = 0.                                        
                   A2(LX)   = 0.                                       
                   Y1(LX)   = 0.                                        
2000  CONTINUE                                                          
!                                                                       
      CALL RASIT52(TRIB)
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      RCP = RHOAIR * CPAIR
      if (prin) then
         print *,'rcp: ',rcp
         print *,'ros: ',rhoair
         print *,'cpair: ',cpair
      endif
      D1 = 1./RA + 1./RB + 1./RD
!xx   TA = ( TGS/RD + TC/RB + TM/RA ) / D1
      TA = ( TGS/RD + TC/RB + TM/RA *bps) / D1
!xx
!                
      HC = RCP * ( TC - TA ) / RB * DTT
      HG = RCP * ( TGS - TA ) / RD * DTT
      if (prin) then
         print *,'rcp: ',rcp
         print *,'ta: ',ta
         print *,'tgs: ',tgs
         print *,'tc: ',tc
         print *,'thm: ',tm
         print *,'rd: ',rd
         print *,'rb: ',rb
         print *,'ra: ',ra
         print *,'d1: ',d1
         print *,'bps: ',bps
         print *,'dtt: ',dtt
      endif

!----------------------------------------------------------------------
!     N.B. FLUXES EXPRESSED IN JOULES M-2
!----------------------------------------------------------------------
!                
      HCDTC = RCP / RB * ( 1./RA + 1./RD ) / D1
      HCDTG = - RCP / ( RB * RD ) / D1
! FOR TM
      HCDTM = - RCP / ( RB * RA ) / D1 * BPS
!                
      HGDTG = RCP / RD * ( 1./RA + 1./RB ) / D1
      HGDTC = - RCP / ( RD * RB ) / D1
! FOR TM
      HGDTM = - RCP / ( RD * RA ) / D1 *BPS
!                
!     RCP = RHOAIR * CPAIR
!----------------------------------------------------------------------
!     MODIFICATION FOR SOIL DRYNESS : HR = REL. HUMIDITY IN TOP LAYER
!----------------------------------------------------------------------
!                
      HRR = HR   
      IF ( FG .LT. .5 ) HRR = 1.
!                
      RCC = RST(1)*FC + 2. * RB
      COC = (1.-WC)/RCC + WC/(2.*RB)
      RG = RST(2)*FG
      RSURF = RSOIL*FG
      COG1 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)*HRR &
     	  + VCOVER(2)/(RSURF+RD+44.)*HRR			      
     COG2 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)      &
     	  + VCOVER(2)/(RSURF+RD+44.)
      COG1 = COG1 + WG/RD * VCOVER(2)
      COG2 = COG2 + WG/RD * VCOVER(2)
!                
      D2 = 1./RA + COC + COG2
      TOP = COC * ETC + COG1 * ETGS + EM/RA
      EA = TOP / D2
!                
      EC = ( ETC - EA ) * COC * RCP/PSY * DTT
!                
      EG = ( ETGS*COG1 - EA*COG2 ) * RCP/PSY * DTT
!                
      DEADTC = GETC * COC / D2
      DEADTG = GETGS * COG1 / D2
!                
      ECDTC = ( GETC - DEADTC ) * COC * RCP / PSY
      ECDTG = - DEADTG * COC * RCP / PSY
!                
      EGDTG = ( GETGS*COG1 - DEADTG*COG2 ) * RCP / PSY
      EGDTC = - DEADTC * COG2 * RCP / PSY
!   FOR QM
      DEADQM = 0.622 * pr /( (0.622+QM)**2 * RA * D2 )    
      ECDQM =        -DEADQM * COC * RCP / PSY
      EGDQM =        -DEADQM * COG2 * RCP / PSY
!   FOR YPDATING TM AND QM
      AK = 1/ RCP / BPS
      AH = 1/ (HLAT*RHOAIR)
      XXT = -DTT * AK * (HGDTM+HCDTM)
      XXQ = -DTT * AH * (EGDQM+ECDQM)
!             
!---------------------------------------------------------------------- 
!                                                                       
!     CALCULATION OF COEFFICIENTS OF TEMPERATURE TENDENCY EQUATIONS     
!        C - CANOPY                                                     
!        G - GROUND                                                     
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      CCODTC = CCX / DTT - RNCDTC + HCDTC + ECDTC                       
      CCODTG = - RNCDTG + HCDTG + ECDTG                                 
      CCORHS = RADT(1) - ( HC + EC ) / DTT                              
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      GCODTG = CG / DTT + TIMCON*CG*2. - RNGDTG + HGDTG + EGDTG         
      GCODTC = - RNGDTC + HGDTC + EGDTC                                 
      GCORHS = RADT(2) - TIMCON*CG*2. * ( TGS -TD ) - ( HG + EG ) / DTT 
!                                                                       
      DENOM = CCODTC * GCODTG - CCODTG * GCODTC                         
!                                                                       
      DTC = ( CCORHS * GCODTG - CCODTG * GCORHS ) / DENOM               
      DTG = ( CCODTC * GCORHS - CCORHS * GCODTC ) / DENOM               
      if (prin) then
         print *,'dtc: ',dtc
      endif
!                                                                       
!---------------------------------------------------------------------- 
!     CHECK IF INTERCEPTION LOSS TERM HAS EXCEEDED CANOPY STORAGE       
!---------------------------------------------------------------------- 
!                                                                       
      ECPOT = ( (ETC - EA) + (GETC - DEADTC)*DTC - DEADTG*DTG )         
      ECI = ECPOT * WC /(2.*RB) * RCP/PSY * DTT                         
      ECIDIF=AMAX1(0.0,(ECI-CAPAC(1)*1.E3*HLAT))                        
      ECI   =AMIN1(ECI,(    CAPAC(1)*1.E3*HLAT))                        
!                                                                       
      EGPOT = ( (ETGS - EA) + (GETGS - DEADTG)*DTG - DEADTC*DTC )       
      EGI = EGPOT * VCOVER(2) * WG/RD * RCP/PSY * DTT                   
      EGIDIF=AMAX1(0.0,(EGI-CAPAC(2)*1.E3*HLAT))                        
      EGI   =AMIN1(EGI,(    CAPAC(2)*1.E3*HLAT))                        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      TGEN = TGS + DTG                                                  
      TCEN = TC + DTC                                                   
      D1 = 1./RA + 1./RB + 1./RD                                        
!xx   TAEN = ( TGEN / RD + TCEN / RB + TM / RA ) / D1                   
      TAEN = ( TGEN / RD + TCEN / RB + TM / RA *bps) / D1                   
!xx
!                                                                       
      HEND = ( TAEN - TM ) * RCP / RA + (ECIDIF + EGIDIF)/DTT           
      Y= TRIB - TAEN                                                    
      I = I + 1                                                         
      HT   = HEND                                                       
      IF ( I .GT. ITRUNK ) GO TO 200                                    
!                                                                       
      CALL NEWTON(TRIB,Y,FINC,NOX,NONPOS,IWALK,LX)                      
      IF(NOX.NE.1)GO TO 2000                                            
!                                                                       
200   CONTINUE                                                          
!     IQIN = IQIN + I                                                   
!     IF (I.GT.10) IQIN1 = IQIN1 + 1                                    
1010  FORMAT(1X,I3,1X,'TR1B,Y,RA,RIB,EGDF',7E11.4)                      
1011  FORMAT(1X,'HEND,HT,Y,TA,TC,TG,ECDF',8E11.4)                       
1012  FORMAT(5X,I10,I5)                                                 
1014  FORMAT(5X,F12.5)                                                  
!                                                                       
!---------------------------------------------------------------------- 
!     EXIT FROM NON-NEUTRAL CALCULATION                                 
!                                                                       
!     EVAPOTRANSPIRATION FLUXES CALCULATED FIRST ( J M-2 )              
!---------------------------------------------------------------------- 
!                                                                       
      HRR = HR                                                          
      IF ( FG .LT. .5 ) HRR = 1.                                        
      RSURF = RSOIL*FG                                                  
!                                                                       
      COCT = (1.-WC)/RCC                                                
      COGT = VCOVER(2) * (1.-WG)/( RG + RD )                            
      COGS1 = (1.-VCOVER(2)) / ( RD + RSURF ) * HRR   & 		 
     	     + VCOVER(2) / ( RD + RSURF + 44.) * HRR     		
     COGS2 = COGS1 / HRR			         		
!    						         		
     ECT = ECPOT * COCT * RCP/PSY * DTT 	         		
!    						         		
     EGT = EGPOT * COGT * RCP/PSY * DTT 	         		
     EGS = (ETGS + GETGS*DTG ) * COGS1  	      &  		
     	   - ( EA + DEADTG*DTG + DEADTC*DTC ) * COGS2		       
      EGS = EGS * RCP/PSY * DTT                                         
      EGSMAX = WWW(1) / 2. * ZDEPTH(1) * POROS * HLAT * 1000.           
      EGIADD = AMAX1( 0., EGS - EGSMAX )                                
      EGS = AMIN1 ( EGS, EGSMAX )                                       
      EGIDIF = EGIDIF + EGIADD                                          
!                                                                       
!---------------------------------------------------------------------- 
!     SENSIBLE HEAT FLUX CALCULATED WITH LATENT HEAT FLUX CORRECTION    
!---------------------------------------------------------------------- 
      HC = HC + (HCDTC*DTC + HCDTG*DTG)*DTT + ECIDIF                    
      HG = HG + (HGDTC*DTC + HGDTG*DTG)*DTT + EGIDIF                    
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     TEST OF DEW CONDITION. LATENT HEAT FLUXES SET TO ZERO IF SIGN     
!     OF FLUX CHANGES OVER TIME STEP.EXCESS ENERGY DONATED TO SENSIBLE  
!     HEAT FLUX.                                                        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      ECF = SIGN( 1., ECPOT )                                           
      EGF = SIGN( 1., EGPOT )                                           
      DEWC = FC * 2. - 1.                                               
      DEWG = FG * 2. - 1.                                               
!                                                                       
      IF(DEWC*ECF.GT.0.0) GO TO 300                                     
      HC = HC + ECI + ECT                                               
      ECI = 0.                                                          
      ECT = 0.                                                          
300   IF(DEWG*EGF.GT.0.0) GO TO 400                                     
      HG = HG + EGS + EGI + EGT                                         
      EGS = 0.                                                          
      EGI = 0.                                                          
      EGT = 0.                                                          
400   CONTINUE                                                          
!                                                                       
      EC = ECI + ECT                                                    
      EG = EGT + EGS + EGI                                              
!                                                                       
!---------------------------------------------------------------------- 
!     ADJUSTMENT OF TEMPERATURES AND VAPOR PRESSURE , CALCULATION OF    
!     SENSIBLE HEAT FLUXES.                                             
!---------------------------------------------------------------------- 
!                                                                       
      TC  = TCEN                                                        
      TGS = TGEN                                                        
      TA  = TAEN                                                        
      EA = EA + DEADTC*DTC + DEADTG*DTG                                 
!                                                                       
      RADT(1) = RADT(1) + RNCDTC*DTC + RNCDTG*DTG                       
      RADT(2) = RADT(2) + RNGDTC*DTC + RNGDTG*DTG 
!
! ** simulated net all-wave radiation **
!
!     sibnet(nmm,ndd,nhh) = RADT(1) + RADT(2)                   
!                                                                       
      CHF = CCX / DTT * DTC                                             
      SHF = CG / DTT * DTG + TIMCON*CG*2. * ( TGS - TD )                
!                                                                       
      ZLWUP = ZLWUP - RNCDTC * DTC / 2.  &                               
                   - RNGDTG * DTG * (1.-VCOVER(1)*(1.-THERMK) )        
!                                                                       
      IF ( TGS .GT. TF ) GO TO 500                                      
      EGS = EG - EGI                                                    
      EGT = 0.                                                          
  500 CONTINUE                                                          
!                                                                       
      VCOVER(2) = RESV2                                                 
      XDD = RESD                                                          
      Z0 = RESZ0                                                        
      RDC = RESRDC                                                      
      RBC = RESRBC                                                      
!                                                                       
      RETURN                                                            
      END                                                               
!=====================================================================  
!                                                                       
      SUBROUTINE RASIT52(TRIB)
!                                                                       
!=======================================================================
!     CUU AND CTT ARE LINEAR  (A SIMPLIFIED VERSION, XUE ET AL. 1991)  
!                                                                      
      use ssib_common                                                 
!      FS(X) = 66.85 * X                                              
!      FT(X) = 0.904 * X                                          
!      FV(X) = 0.315 * X                                              
!                                                                      
!     CU AND CT ARE THE FRICTION AND HEAT TRANSFER COEFFICIENTS.    
!     CUN AND CTN ARE THE NEUTRAL FRICTION AND HEAT TRANSFER         
!     COEFFICIENTS.                                                     
!                                                                       
      G2= 0.75                                                          
      G3= 0.75                                                          
      Z22 = Z2                                                          
      ZL = Z2 + 11.785 * Z0                                             
!xx
      if(zwind.le.D.or.zl.le.D) D=min(zwind,zl)-0.1
!xx
      Z2 = D + Z0                                                       
      CUNI = ALOG((ZWIND-D)/Z0)/VKC                                     
      IF (ZL.LT.ZWIND) THEN                                             
         XCT1 = ALOG((ZWIND-D)/(ZL-D))                                  
         XCT2 = ALOG((ZL-D)/(Z2-D))                                     
         XCTU2 = ALOG((ZL-D)/(Z22-D))                                   
         CTNI = (XCT1 + G3 * XCT2) / VKC                                
      ELSE                                                              
         XCT2 =  ALOG((ZWIND-D)/(Z2-D))                                 
         XCTU2 =  ALOG((ZWIND-D)/(Z22-D))                               
         CTNI = G3 * XCT2 /VKC                                          
      END IF                                                            
!        NEUTRAL VALUES OF USTAR AND VENTMF                             
!                                                                       
         USTARN=UM/CUNI                                                 
         VENTN =RHOAIR   /CTNI*USTARN                                   
      IF (ZL.LT.ZWIND) THEN                                             
         U2 = UM - 1. / VKC * USTARN * (XCT1 + G2 *  &                   
                XCTU2)                                                 
      ELSE                                                              
         U2 = UM - 1. / VKC * USTARN * G2 * XCTU2                      
      END IF                                                            
!xx
      if(u2.lt.0.01) u2=0.01
!xx
!                                                                       
!     STABILITY BRANCH BASED ON BULK RICHARDSON NUMBER.                 
!                                                                       
!xx
      THM=TM*bps
      THVGM= TRIB-THM
!xx   THVGM= TRIB-TM
      IF (TA.EQ.0.) THVGM = 0.                                          
!xx
      RIB  = -THVGM*GRAV*(ZWIND-D) / (THM*(UM-U2)**2)
!xx   RIB  = -THVGM*GRAV*(ZWIND-D) / (TM*(UM-U2)**2)
      RIB  = MAX(-10.E0 ,RIB)               
      RIB  = MIN( .1643E0 ,RIB)                                        
!                                                                       
!     NON-NEUTRON CORRECTION  (SEE XUE ET AL(1991))                     
      IF(RIB.LT.0.0)THEN                                                
         GRIB = +RIB                                                    
         GRZL = +RIB*(ZL-D)/(ZWIND-D)                                   
         GRZ2 = +RIB*(Z2-D)/(ZWIND-D)                                   
         FVV =  SSIB_FV(GRIB)                                                
         IF (ZL.LT.ZWIND) THEN                                          
             FTT = SSIB_FT(GRIB) + (G3-1.) * SSIB_FT(GRZL) - G3 * SSIB_FT(GRZ2)        
         ELSE                                                           
             FTT = G3*(SSIB_FT(GRIB) - SSIB_FT(GRZ2))                             
         END IF                                                         
         CUI = CUNI + FVV                                               
         CTI = CTNI + FTT                                               
      ELSE                                                              
         RZL = RIB/(ZWIND-D)*(ZL-D)                                     
         RZ2 = RIB/(ZWIND-D)*(Z2-D)                                     
         FVV = SSIB_FS(RIB)                                                  
         IF (ZL.LT.ZWIND) THEN                                          
             FTT = SSIB_FS(RIB) + (G3-1) * SSIB_FS(RZL) - G3 * SSIB_FS(RZ2)            
         ELSE                                                           
             FTT = G3 * (SSIB_FS(RIB) - SSIB_FS(RZ2))                             
         END IF                                                         
 312     CUI = CUNI + FVV                                               
         CTI = CTNI + FTT                                               
      ENDIF                                                             
 310  CONTINUE                                                          
!                                                                       
      USTAR =UM/CUI                                                     
!rr
      XAKMS=USTAR/CUI
!rr
      RAF = CTI / USTAR                                                 
      IF (RAF.LT.0.80) RAF = 0.80                                       
!                                                                       
      RA  = RAF                                                         
!                                                                       
      UEST  = USTAR                                                     
      DRAG = RHOAIR * UEST*UEST                                         
      Z2 = Z22                                                          
 1010 FORMAT(1X,'RIB,CTI,CUI,CTN,CUN',I10,7E10.3)                       
 1011 FORMAT(1X,'RIB,RAF,USTAR,UM,U2',7E10.3)                           
      RETURN                                                            
      END                                                               

!=======================================================================
!                                                                     
      SUBROUTINE UPDAT1(prin)
!                                                         12 AUGUST 2000
!=======================================================================
!                                                                       
!     UPDATING OF SOIL MOISTURE STORES AND INTERCEPTION CAPACITY        
!                                                                       
!-----------------------------------------------------------------------
!      include 'COMSSIB.ssib.com' 
!      DIMENSION WWW(3), CAPAC(2),EF(3),SNOWW(2)
!      DIMENSION ROOTD(2), ZDEPTH(3), ROOTP(3)
      use ssib_common
      logical prin
       real :: temw(3)
       real :: temwp(3)
       real :: temwpp(3)                                
       real :: aaa(2)
       real :: bbb(2)
       real :: ccc(2)
       real :: qqq(2)                            
       real :: ef(3)                            
!                                                                       
!---------------------------------------------------------------------- 
!     EVAPORATION LOSSES ARE EXPRESSED IN J M-2 : WHEN DIVIDED BY       
!     ( HLAT*1000.) LOSS IS IN M M-2                                    
!     MASS TERMS ARE IN KG M-2 DT-1                                     
!---------------------------------------------------------------------- 
!                                                                       
      SNOFAC = HLAT / ( HLAT + SNOMEL /1000. )                          
      FACKS = 1.                                                        
      IF ( (TC-DTC) .LE. TF ) FACKS = SNOFAC                            
      IF ( (ECT+ECI) .GT. 0.) GO TO 100                                 
      ECI = ECT + ECI                                                   
      ECT = 0.                                                          
      FACKS = 1. / FACKS                                                
100   CAPAC(1)=CAPAC(1) - ECI*FACKS/HLAT/1000.                          
!                                                                       
      ECMASS = ( ECT + ECI * FACKS ) / HLAT                             
!                                                                       
      FACKS = 1.                                                        
      IF ( (TGS-DTG) .LE. TF ) FACKS = SNOFAC                           
      IF ( (EGT+EGI) .GT. 0. ) GO TO 200                                
      EGI = EGT + EGI                                                   
      EGT = 0.                                                          
      FACKS = 1. / FACKS                                                
200   CAPAC(2)=CAPAC(2) - EGI*FACKS/HLAT/1000.                          
!                                                                       
      EGMASS = ( EGT + EGS + EGI * FACKS ) / HLAT                       
!                                                                       
      ETMASS = ECMASS + EGMASS                                          
!                                                                       
      HFLUX = ( HC + HG ) / DTT                                         

      if (prin) then
         print *,'etmass: ',(etmass/1000.)
         print *,'capac_1: ',capac(1)
         print *,'capac_2: ',capac(2)
      endif

!---------------------------------------------------------------------- 
!      DUMPING OF SMALL CAPAC VALUES ONTO SOIL SURFACE STORE            
!---------------------------------------------------------------------- 
!                                                                       
      DO 1000 IVEG = 1, 2                                               
      IF ( CAPAC(IVEG) .GT. 0.000001 ) GO TO 300                        
      FILTR = FILTR + CAPAC(IVEG)                
      WWW(1) = WWW(1) + CAPAC(IVEG) / ( POROS*ZDEPTH(1) )               
      CAPAC(IVEG) = 0.                                                  
300   CONTINUE                                                          
1000  CONTINUE                                                          
!---------------------------------------------------------------------- 
!     SNOWMELT / REFREEZE CALCULATION                                   
!---------------------------------------------------------------------- 
!=======================================================================
!                                                                       
!     CALCULATION OF SNOWMELT AND MODIFICATION OF TEMPERATURES          
!     N.B. THIS VERSION DEALS WITH REFREEZING OF WATER                  
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      DO 7000 IVEG = 1, 2                                               
!                                                                       
      CCT = CCX                                                         
      TS = TC                                                           
      DTS = DTC                                                         
      FLUX = CHF                                                        
      IF ( IVEG .EQ. 1 ) GO TO 7100                                      
      CCT = CG                                                          
      TS = TGS                                                          
      DTS = DTG                                                         
      FLUX = CCT * DTG / DTT                                            
                                                                        
7100  CONTINUE                                                          
!                                                                       
      TTA = TS - DTS                                                    
      TTB = TS                                                          
      if (prin) print *,'tta: ',tta
      SNOWW(IVEG) = 0.                                                  
      IF ( TTA .LE. TF ) SNOWW(IVEG) = CAPAC(IVEG)                      
      CAPAC(IVEG) = CAPAC(IVEG) - SNOWW(IVEG)                           
      IF ( TTA .GT. TF .AND. TTB .GT. TF ) GO TO 7200                    
      IF ( TTA .LE. TF .AND. TTB .LE. TF ) GO TO 7200                    
!                                                                       
      DTF = TF - TTA                                                    
      DTIME1 = CCT * DTF /  FLUX                                        
      HF = FLUX*(DTT-DTIME1)                                            
      FCAP = - CAPAC(IVEG)  * SNOMEL                                    
      SPWET = AMIN1( 5. , SNOWW(IVEG) )                                 
      IF ( DTS .GT. 0. ) FCAP =  SPWET * SNOMEL                         
      DTIME2 = FCAP / FLUX                                              
      DTF2 =   FLUX * (DTT-DTIME1-DTIME2)/CCT                           
      TN = TF + DTF2                                                    
      TS = TF - 0.1                                                     
      IF (ABS(HF) .GE.ABS(FCAP) ) TS = TN                               
      CHANGE = HF                                                       
      IF (ABS(CHANGE) .GE.ABS(FCAP) ) CHANGE = FCAP                     
!                                                                       
      CHANGE = CHANGE / SNOMEL                                          
!rr
      IF (CHANGE.GT.0.0) SMELT=CHANGE+SMELT
!rr
      SNOWW(IVEG) = SNOWW(IVEG) - CHANGE                                
      CAPAC(IVEG) = CAPAC(IVEG) + CHANGE                                
!                                                                       
      IF ( IVEG .EQ. 1 ) TC = TS                                        
      IF ( IVEG .EQ. 2 ) TGS = TS                                       
      IF ( SNOWW(IVEG) .LT. 0.00001 ) GO TO 7200                         
      ZMELT = 0.                                                        
!     modified to force water into soil. Xue Feb. 1994
      ZMELT = CAPAC(IVEG)                             
!     IF ( TD .GT. TF ) ZMELT = CAPAC(IVEG)                             
      FILTR =  FILTR+ ZMELT                    
      WWW(1) = WWW(1) + ZMELT / ( POROS * ZDEPTH(1) )                   
!     IF ( TD .LE. TF )  ROFF = ROFF + CAPAC(IVEG)   
      CAPAC(IVEG) = 0.                                                  
7200   CONTINUE                                                          
!                                                                       
      CAPAC(IVEG) = CAPAC(IVEG) + SNOWW(IVEG)                           
!                                                                       
7000  CONTINUE                                                          
!                                                                       
      FLUXEF = SHF - CCT*DTG/DTT                                        
      TD = TD + FLUXEF / ( CG * 2. * SQRT ( PIE*365. ) ) * DTT          
!
!       *** LOAD PILPS DATA
!
!     if (change .gt. 0) snm(istat)=snm(istat)+(change*1000.)
      change=0.0
!                                                                       
!---------------------------------------------------------------------- 
!     BARE SOIL EVAPORATION LOSS                                        
!---------------------------------------------------------------------- 
!                                                                       
      FILTR = FILTR - EGS / HLAT / 1000.       
      WWW(1) = WWW(1) - EGS / HLAT / 1000. / ( POROS * ZDEPTH(1) )      
!                                                                       
!---------------------------------------------------------------------- 
!   EXTRACTION OF TRANSPIRATION LOSS FROM ROOT ZONE                     
!---------------------------------------------------------------------- 
!                                                                       
      DO 2000 IVEG = 1, 2                                               
!                                                                       
      IF ( IVEG .EQ. 1 ) ABSOIL = ECT / HLAT / 1000.                    
      IF ( IVEG .EQ. 2 ) ABSOIL = EGT / HLAT / 1000.                    
!                                                                       
! Gustavo, what does this code do?  nroot and rootp never defined!
!      if (NROOT.EQ.1) then
      EF(2) = 0.                                                        
      EF(3) = 0.                                                        
      TOTDEP = ZDEPTH(1)                                                
!                                                                       
      DO 3000 IL = 2, 3                                                 
      TOTDEP = TOTDEP + ZDEPTH(IL)                                      
!                                                                       
!     DIV = AMAX1 ( 1., ( PHSOIL(IL) - PHL(IVEG) ) )                    
!                                                                       
      IF ( ROOTD(IVEG) .LT. TOTDEP ) GO TO 400                          
!                                                                       
      EF(IL) = ZDEPTH(IL) / ROOTD(IVEG)                                 
      GO TO 500                                                         
!                                                                       
400   CONTINUE                                                          
      EF(IL) = ROOTD(IVEG) - TOTDEP + ZDEPTH(IL)                        
      EF(IL) = EF(IL) / ROOTD(IVEG)                                     
      GO TO 600                                                         
!                                                                       
500   CONTINUE                                                          
3000  CONTINUE                                                          
!                                                                       
600   EFT = EF(2) + EF(3)                                               
!xx
      EFT = MAX(EFT,0.1E-5)
!xx
      EF(2) = EF(2) / EFT                                               
      EF(3) = EF(3) / EFT                                               
!
      DO 4000 IL = 2, 3                                                 
      WWW(IL) = WWW(IL) - ABSOIL * EF(IL) / ( POROS * ZDEPTH(IL) )      
4000  CONTINUE                                                          
! Gustavo, what does this code do?  nroot and rootp never defined!
!      else
!      ef(1) = rootp(1)
!      ef(2) = rootp(2)
!      ef(3) = rootp(3)
!      DO 4004 IL = 1, 3                                                 
!      WWW(IL) = WWW(IL) - ABSOIL * EF(IL) / ( POROS * ZDEPTH(IL) )      
!4004  CONTINUE                                                          
!      end if
!                                                                       
2000  CONTINUE                                                          
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     CALCULATION OF INTERFLOW, INFILTRATION EXCESS AND LOSS TO         
!     GROUNDWATER .  ALL LOSSES ARE ASSIGNED TO VARIABLE 'ROFF' .       
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      DO 5000 IL = 1, 2                                                 
      IF ( WWW(IL) .GT. 0. ) GO TO 700                                  
      WWW(IL+1) = WWW(IL+1) + WWW(IL) * ZDEPTH(IL)/ZDEPTH(IL+1)         
      WWW(IL) = 0.                                                      
700   CONTINUE                                                          
5000  CONTINUE                                                          
!                                                                       
!     IF ( TD .LT. TF ) GO TO 800                                       
!  
!=======================================================================        
!    calculation of interflow, infiltration excess and loss to                  
!    groundwater .  all losses are assigned to variable 'roff' .                
!----------------------------------------------------------------------         
!                                                                               
      do 8000 i = 1, 3                                                          
!                                                                               
      temw(i)   = amax1( 0.03, www(i) )                                         
      temwp(i)  = temw(i) ** ( -bee )                                           
      temwpp(i) = amin1( 1., temw(i)) ** ( 2.*bee+ 3. )                         
8000  continue                                                                  
!                                                                               
!-----------------------------------------------------------------------        
!                                                                               
!    calculation of gravitationally driven drainage from w(3) : taken           
!    as an integral of time varying conductivity.addition of liston             
!    baseflow term to original q3g to insure flow in                            
!    dry season. modified liston baseflow constant scaled                       
!    by available water.                                                        
!                                                                               
!     q3g (q3) : equation (62) , SE-86                                          
!                                                                               
!-----------------------------------------------------------------------        
!                                                                               
      pows = 2.*bee+2.                                                          
      q3g = temw(3)**(-pows) + satco/zdepth(3)/poros*slope*pows*dtt             
      q3g = q3g ** ( 1. / pows )                                                
      q3g = - ( 1. / q3g - www(3) ) * poros * zdepth(3) / dtt                   
      q3g = amax1( 0., q3g )                                                    
      q3g = amin1( q3g, www(3)*poros*zdepth(3)/dtt )                            
      if (prin) then
         print *,'roff: ',roff
         print *,'q3g: ',q3g
         print *,'temw_3: ',temw(3)
         print *,'pows: ',pows
         print *,'bee: ',bee
         print *,'satco: ',satco
         print *,'zdepth_3: ',zdepth(3)
         print *,'poros: ',poros
         print *,'slope: ',slope
      endif

      q3g = q3g + (0.002*poros*zdepth(3)*0.5/86400.*www(3)/dtt)

!----------------------------------------------------------------------         
!                                                                               
!    calculation of inter-layer exchanges of water due to gravitation           
!    and hydraulic gradient. the values of w(x) + dw(x) are used to             
!    calculate the potential gradients between layers.                          
!    modified calculation of mean conductivities follows ME-82 ), 
!    reduces recharge flux to top layer.                      
!                                                                               
!      dpdw           : estimated derivative of soil moisture potential         
!                       with respect to soil wetness. assumption of             
!                       gravitational drainage used to estimate likely          
!                       minimum wetness over the time step.                     
!                                                                               
!      qqq  (q     )  : equation (61) , SE-86                                   
!             i,i+1                                                             
!            -                                                                  
!      avk  (k     )  : equation (4.14) , ME-82                                 
!             i,i+1                                                             
!                                                                               
!----------------------------------------------------------------------         
!                                                                               
      wmax = amax1( www(1), www(2), www(3), 0.05 )                              
      wmax = amin1( wmax, 1. )                                                  
      pmax = wmax**(-bee)                                                       
      wmin = (pmax-2./( phsat*(zdepth(1)+2.*zdepth(2)+zdepth(3)))) &             
             **(-1./bee)                                                       
      wmin = amin1( www(1), www(2), www(3), wmin )                              
      wmin = amax1( wmin, 0.02 )                                                
      pmin = wmin**(-bee)                                                       
      dpdw = phsat*( pmax-pmin )/( wmax-wmin )                                  
!                                                                               
      do 8200 i = 1, 2                                                          
!                                                                               
      rsame = 0.                                                                
      avk  = temwp(i)*temwpp(i) - temwp(i+1)*temwpp(i+1)                        
      div  = temwp(i+1) - temwp(i)                                              
      if ( abs(div) .lt. 1.e-6 ) rsame = 1.                                     
      avk = satco*avk / ( ( 1. + 3./bee ) * div + rsame )                       
      avkmin = satco * amin1( temwpp(i), temwpp(i+1) )                          
      avkmax = satco * amax1( temwpp(i), temwpp(i+1) )*1.01                     
      avk = amax1( avk, avkmin )                                                
      avk = amin1( avk, avkmax )                                                
!                                                                               
!-----------------------------------------------------------------------        
!     conductivities and base flow reduced when temperature drops below         
!     freezing.                                                                 
!-----------------------------------------------------------------------        
!                                                                               
      tsnow = amin1 ( tf-0.01, tgs ) 
      areas = amin1 (0.999,13.2*snoww(2))
      tgg = tsnow*areas + tgs*(1.-areas)
      ts    = tgg*(2-i) + td*(i-1)                                              
      props = ( ts-(tf-10.) ) / 10.                                             
!     props = 1.+5*(ts-tf)                                             
      props = amax1( 0.05, amin1( 1.0, props ) )                                
      avk  = avk * props                                                        
      q3g  = q3g * props                                                        
!                                                                               
!-----------------------------------------------------------------------        
!     backward implicit calculation of flows between soil layers.               
!-----------------------------------------------------------------------        
!                                                                               
      dpdwdz = dpdw * 2./( zdepth(i) + zdepth(i+1) )                            
      aaa(i) = 1. + avk*dpdwdz*( 1./zdepth(i)+1./zdepth(i+1) ) &		 
     		 *dtt/poros				         		
      bbb(i) =-avk *   dpdwdz * 1./zdepth(2)*dtt/poros	         		
      ccc(i) = avk * ( dpdwdz * ( www(i)-www(i+1) ) + 1. +      & 		
     		(i-1)*dpdwdz*q3g*1./zdepth(3)*dtt/poros )		       
8200  continue                                                                  
!                                                                               
      denom  = ( aaa(1)*aaa(2) - bbb(1)*bbb(2) )                                
      rdenom = 0.                                                               
      if ( abs(denom) .lt. 1.e-6 ) rdenom = 1.                                  
      rdenom = ( 1.-rdenom)/( denom + rdenom )                                  
      qqq(1)   = ( aaa(2)*ccc(1) - bbb(1)*ccc(2) ) * rdenom                     
      qqq(2)   = ( aaa(1)*ccc(2) - bbb(2)*ccc(1) ) * rdenom                     
!                                                                               
!-----------------------------------------------------------------------        
!     update wetness of each soil moisture layer due to layer interflow         
!        and base flow.                                                         
!-----------------------------------------------------------------------        
!                                                                               
      www(3) = www(3) - q3g*dtt/(poros*zdepth(3))                               
      roff = roff + q3g * dtt                                                   
!                                                                               
      do 8300 i = 1, 2                                                          
!                                                                               
      qmax   =  www(i)   * (poros*zdepth(i)  /dtt)                              
      qmin   = -www(i+1) * (poros*zdepth(i+1)/dtt)                              
      qqq(i) = amin1( qqq(i),qmax)                                              
      qqq(i) = amax1( qqq(i),qmin)                                              
      www(i)   =   www(i)   - qqq(i)/(poros*zdepth(i)  /dtt)                    
      www(i+1) =   www(i+1) + qqq(i)/(poros*zdepth(i+1)/dtt)                    
8300  continue                                                                  
!
!       *** LOAD water flow & root-zone drainage PILPS DATA
      soildif=soildif+qqq(1)*dtt*1000.
      soildra=soildra+q3g*dtt*1000.
!
      do 8400 i = 1, 3                                                          
      excess = amax1(0.,(www(i) - 1.))                                          
      www(i) = www(i) - excess                                                  
      roff   = roff   + excess * poros*zdepth(i)                                
!
!       *** LOAD IN as root-drainage for PILPS
      if (i.lt.2) then
        RNOFFS= RNOFFS+ 1000.*excess*POROS*ZDEPTH(I)
      else
        RNOFFB= RNOFFB+ 1000.*excess*POROS*ZDEPTH(I)
      endif
8400  continue                                                                  
!                                                                               
!-----------------------------------------------------------------------        
!     prevent negative values of www(i)                                         
!-----------------------------------------------------------------------        
!                                                                               
      do 8402 i = 1,2                                                           
      deficit   = amax1 (0.,(1.e-12 - www(i)))                                  
      if (i.eq.1) soildif=soildif-deficit*		&
     		 zdepth(1)*poros			 
     www (i)   = www(i) + deficit			 			
     www (i+1) = www(i+1) - deficit * zdepth(i) / zdepth  (i+1) 		
 8402 continue						 			
     www(3)    = amax1(www(3),1.E-12)  		 			
! 							 
800  CONTINUE						 		
!    							 		
     IF (WWW(1) .GT.1.) then 				 
     	 WWW(2) = WWW(2) + (WWW(1)-1.) * ZDEPTH(1)/	&
     				  ZDEPTH(2)		 		
     	 soildif=soildif+(www(1)-1.)*ZDEPTH(1)		&
     			*poros*1000.			 
     	 WWW(1) = 1.					 
     end if						 
     If (WWW(2) .GT.1.) WWW(3) = WWW(3) + (WWW(2)-1.) * &
     				 ZDEPTH(2) / ZDEPTH(3) 
!							 
!      *** LOAD IN AS PILP ROOT DRAINAGE		 
     IF (WWW(2) .GT.1.) WWW(2) = 1.			 		
     if (WWW(3) .GT.1.) then 				 
     	 ROFF	= ROFF + (WWW(3)-1.)*POROS*ZDEPTH(3)	 
     	 RNOFFB=RNOFFB+((WWW(3)-1.)*ZDEPTH(3)*		&
     		POROS*1000.)
          WWW(3) = 1.                                    
      end if
!                                                                       
      RETURN                                                            
      END                                                               

!***********************************************************************
      function fun_astro(y)
!***********************************************************************
      real pi, one, two, dg2rd, yrlen, ecc, per
      real tempdmm, y

      parameter(pi = 3.1415926535898)
      parameter(dg2rd = pi/180.)
      parameter(yrlen = 365.25)
      parameter(ecc = .0167, per = 102.0*dg2rd)
      parameter(one = 1., two = 2.)

      tempdmm = (two * pi/ ((one - ecc ** 2) ** 1.5)) * (one / yrlen)  &
                * (one - ecc * cos(y - per)) ** 2

      fun_astro = tempdmm

      return
      end

