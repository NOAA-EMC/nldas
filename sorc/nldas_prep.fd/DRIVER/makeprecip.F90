!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. CONTINENTAL-SCALE
!  D                                      L  LAND MODELING AND DATA
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  ASSIMILATION PROJECT.
!  S                                      A  THIS IS THE GSFC-LDAS CODE.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  HTTP://LDAS.GSFC.NASA.GOV
!
!   GSFC - NCEP - OH - PRINCETON - WASHINGTON - RUTGERS
!
!=========================================================================
! MAKEPRECIP.F:
!
! DESCRIPTION:
!  OPENS, READS, AND INTERPOLATES STAGEIV, AND HIGGINS CPC
!  PRECIPITATION FORCING.  FORMS A MERGED PRECIPITATION PRODUCT
!  FOR USE IN LDAS DRIVER.
!
!!!!!!!  WARNING--Can't be used with Timesteps greater than 1 hour  !!!!!!!
!
! REVISION HISTORY:
!  18 FEB 2000: BRIAN COSGROVE; INITIAL CODE
!               SYNTHESIS OF DIFFERENT BASE CODES (FROM CURTIS
!               MARSHALL OF NCEP, PAUL HOUSER, JARED ENTIN) INTO
!               FUNCTIONAL UNIT.  ADDED CODE TO CHECK FOR PRE-EXISTANCE OF
!               FILES, ADDED CODE FOR ACTIVATION AT PROPER DRIVER TIMES AS
!               NEEDED.
!  11 Apr 2000: Brian Cosgrove; Changed code so that it uses the Forcing
!               land/sea mask (with inland water filled in)
!  11 May 2000: Brian Cosgrove; Changed code so it dynamically figures
!               out the length of the Merged Precip Product output
!               directory without the hardwiring
!  20 Jun 2000: Brian Cosgrove; Changed code so that sum and weightsum
!               are initialized to 0....old values of these variables
!               were persisting and causing bad values
!  14 Nov 2000: Brian Cosgrove; Changed LOGICAL *1 to LOGICAL
!   4 Jan 2001: Brian Cosgrove; Fixed way interpolation file is used
!               so that no interpolation is done if dist > 999
!   5 Sep 2001: Brian Cosgrove; Added code to use precip weighting mask
!               and to adjust convective precip according to covective
!               versus total precip ratio applied to total observed
!               precip amount
!   6 Feb 2002: Brian Cosgrove; Changed name of stageiv data read in
!  15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new
!                GRIB libraries
!   13 May 2003: Brian Cosgrove; Changed name of stageii data read in
!                from mul4 file names to ST2ml file name convention
!                Also, files are now called Stage II whereas they
!                used to be called Stage IV...same contents, just
!                a different name from NOAA
!   23 Jan 2004: Brian Cosgrove; Modified code to use PRISM adjusted
!                1/8th degree Unified Higgins precip data.
!                If user requests to use Unified data in the card
!                file, makeprecip.f will now first try to use the
!                PRISM Unified data, and then will use the normal
!                UNIFIED data if PRISM is unavailable.  PRISM and
!                non-PRISM data is kept in same UNIFIED directory.
!                Holes in PRISM are filled in with non-prism data
!                or CPC data if non-prism data is unavailable.
!   9 Nov 2006:  Brian Cosgrove; Modified code to use Schaake's PRISM
!   adjusted daily gauge product.  Holes are filled in with edas or stage2
!
!  Need something that records when EDAS precip has been used
!  to fill in a field where the Stage IV was missing.
!  (Perhaps use the FSOURCE array?)
!
!  Fix the program so that it has a fall back when there is no CPC gage data file.
!
!  1 May 2012: Youlong Xia: modifed code to use 0.5 degree CPC operational
!  Global precipitation to replace CPC US-Mexico precipitation as latter
!  stopped in March 2010  
!=========================================================================

SUBROUTINE MAKEPRECIP(LDAS,GRID)

use ldas_module     ! LDAS NON-MODEL-SPECIFIC 1-D VARIABLES
use grid_module     ! LDAS NON-MODEL-SPECIFIC GRID VARIABLES
IMPLICIT NONE

type (ldasdec) ldas
type (griddec) grid(ldas%nc,ldas%nr)
INTEGER C,R,LYR,LMO,LDA,LHR,LMN,LSS,DOY1,TS1,FILEINFO(24)
INTEGER HYR,HMO,HDA,HHR,HMN,HSS,HDOY1,HGMT1,DUMMYINT
INTEGER LLYR,LLMO,LLDA,LLHR,LLMN,LLSS,LTS1,LONGYEAR,COUNTH,Z
INTEGER IT,counti,try
REAL GMT1
REAL*8 TIME1,HTIME1

INTEGER NLDAS, NX, NY, BOUND, JPDS(25), IP, IPOPT(20), I, N,&
        LUGB, IRET, LUGI, J, JGDS(22), KF, K, KPDS(200), KM,&
        KGDSUNI(200),NUNI,IPUNI,IPUNIOPT(20),KMUNI,IBIUNI,&
        KGDSEDAS(200), KGDSS4(200), IBI, IBO, NO, DATE(25),&
        KPDSOUT(24,25), X, Y, COUNT, DIST, XX, YY, NSTAGE4,&
        KGDSLDAS(200), GDS(22), LENGDS, LDASGRD, JRET,NHX,NHY,&
        NPRISM,NHPD,ncmor


PARAMETER (NLDAS = 103936, NX = 464, NY = 224, BOUND=16,&
           NSTAGE4 = 1020800,NHX=720,NHY=360,NUNI=259200,&
           NPRISM=103936,NHPD=693,ncmor=8159252)


REAL HIGGINSUNIFIED(NX,NY),HIGGINSTEMP1D(NUNI),&
     PRISMTEMP1D(NPRISM),PRISMUNIFIED(NX,NY)
REAL NEWMASK(NPRISM)

REAL HIGGINSTEMP2D(NHX,NHY),xrat

REAL S4IN(NSTAGE4), S4OUT(NLDAS), RLAT(NLDAS),&
     RLON(NLDAS), S42DTEMP(NX,NY), S42D(NX,NY),S42DTEMPTEMP(NX,NY), &
     S42D24(24,NX,NY),&
     S4(24,NLDAS), SUM(NLDAS), WEIGHT(24,NLDAS), CPC(NLDAS),&
     HOURPREC(NLDAS), HOURS4(NLDAS), WEIGHTSUM(NLDAS),&
     CPCTMP(NX,NY), WTSUMTMP(NX,NY), WT2DTMP(24,NX,NY),&
     WT2D(24,NX,NY), S4SAVE(24,NLDAS),FINALPRECIP(24,NX,NY),&
     MASKGRADS(NX,NY),GAGE(NX,NY),HIGGINSLDAS1D(NLDAS),&
     S4SAVE2D(24,NX,NY),HPDIN(NHPD),HPDOUT(NLDAS),&
     HPD2DTEMP(NX,NY),&
     HPD2D(NX,NY),HPD2D24(24,NX,NY),HPD(24,NLDAS),&
     HPDSAVE(24,NLDAS),HPDSAVE2D(24,NX,NY),&
     CMORIN(NCMOR),CMOROUT(NLDAS),hourout(nldas),&
     CMOR2DTEMP(NX,NY),&
     CMOR2D(NX,NY),CMOR2D24(24,NX,NY),CMOR(24,NLDAS),&
     CMORSAVE(24,NLDAS),CMORSAVE2D(24,NX,NY)
     real briantemp(464,224)
     real :: lat(ldas%nc,ldas%nr),lon(ldas%nc,ldas%nr)
     real :: gmtm, tmp, tmpcz, fill, doyd, doy,tmp24(24)
     real :: zw1, zw2, czb, cze, czm, undef
     real :: angle, deg2rad
     integer rad,ifound,irad,i1,j1,ix,jx
     integer p_angle, c_angle, numarcs,badprecip(24),flag


character*1  cmorinchar(4948,1649),timestamp(4948,1649)
character*1  staid(4948,1649),testout1(4948,1649)
real realprecip(4948,1649),realprecip2(4948,1649)


LOGICAL*1 LB(NHPD), LBC(NCMOR),LO(NLDAS),LBUNI(NUNI)
LOGICAL*1 LOcmor(nldas),LOhpd(nldas),LOfinal(nldas),lofinalxy(nx,ny)
LOGICAL*1 MASKTEMP(NX,NY), FOUND, S4MASK(24,NLDAS)
LOGICAL*1 HPDMASK(24,NLDAS),CMORMASK(24,NLDAS)


LOGICAL*1 LEAP, S4AVAIL(24),LBNOSTAR(NHPD),HPDAVAIL(24),cmoravail(24)
LOGICAL*1 LBNOSTARC(NCMOR),lbs4(nstage4)

CHARACTER  FILENAM*255, CENT*2, CYR*2, CMON*2, CDAY*2, &
           CLUGB*2, CHOUR*2

INTEGER CPCFLAG,II,JJ,DUMDOY,YEAR,MONTH,DAY,HOUR,MINUTE
INTEGER SECOND,TS,ZEROI1,ZEROI2,ZEROI3,UNIFIEDFLAG,GAGEFLAG
INTEGER SYEAR,SMONTH,SDAY,SHOUR,SMINUTE,SSECOND,SDOY,STS,SLYR
INTEGER TYR,TMO,TDA,THR,TMN,TSS,TTS,TDOY,STATUS,NEEDPRECIP
INTEGER S4FLAG(24),CALLTYPE,DIRLEN,prismflag
REAL *8 TTIME
REAL TGMT
REAL *8 DUMTIME,TIMESTEP1
REAL DUMGMT,SGMT
CHARACTER   NAME1*80,NAME2*80,TNAME*80,CNAME*140
CHARACTER*1 FBASE1(80),FBASEC(100),fdir33(19)
CHARACTER*1 FDIR1(27),FDIR2(25),FDIR3(29),MERGEH2(40),FDIRC(68)
CHARACTER*100 NAMES(24),TNAMES(24),HNAME,UNIFIEDNAME,hnamemex
CHARACTER*100 PRISMNAME,NAMES2

INTEGER ICE,LI,IA,JB
integer mergeprecipx,mergeprecipy,s4precipx,s4precipy
integer ire

!        HIGGINSUNIFIED=LDAS%UDEF
!       HIGGINSTEMP1D=LDAS%UDEF
!        PRISMTEMP1D=LDAS%UDEF
!       PRISMUNIFIED=LDAS%UDEF
!        HIGGINSTEMP2D=LDAS%UDEF
!       S4IN=LDAS%UDEF
!       S4OUT=LDAS%UDEF
!       RLAT=LDAS%UDEF
!        RLON=LDAS%UDEF
!       S42DTEMP=LDAS%UDEF
!       S42D=LDAS%UDEF
!        S42D24=LDAS%UDEF
!        S4=LDAS%UDEF
!       SUM=LDAS%UDEF
!       WEIGHT=LDAS%UDEF
!       CPC=LDAS%UDEF
!        HOURPREC=LDAS%UDEF
!       HOURS4=LDAS%UDEF
!       WEIGHTSUM=LDAS%UDEF
!        CPCTMP=LDAS%UDEF
!       WTSUMTMP=LDAS%UDEF
!       WT2DTMP=LDAS%UDEF
!        WT2D=LDAS%UDEF
!       S4SAVE=LDAS%UDEF
!       FINALPRECIP=LDAS%UDEF
!        MASKGRADS=LDAS%UDEF
!       GAGE=LDAS%UDEF
!       HIGGINSLDAS1D=LDAS%UDEF
!        S4SAVE2D=LDAS%UDEF

!  SOME CONSTANTS AND SUCH
!

LUGI=0
ZEROI1=0
ZEROI2=0
ZEROI3=0
TTIME=LDAS%TIME
TDOY=LDAS%DOY
TGMT=LDAS%GMT
TYR=LDAS%YR
TMO=LDAS%MO
TDA=LDAS%DA
THR=LDAS%HR
TMN=LDAS%MN
TSS=LDAS%SS
TTS=(-60*60)
NEEDPRECIP=0        !ASSUME DON'T NEED TO MAKE NEW PRECIP
CALLTYPE=0          !INITIALIZED AT 0
                    !1 MEANS NEED TO GENERATE PRECIP USING TODAYS DATE
                    !2 MEANS NEED TO GENERATE PRECIP USING TOMORROWS DATE
DO I=1,24
   FILEINFO(i)=0    !INITIALIZED AT 0
                    !0 MEANS FILES CONTAIN STAGEIV WEIGHTED W/24hr HIGGINS UNIFIED GAGE
                    !1 MEANS FILES CONTAIN STAGEIV WEIGHTED W/24hr CPC GAGE
                    !2 MEANS EDAS DATA WAS USED TO FILL IN MISSING STAGEIV DATA
                    !AND WEIGHTED BY HIGGINS UNIFIED DATA
                    !3 MEANS EDAS DATA WAS USED TO FILL IN MISSING STAGEIV DATA
                    !AND WEIGHTED BY 24hr CPC GAGE
                    !4 MEANS NO GAGE AVAILABLE (or used), SO UNWEIGHTED STAGEIV USED
                    !5 MEANS NO GAGE AVAILABLE, AND EDAS USED TO FILL
                    !IN MISSING STAGEIV DATA, SO UNWEIGHTED EDAS DATA USED
   S4FLAG(i)=0      !INITIALIZED AT 0
                    !0 MEANS NO STAGEIV DATA IS AVAILABLE
                    !1 MEANS STAGEIV DATA WAS USED
ENDDO

!       open (unit=89,file='testc.bin',form='unformatted')
!       open (unit=88,file='test.bin',form='unformatted')

DO N=1,NLDAS
   SUM(N)=0.0

   DO II = 1, 24
      S4SAVE(II,N)=LDAS%UDEF
   END DO
ENDDO

!  Determine length of directory to make for final merged precip product
DIRLEN=0

OPEN(72,FILE='templen',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
WRITE(72,26,REC=1) LDAS%MERGEH
READ(72,22,REC=1) (MERGEH2(I),I=1,80)
CLOSE (72)
DO I=1,80
   IF(MERGEH2(I).NE.(' '))DIRLEN=DIRLEN+1
ENDDO

!  Now add seven to directory length to account for the four digit
!  year and two digit month that will be added to the base directory
!  name (and also a back slash).
DIRLEN=DIRLEN+7

!  If model is on first timestep, or model hour is 15Z
!  or 16Z, then set (reset) LDAS%PCPCOUNT to 0 so
!  that IF THEN statements are triggered below to regenerate
!  precipation files if the REGENERATE flag is set to '1'.
IF ((LDAS%TSCOUNT.EQ.1).OR.(LDAS%HR.EQ.15).OR.(LDAS%HR.EQ.16) ) LDAS%PCPCOUNT=0

!  If first timestep and current hour is less than 12,
!  or if first timestep, hour equals 12 and minutes equal 0,
!  then generate the Higgins CPC filename as well as the 24
!  Merged precip product filenames
IF ( ((LDAS%TSCOUNT.EQ.1).AND.(LDAS%HR.LT.12)).OR. &
     ((LDAS%TSCOUNT.EQ.1).AND.(LDAS%HR.EQ.12) .AND. &
     (LDAS%MN.EQ.0)) ) THEN

   CALLTYPE=1

   DO COUNT=1,24
      TYR=LDAS%YR
      TMO=LDAS%MO
      TDA=LDAS%DA
      THR=12
      TMN=0
      TSS=0
      TTS=(-3600*(COUNT-1))

     CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)

      !  If on first hour of loop, form the Higgins CPC filename
      51   FORMAT(80A1)
      52   FORMAT(80A1)
      53   FORMAT(A80)
      56   FORMAT(A80)
      58   FORMAT(A1,I4,I2,A1,I4,I2,I2,A3)
      59   FORMAT(19A1)

      IF (COUNT.EQ.1) THEN

         OPEN(90,FILE='TEMPmex',FORM='FORMATTED',ACCESS='DIRECT', &
              RECL=80, status='unknown')
         WRITE(90,58,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.ll'
         READ(90,59,REC=1) FDIR33

         DO I=1,19
            IF(FDIR33(I).EQ.(' '))FDIR33(I)='0'
         ENDDO

         WRITE(90,56,REC=1) LDAS%HIGGINS
         READ(90,52,REC=1) (FBASE1(I),I=1,80)
         C=0

         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

	 WRITE(90,51,REC=1)(FBASE1(I),I=1,C), (FDIR33(I),I=1,19)
         READ(90,53,REC=1)HNAMEmex
         CLOSE(90)

      ENDIF ! IF COUNT == 1

      !  If on first hour of loop, form the Higgins UNIFIED filename
      41   FORMAT(80A1)
      42   FORMAT(80A1)
      43   FORMAT(A80)
      46   FORMAT(A80)
      48   FORMAT(A1,I4,I2,A1,I4,I2,I2,A13)
      49   FORMAT(29A1)

      68   FORMAT(A1,I4,I2,A1,I4,I2,I2,A13)
      69   FORMAT(29A1)

      IF (COUNT.EQ.1) THEN

         OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

         WRITE(90,68,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.higg.unified'
         READ(90,69,REC=1)FDIR3

         DO I=1,29
            IF(FDIR3(I).EQ.(' '))FDIR3(I)='0'
         ENDDO

         WRITE(90,46,REC=1) LDAS%HIGGINSUNI
         READ(90,42,REC=1) (FBASE1(I),I=1,80)

         C=0
         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

         WRITE(90,41,REC=1)(FBASE1(I),I=1,C), (FDIR3(I),I=1,29)
         READ(90,43,REC=1)UNIFIEDNAME
         CLOSE(90)
      ENDIF

      !  If on first hour of loop, form the Schaake PRISM filename

      IF (COUNT.EQ.1) THEN

         OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

         WRITE(90,48,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.higgin.prism'
         READ(90,49,REC=1)FDIR3

         DO I=1,29
            IF(FDIR3(I).EQ.(' '))FDIR3(I)='0'
         ENDDO

         WRITE(90,46,REC=1) LDAS%HIGGINSUNI
         READ(90,42,REC=1) (FBASE1(I),I=1,80)

         C=0

         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

         WRITE(90,41,REC=1)(FBASE1(I),I=1,C), (FDIR3(I),I=1,29)
         READ(90,43,REC=1)PRISMNAME
         CLOSE(90)
      ENDIF

      !  Form the Merged precip product filename
      21   FORMAT(80A1)
      22   FORMAT(80A1)
      23   FORMAT(A80)
      26   FORMAT(A80)
      28   FORMAT(A1,I4,I2,A1,I2,I2,I4,A1,I2,A6)
      29   FORMAT(25A1)

      OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
      WRITE(90,28,REC=1)'/',TYR,TMO,'/',TMO,TDA,TYR,'.',THR,'.MERGE'
      READ(90,29,REC=1)FDIR2

      DO I=1,25
         IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO

      WRITE(90,26,REC=1) LDAS%MERGEH
      READ(90,22,REC=1) (FBASE1(I),I=1,80)

      C=0

      DO I=1,80
         IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      WRITE(90,21,REC=1)(FBASE1(I),I=1,C),(FDIR2(I),I=1,25)
      READ(90,23,REC=1)TNAME

      CLOSE(90)

      TNAMES(COUNT)=TNAME

      CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)

      !  Open up the newly formed Merged Precip Product filenames
      !  to see if they already exist.
      !  If they do, then don't need to generate precip.
      !  If they don't, set NEEDPRECIP flag to 1 for precip
      !  generation which occurs farther down in the program


      OPEN (UNIT=70,FILE=TNAMES(COUNT),STATUS='OLD',IOSTAT=STATUS)
      IF(STATUS.NE.0) NEEDPRECIP=1

      !  Regardless of whether or not files already exist,
      !  set NEEDPRECIP equal to '1' if the REGENERATE
      !  flag is set to '1'.  This will force
      !  it to regenerate the hourly precipitation files
      !  on the first timestep of its run even if they
      !  already exist.

      IF(LDAS%REGENERATE.EQ.1) NEEDPRECIP=1
      CLOSE(70)

   ENDDO ! ENDO COUNT=1,24

ENDIF ! ENDIF TS1 & HR <= 12 & MIN == 0

!  If first timestep and current hour is greater than 12,
!  or if first timestep, hour equals 12 and minutes greater than 0,
!  then generate the Higgins CPC filename as well as the 24
!  Merged precip product filenames
IF ( ((LDAS%TSCOUNT.EQ.1).AND.(LDAS%HR.GT.12)).OR. &
   ( (LDAS%TSCOUNT.EQ.1).AND.(LDAS%HR.EQ.12).AND. &
   (  LDAS%MN.GT.0 )) ) THEN

   CALLTYPE=2

   DO COUNT=1,24

      TYR=LDAS%YR
      TMO=LDAS%MO
      TDA=LDAS%DA
      THR=12
      TMN=0
      TSS=0

      TTS=60*60*24

      CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)
      TTS=(-3600*(COUNT-1))
      CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)

      !  If on first hour of loop, form the Higgins CPC filename
      IF (COUNT.EQ.1) THEN

         OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

         WRITE(90,58,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.ll'
         READ(90,59,REC=1)FDIR3

         DO I=1,19

            IF(FDIR3(I).EQ.(' '))FDIR3(I)='0'

         ENDDO

         WRITE(90,56,REC=1) LDAS%HIGGINS
         READ(90,52,REC=1) (FBASE1(I),I=1,80)

         C=0

         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

         WRITE(90,51,REC=1)(FBASE1(I),I=1,C), (FDIR3(I),I=1,19)
         READ(90,53,REC=1)HNAMEmex
         CLOSE(90)
      ENDIF

      !  Form the Higgins UNIFIED name
      IF (COUNT.EQ.1) THEN
         OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
         WRITE(90,68,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.higg.unified'
         READ(90,69,REC=1)FDIR3

         DO I=1,29
            IF(FDIR3(I).EQ.(' '))FDIR3(I)='0'
         ENDDO

         WRITE(90,46,REC=1) LDAS%HIGGINSUNI
         READ(90,42,REC=1) (FBASE1(I),I=1,80)
         C=0

         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

         WRITE(90,41,REC=1)(FBASE1(I),I=1,C), (FDIR3(I),I=1,29)
         READ(90,43,REC=1)UNIFIEDNAME
         CLOSE(90)

      ENDIF

      !  Form the Higgins UNIFIED PRISM name
      IF (COUNT.EQ.1) THEN

         OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
         WRITE(90,48,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.higgin.prism'
         READ(90,49,REC=1)FDIR3

         DO I=1,19
            IF(FDIR3(I).EQ.(' '))FDIR3(I)='0'
         ENDDO

         WRITE(90,46,REC=1) LDAS%HIGGINSUNI
         READ(90,42,REC=1) (FBASE1(I),I=1,80)
         C=0

         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

         WRITE(90,41,REC=1)(FBASE1(I),I=1,C), (FDIR3(I),I=1,29)
         READ(90,43,REC=1)PRISMNAME

         CLOSE(90)

      ENDIF

      !  Form the Merged precip product filename

      OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      WRITE(90,28,REC=1)'/',TYR,TMO,'/',TMO,TDA,TYR,'.',THR,'.MERGE'
      READ(90,29,REC=1)FDIR2

      DO I=1,25
         IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO

      WRITE(90,26,REC=1) LDAS%MERGEH
      READ(90,22,REC=1) (FBASE1(I),I=1,80)

      C=0

      DO I=1,80
         IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      WRITE(90,21,REC=1)(FBASE1(I),I=1,C), (FDIR2(I),I=1,25)
      READ(90,23,REC=1)TNAME

      CLOSE(90)
      TNAMES(COUNT)=TNAME

      CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)

      !  Open up the newly formed Merged Precip Product filenames
      !  to see if they already exist.
      !  If they do, then don't need to generate precip.
      !  If they don't, set NEEDPRECIP flag to 1 for precip
      !  generation which occurs farther down in the program


      OPEN (UNIT=70,FILE=TNAMES(COUNT),STATUS='OLD',IOSTAT=STATUS)

      IF(STATUS.NE.0) NEEDPRECIP=1

      !  Regardless of whether or not files already exist,
      !  set NEEDPRECIP equal to '1' if the REGENERATE flag
      !  is set to '1'.  This will force
      !  it to regenerate the hourly precipitation files
      !  on the first timestep of its run even if they
      !  already exist.
      IF(LDAS%REGENERATE.EQ.1) NEEDPRECIP=1

      CLOSE(70)

   ENDDO ! COUNT = 1, 24

ENDIF ! TS == 1 & HR >= 12 & MIN = 0

!  If NOT first timestep and current hour equals 12,
!  and minutes greater than 0, or if NOT first timestep
!  and current hour equals 13, then generate the Higgins
!  CPC filename as well as the 24 Merged precip filenames

IF ( ((LDAS%TSCOUNT.NE.1).AND.(LDAS%HR.EQ.12).AND. &
   ( LDAS%MN.GT.0)).OR.((LDAS%TSCOUNT.NE.1).AND. &
   ( LDAS%HR.EQ.13)) ) THEN

   CALLTYPE=2

   DO COUNT=1,24

      TYR=LDAS%YR
      TMO=LDAS%MO
      TDA=LDAS%DA
      THR=12
      TMN=0
      TSS=0

      TTS=60*60*24
      CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)
      TTS=(-3600*(COUNT-1))
      CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)

      ! If on first hour of loop, form the Higgins CPC filename

      IF (COUNT.EQ.1) THEN

         OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

         WRITE(90,58,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.ll'
         READ(90,59,REC=1)FDIR3

         DO I=1,19
            IF(FDIR3(I).EQ.(' '))FDIR3(I)='0'
         ENDDO

         WRITE(90,56,REC=1) LDAS%HIGGINS
         READ(90,52,REC=1) (FBASE1(I),I=1,80)
         C=0

         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

         WRITE(90,51,REC=1)(FBASE1(I),I=1,C), (FDIR3(I),I=1,19)
         READ(90,53,REC=1)HNAMEmex
         CLOSE(90)

      ENDIF

      !  Form the HIGGINS UNIFIED name
      IF (COUNT.EQ.1) THEN

         OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
         WRITE(90,68,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.higg.unified'
         READ(90,69,REC=1)FDIR3

         DO I=1,29
            IF(FDIR3(I).EQ.(' '))FDIR3(I)='0'
         ENDDO

         WRITE(90,46,REC=1) LDAS%HIGGINSUNI
         READ(90,42,REC=1) (FBASE1(I),I=1,80)

         C=0

         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

         WRITE(90,41,REC=1)(FBASE1(I),I=1,C), (FDIR3(I),I=1,29)
         READ(90,43,REC=1)UNIFIEDNAME

         CLOSE(90)

      ENDIF

      !  Form the HIGGINS UNIFIED PRISM name
      IF (COUNT.EQ.1) THEN

         OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

         WRITE(90,48,REC=1)'/',TYR,TMO,'/',TYR,TMO,TDA,'.higgin.prism'

         READ(90,49,REC=1)FDIR3

         DO I=1,29
            IF(FDIR3(I).EQ.(' '))FDIR3(I)='0'
         ENDDO

         WRITE(90,46,REC=1) LDAS%HIGGINSUNI
         READ(90,42,REC=1) (FBASE1(I),I=1,80)

         C=0
         DO I=1,80
            IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
         ENDDO

         WRITE(90,41,REC=1)(FBASE1(I),I=1,C), (FDIR3(I),I=1,29)
         READ(90,43,REC=1)PRISMNAME
         CLOSE(90)

      ENDIF

      !  Form the Merged precip product filename

      OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
      WRITE(90,28,REC=1)'/',TYR,TMO,'/',TMO,TDA,TYR,'.',THR,'.MERGE'
      READ(90,29,REC=1) FDIR2

      DO I=1,25
         IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO

      WRITE(90,26,REC=1) LDAS%MERGEH
      READ(90,22,REC=1) (FBASE1(I),I=1,80)

      C=0
      DO I=1,80
         IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      WRITE(90,21,REC=1)(FBASE1(I),I=1,C), (FDIR2(I),I=1,25)
      READ(90,23,REC=1)TNAME

      CLOSE(90)
      TNAMES(COUNT)=TNAME

      CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)

      !  Open up the newly formed Merged Precip Product filenames
      !  to see if they already exist.
      !  If they do, then don't need to generate precip.
      !  If they don't, set NEEDPRECIP flag to 1 for precip
      !  generation which occurs farther down in the program


      OPEN (UNIT=70,FILE=TNAMES(COUNT),STATUS='OLD',IOSTAT=STATUS)

      IF(STATUS.NE.0) NEEDPRECIP=1

      !  Regardless of whether or not files already exist,
      !  set NEEDPRECIP equal to '1' if the REGENERATE flag is
      !  set to 1 and if hourly precipitation
      !  files have not been regenerated already...ie, LDAS%PCPCOUNT
      !  would equal 1 if the files had already been regenerated.
      !  This will force the model to regenerate the hourly precipitation
      !  files even if they already exist.
      !  Set LDAS%PCPCOUNT to '1' after NEEDPRECIP is set to '1' so that
      !  this regeneration project will not be triggered the next time
      !  around.  LDAS%PCPCOUNT is reset to 0 during the hours of 15Z
      !  and 16Z in preparation for the next time this current loop
      !  is activated and precipitation must once again be regenerated.

      IF((LDAS%REGENERATE.EQ.1).AND.(LDAS%PCPCOUNT.EQ.0)) NEEDPRECIP=1
      LDAS%PCPCOUNT=1

      CLOSE(70)

   ENDDO ! COUNT = 1, 24

ENDIF ! TS != 1 & HR = 12 OR HR = 13 & MIN = 0


!  If NEDPRECIP flag has been set to 1, then generate the
!  required merged precip product by reading in STAGEIV precip
!  and 24 Gage precip, deriving hourly weights from the STAGEIV
!  precip, and applying these weights to the Gage precip.
IF (NEEDPRECIP.EQ.1) THEN

   JPDS = -1
   JPDS(5) = 61
   JGDS = -1
   IP = 3

   DO I=1,20
     IPOPT(I)=0
   ENDDO

   IPOPT(1) = 4
   IPOPT(2) = -1
   LDASGRD = 110
   IBI = 1
   KM = 1
   TS=60*60
   STS=LDAS%TS
   DUMTIME=LDAS%TIME
   DUMDOY=LDAS%DOY
   DUMGMT=LDAS%GMT
   I=1
   J=0

   SYEAR=LDAS%SYR
   SMONTH=LDAS%SMO
   SDAY=LDAS%SDA
   SHOUR=LDAS%SHR
   SMINUTE=LDAS%SMN
   SSECOND=LDAS%SSS
   SDOY=LDAS%SDOY


   DO II=1,24

      LYR=LDAS%YR
      HYR=LDAS%YR
      LMO=LDAS%MO
      HMO=LDAS%MO
      LDA=LDAS%DA
      HDA=LDAS%DA
      LHR=13
      HHR=12
      LMN=0
      HMN=0
      LSS=0
      HSS=0

      IF (CALLTYPE.EQ.1) THEN

         TTS=-60*60*24
         CALL TICK(TIME1,DOY1,GMT1,LYR,LMO,LDA,LHR,LMN,LSS,TTS)
         CALL TICK(HTIME1,HDOY1,HGMT1,HYR,HMO,HDA,HHR,HMN,HSS,TTS)

      ENDIF

      !  OPEN (UNIT=88,FILE='MERGE.OUT',FORM='UNFORMATTED')

      TS1=(3600*(II-1))
!      TS1=(3600*(II))

      CALL TICK(TIME1,DOY1,GMT1,LYR,LMO,LDA,LHR,LMN,LSS,TS1)
      CALL TICK(HTIME1,HDOY1,HGMT1,HYR,HMO,HDA,HHR,HMN,HSS,TS1)

      IF (MOD((LHR*3),2).EQ.0) THEN
         try=34
      ELSE
         try=35
      ENDIF

      CALL GETETAMOD(LDAS,GRID,LYR,LMO,LDA,LHR,LMN,try)

      DO R=1,LDAS%NR
         DO C=1,LDAS%NC

            IF ((GRID(C,R)%FMASK.GE.1.0).AND.(GRID(C,R)%PRECIP(1).NE.LDAS%UDEF)) THEN

               IF (GRID(C,R)%PRECIP(1).LT.0.0) THEN

                  write (*,'(A24,F11.8,A6,I3,A1,I3,A21,I2)') &
                         'MERGE ETA OUT OF RANGE (', &
                         GRID(C,R)%PRECIP(1),') at (',c,',',r, &
                         ') CORRECTING TO 0 at ',II
                  GRID(C,R)%PRECIP(1)=0.0

               ENDIF

               IF (GRID(C,R)%PRECIP(1).GT.0.03888) THEN

                  write (*,'(A24,F11.8,A6,I3,A1,I3,A34,I2)') &
                         'MERGE ETA OUT OF RANGE (', &
                         GRID(C,R)%PRECIP(1),') at (',c,',',r, &
                         ') CORRECTING TO 0.03888 mm/sec at ',II
                  GRID(C,R)%PRECIP(1)=0.03888

                ENDIF

            ENDIF

            GRID(C,R)%ETAPRECIP(II)=GRID(C,R)%PRECIP(1)

         ENDDO ! C=1, NC

      ENDDO ! R=1, NR

      LDAS%STAGEIVTIME=TIME1

      !  GENERATE FILENAME FOR STAGE 2 PRECIPITATION

      91   FORMAT(80A1)
      92   FORMAT(80A1)
      93   FORMAT(A80)
      96   FORMAT(A80)
      ! 98   FORMAT(A1,I4,I2,A6,I2,I2,I2,I2,A4)
      ! 98   FORMAT(A1,I4,I2,A6,I4,I2,I2,A1,I2)

      90   FORMAT(A1,I4,I2,A6,I4,I2,I2,I2,A4)
      99   FORMAT(27A1)

      OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80,status='unknown')

      IF(LYR.LT.2000) SLYR=LYR-1900
      IF(LYR.GE.2000) SLYR=LYR-2000

      !  WRITE(90,90,REC=1)'/',LYR,LMO,'/mul4.',LYR,LMO,LDA,'.',LHR
      WRITE(90,90,REC=1)'/',LYR,LMO,'/ST2ml',LYR,LMO,LDA,LHR,'.Grb'
      !  WRITE(90,90,REC=1)'/',LYR,LMO,'/USAml',LMO,LDA,SLYR,
      !  &  LHR,'.Grb'

      READ(90,99,REC=1)FDIR1

      DO I=1,27
         IF(FDIR1(I).EQ.(' '))FDIR1(I)='0'
      ENDDO

      WRITE(90,96,REC=1) LDAS%STAGEIV
      READ(90,92,REC=1) (FBASE1(I),I=1,80)

      C=0

      DO I=1,80
         IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      WRITE(90,91,REC=1)(FBASE1(I),I=1,C), (FDIR1(I),I=1,27)
      READ(90,93,REC=1)NAMES2

      CLOSE(90)
      ! GENERATE FILENAME FOR HPD PRECIPITATION
      ! 98   FORMAT(A1,I4,I2,A6,I2,I2,I2,I2,A4)
      ! 98   FORMAT(A1,I4,I2,A6,I4,I2,I2,A1,I2)

      98   FORMAT(A1,I4,I2,A5,I4,I2,I2,A1,I2,A4)

      OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      IF(LYR.LT.2000) SLYR=HYR-1900
      IF(LYR.GE.2000) SLYR=HYR-2000

      !  WRITE(90,98,REC=1)'/',LYR,LMO,'/mul4.',LYR,LMO,LDA,'.',LHR
      !  WRITE(90,98,REC=1)'/',LYR,LMO,'/ST2ml',LYR,LMO,LDA,LHR,'.Grb'
      WRITE(90,98,REC=1)'/',HYR,HMO,'/HPD.',HYR,HMO,HDA,'.',HHR,'.bin'
      !  WRITE(90,98,REC=1)'/',LYR,LMO,'/USAml',LMO,LDA,SLYR,
      !  &  LHR,'.Grb'

      READ(90,99,REC=1)FDIR1

      DO I=1,27
         IF(FDIR1(I).EQ.(' '))FDIR1(I)='0'
      ENDDO

      WRITE(90,96,REC=1) LDAS%HPD
      READ(90,92,REC=1) (FBASE1(I),I=1,80)

      C=0
      DO I=1,80
         IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      WRITE(90,91,REC=1)(FBASE1(I),I=1,C), (FDIR1(I),I=1,27)
      READ(90,93,REC=1)NAME1

      CLOSE(90)

      !  GENERATE FILENAME FOR CMORPH PRECIPITATION
      !  advt-8km-intrp-prim-sat-spat-2lag-2.5+5dovlp8kmIR-2002120410
      31   FORMAT(140A1)
      32   FORMAT(80A1)
      33   FORMAT(A140)
      36   FORMAT(A80)
      ! 28   FORMAT(A1,I4,I2,A6,I2,I2,I2,I2,A4)
      ! 28   FORMAT(A1,I4,I2,A6,I4,I2,I2,A1,I2)
      38   FORMAT(A1,I4,I2,A51,I4,I2,I2,I2)
      ! 28   FORMAT(A1,I4,I2,A5,I4,I2,I2,A1,I2,A4)
      39   FORMAT(68A1)
      OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=140)

      IF(LYR.LT.2000) SLYR=HYR-1900
      IF(LYR.GE.2000) SLYR=HYR-2000

      !  WRITE(90,98,REC=1)'/',LYR,LMO,'/mul4.',LYR,LMO,LDA,'.',LHR
      !  WRITE(90,98,REC=1)'/',LYR,LMO,'/ST2ml',LYR,LMO,LDA,LHR,'.Grb'
      !  WRITE(90,98,REC=1)'/',HYR,HMO,'/HPD.',HYR,HMO,HDA,'.',HHR,'.bin'
      WRITE(90,38,REC=1) '/',HYR,HMO, &
                         '/advt-8km-intrp-prim-sat-spat-2lag-2.5+5dovlp8kmIR-', &
                         HYR,HMO,HDA,HHR

      !  WRITE(90,98,REC=1)'/',LYR,LMO,'/USAml',LMO,LDA,SLYR,
      !  &  LHR,'.Grb'

      READ(90,39,REC=1)FDIRC

      DO I=1,68
         IF(FDIRC(I).EQ.(' '))FDIRC(I)='0'
      ENDDO

      WRITE(90,36,REC=1) LDAS%CMORPH
      READ(90,32,REC=1) (FBASE1(I),I=1,80)

      C=0
      DO I=1,80
         IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      WRITE(90,31,REC=1)(FBASE1(I),I=1,C), (FDIRC(I),I=1,68)
      READ(90,33,REC=1)CNAME

      CLOSE(90)


      ! GENERATE FILENAME FOR MERGED PRECIPITATION
      81   FORMAT(80A1)
      82   FORMAT(40A1)
      83   FORMAT(A80)
      86   FORMAT(A40)
      88   FORMAT(A1,I4,I2,A1,I2,I2,I4,A1,I2,A6)
      89   FORMAT(25A1)

      OPEN(90,FILE='MPTEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      WRITE(90,88,REC=1)'/',LYR,LMO,'/',LMO,LDA,LYR,'.',LHR,'.MERGE'
      READ(90,89,REC=1)FDIR2

      DO I=1,25
         IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO

      WRITE(90,86,REC=1) LDAS%MERGEH
      READ(90,82,REC=1) (FBASE1(I),I=1,40)

      C=0
      DO I=1,40
         IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      WRITE(90,81,REC=1)(FBASE1(I),I=1,C), (FDIR2(I),I=1,25)
      READ(90,83,REC=1)NAME2

      CLOSE(90)
      NAMES(II)=NAME2

      !  MAKE A GDS FOR THE OUTPUT GRID (GDS=Grid Definition Section)
      CALL MAKGDS (LDASGRD,KGDSLDAS,GDS,LENGDS,IRET)

! initialize lo bitmaps and badprecip array
	do i=1,nldas
	locmor(i)=.false.
	lohpd(i)=.false.
        lo(i)=.false.
        lofinal(i)=.false.
	enddo
	do j=1,ny
	do i=1,nx
	lofinalxy(i,j)=.false.
	enddo
	enddo
	do i=1,24
	badprecip(i)=0
	enddo
	flag=0

      IF (((ldas%yr.eq.2002).and.(ldas%mo.lt.12)).or.(ldas%yr.lt.2002)) THEN

         !  READ IN THE HPD GRID AND INTERPOLATE TO LDAS GRID.
         !  ELSE, PROCEED AND FILL HOLES IN
         !  HPD GRID WITH EDAS VALUES BELOW.

         !  Read in hpd files
         HPDIN=9999
         OPEN(unit=38,file=name1,form='unformatted',status='old',ERR=97)

         READ(38) HPDIN
 97    continue

        HPDAVAIL(II)=.FALSE.
!	print *,'read in hpd and set to false',HPDAVAIL(II)
! set HPDAVAIL to true if find defined data
         DO i=1,NHPD
            if (HPDIN(i).lt.9000)  HPDAVAIL(II)=.TRUE.
            !  HPDIN(i)=HPDIN(i)*1000
         ENDDO
!	print *,'hpd now is ',HPDAVAIL(II)



         DO I=1,200
            KGDSS4(I)=0
         ENDDO

         ! xdef   33 linear   -140.000  2.500
         ! ydef    21 linear  20.000  2.000

         KGDSS4(1)=0
         KGDSS4(2)=33
         KGDSS4(3)=21
         KGDSS4(4)=20000
         KGDSS4(5)=-140000
         KGDSS4(6)=128
         KGDSS4(7)=60000
         KGDSS4(8)=-60000
         KGDSS4(9)=2500
         KGDSS4(10)=2000
         KGDSS4(11)=64
         KGDSS4(12)=0
         KGDSS4(13)=0
         KGDSS4(14)=0
         KGDSS4(15)=0
         KGDSS4(16)=0
         KGDSS4(17)=0
         KGDSS4(18)=0
         KGDSS4(19)=0
         KGDSS4(20)=255
         KGDSS4(21)=0
         KGDSS4(22)=0
         IP =  3

         DO I=1,20
            IPOPT(I)=-1
         ENDDO
	IPOPT(1)=4

         KMUNI = 1
         IBIUNI = 1

         LB=.TRUE.

         CALL IPOLATES (IP,IPOPT,KGDSS4,KGDSLDAS,NHPD,NLDAS,KM,IBI, &
                        LB,HPDIN,NO,RLAT,RLON,IBO,LOhpd,hourOUT,IRET)

         DO i=1,nldas
            hourout(i)=hourout(i)*25.4
         ENDDO

         !print *,'hpd iret=',iret
         !open (unit=88,file='out.bin',form='unformatted')
         !write(88) hourout
         !close(88)
         !stop

         IF (.not.hpdavail(ii)) THEN
!	print *,'setting hourout to eta'
            COUNT = 0
            DO Y = 1, NY
               DO X = 1, NX
!                  S4(II,X+COUNT)=GRID(X,Y)%ETAPRECIP(II)
                  hourout(X+COUNT)=GRID(X,Y)%ETAPRECIP(II)
               END DO
               COUNT = COUNT + NX
            END DO
         ENDIF

            COUNT = 0
            DO Y = 1, NY
               DO X = 1, NX
                  if (hourout(x+count).lt.0) then
                  lohpd(x+count)=.false.
                  endif
               END DO
               COUNT = COUNT + NX
            END DO



      ELSE

         !  READ IN THE CMORPH GRID AND INTERPOLATE TO LDAS GRID.
         !  ELSE, PROCEED AND FILL HOLES IN
         !  CMORPH GRID WITH EDAS VALUES BELOW.

         ! read in cmorph files
         cmorAVAIL(ii)=.FALSE.
         open(unit=38,file=cname, status='old',access='direct', &
              form='unformatted',recl=4948*1649*3,err=133)

         read (38,rec=1) cmorinchar, timestamp, staid
         ! flipping data from 60N->60S to 60S->60N LIS standard

         r=1649

         do i = 1,1649
            do j = 1,4948
               testout1(j,r) = cmorinchar(j,i)
            enddo
            r = r - 1
         enddo

         do i = 1,1649
            do j = 1,4948
               if(ichar(testout1(j,i)).eq.255) then
                  realprecip(j,i) = -1
               else
                  realprecip(j,i) = ichar(testout1(j,i))*0.2
               endif
            enddo
         enddo

         counti=1

         do j = 1,1649
            do i=1,4948
               cmorin(counti)=realprecip(i,j)
               counti=counti+1
            enddo
         enddo

         ! read in cmorph files for second half of hour
         read (38,rec=2) cmorinchar, timestamp, staid

         !
         ! flipping data from 60N->60S to 60S->60N LIS standard
         !
         r=1649
         do i = 1,1649
            do j = 1,4948
               testout1(j,r) = cmorinchar(j,i)
            enddo
            r = r - 1
         enddo

         do i = 1,1649
            do j = 1,4948
               if(ichar(testout1(j,i)).eq.255)then
                  realprecip2(j,i) = -1
               else
                  realprecip2(j,i) = ichar(testout1(j,i))*0.2
               endif
            enddo
         enddo

         counti=1
         do j = 1,1649
            do i=1,4948
               ! add up two cmorph fields to create hourly total
               cmorin(counti)=realprecip2(i,j)+realprecip(i,j)
               counti=counti+1
            enddo
         enddo

         !  open (unit=44,file='cmorph.bin',status='unknown',
         !  &  form='unformatted')
         !  write(44)cmorin
         !  close (44)
         !  stop

         cmorAVAIL(ii)=.TRUE.

         DO I=1,200
            KGDSS4(I)=0
         ENDDO

         KGDSS4(1)=0
         KGDSS4(2)=4948
         KGDSS4(3)=1649
         KGDSS4(4)=-59963
!         KGDSS4(5)=-179963
         KGDSS4(5)=37
         KGDSS4(6)=128
         KGDSS4(7)=59963
!         KGDSS4(8)=179963
         KGDSS4(8)=359963
         KGDSS4(9)=72
         KGDSS4(10)=72
         KGDSS4(11)=64
         KGDSS4(12)=0
         KGDSS4(13)=0
         KGDSS4(14)=0
         KGDSS4(15)=0
         KGDSS4(16)=0
         KGDSS4(17)=0
         KGDSS4(18)=0
         KGDSS4(19)=0
         KGDSS4(20)=255
         KGDSS4(21)=0
         KGDSS4(22)=0
         IP =  3

         DO I=1,20
            IPOPT(I)=-1
         ENDDO
        IPOPT(1)=4


         KM = 1
         IBI = 1
         LBC=.TRUE.

         CALL IPOLATES (IP,IPOPT,KGDSS4,KGDSLDAS,NCMOR,NLDAS,KM,IBI, &
                        LBC,CMORIN,NO,RLAT,RLON,IBO,LOcmor,hourOUT,IRET)

         !  write(89) cmorin
         !  write(88) s4out

         133  CONTINUE

         IF (.not.cmoravail(ii)) THEN

            COUNT = 0
            DO Y = 1, NY
               DO X = 1, NX
!                  S4(II,X+COUNT)=GRID(X,Y)%ETAPRECIP(II)
                  hourout(X+COUNT)=GRID(X,Y)%ETAPRECIP(II)
               END DO
               COUNT = COUNT + NX
            END DO
         ENDIF


            COUNT = 0
            DO Y = 1, NY
               DO X = 1, NX
	          if (hourout(x+count).lt.0) then
                  locmor(x+count)=.false.
	          endif
               END DO
               COUNT = COUNT + NX
            END DO




         ! endif the date check for hpd versus cmorph data usage

!	if (ii.eq.12) then
!               open (unit=89,file='out.bin',form='unformatted', &
!                status='unknown')
!               write(89) hourout
!               close(89)
!               stop
!	endif
      ENDIF

      ! now read in stage2 radar data over CONUS if 199607 or later
      do n=1,nldas
         lo(n)=.false.
      enddo

      IF (((ldas%yr.eq.1996).and.(ldas%mo.ge.7)).or.(ldas%yr.gt.1996)) THEN
         !  READ IN THE STAGE IV GRID AND INTERPOLATE TO LDAS GRID.  IF ENTIRE
         !  STAGE IV GRID IS MISSING (I.E., THE HOURLY FILE IS MISSING) THEN
         !  REPLACE GRID WITH EDAS GRID.  ELSE, PROCEED AND FILL HOLES IN
         !  STAGE IV GRID WITH EDAS VALUES BELOW.

         IF (MOD((LHR*3),2).EQ.0) THEN
            LUGB=34
         ELSE
            LUGB=35
         ENDIF

         J=0
         WRITE (CLUGB, '(I2.2)') LUGB
         CALL BAOPENR (LUGB, NAMEs2, IRET)
         CALL GETGB (LUGB,LUGI,NSTAGE4,J,JPDS,JGDS,KF,K,KPDS,KGDSS4, &
                     LBs4,S4IN,IRET)

         ! check for files that open successfully but which are all zero
         IF (iret.eq.0) THEN

            iret=1
            DO y=1,nstage4
               IF (s4in(y).gt.0) iret=0
            ENDDO

         ENDIF


         IF (IRET.NE.0) THEN
            S4AVAIL(II) = .FALSE.
            COUNT =0

            DO Y = 1, NY
               DO X = 1, NX
                  ! S4(II,X+COUNT)=GRID(X,Y)%ETAPRECIP(II)
                  ! S4(II,X+COUNT)=hourout(x+count)
                  S42D(X,Y)=hourout(x+count)/3600.0
               END DO
               count=count+nx
            END DO

            !DO i = 1,nldas
            !   S4out(i)=hourout(i)
            !END DO

         ELSE

            S4AVAIL(II) = .TRUE.

         DO I=1,20
            IPOPT(I)=-1
         ENDDO
        IPOPT(1)=4



!	if (ii.eq.12)  then
!	open (unit=55,file='s4in.bin',form='unformatted',status='unknown')
!	write(55) s4in
!	close(55)	
!	ENDIF
            CALL IPOLATES (IP,IPOPT,KGDSS4,KGDSLDAS,NSTAGE4,NLDAS,KM,IBI, &
                           LBs4,S4IN,NO,RLAT,RLON,IBO,LO,S4OUT,IRET)

!        if (ii.eq.11)  then
!        open (unit=55,file='s4out.bin',form='unformatted',status='unknown')
!        write(55) s4out
!        close(55)
!	do y=1,nldas
!	print *,'y,s4out',y,s4out(y)
!	enddo
!        ENDIF

         END IF
         CALL BACLOSE (LUGB,JRET)

         !  SAVE S4OUT AND MASK FOR LATER OUTPUT IF STAGE IV AVAILABLE.


         DO N = 1, NLDAS

            IF ((LO(N)).AND.(S4OUT(N).NE.LDAS%UDEF)) THEN

               IF (S4OUT(N).LT.0.0) THEN
                  write (*,'(A17,F11.8,A6,I3,A1,I3,A21,I2)') &
                           'S4 OUT OF RANGE (', &
                            S4OUT(N),') at (',c,',',r, &
                            ') CORRECTING TO 0 at ',II
                  S4OUT(N)=0.0
               ENDIF

               IF (S4OUT(N).GT.140.0) THEN

                  write (*,'(A17,F11.8,A6,I3,A1,I3,A25,I2)') &
                         'S4 OUT OF RANGE (', &
                          S4OUT(N),') at (',c,',',r, &
                          ') CORRECTING TO 140.0 at ',II
                  S4OUT(N)=140.0
               ENDIF
            ENDIF

            IF (S4AVAIL(II)) THEN
               S4SAVE(II,N)=S4OUT(N)
            ENDIF

         END DO

         COUNT = 0

         DO Y = 1, NY

            DO X = 1, NX
               S4SAVE2D(II,X,Y) = (S4SAVE(II,X+COUNT))

            END DO
            COUNT = COUNT + NX

         END DO

      ENDIF ! endif the date based stage 2 retreival
	

      DO N = 1, NLDAS
         IF ((.not.s4avail(ii)).or.(.not.lo(n))) s4out(n)=hourout(n)
!         IF ((.not.s4avail(ii)).or.(.not.lo(n))).and.   &
!         ((lohpd(x+count)).or.(locmor(x+count))) s4out(n)=hourout(n)

      ENDDO

      IF (S4AVAIL(II)) THEN

         DO N = 1, NLDAS
            S4MASK(II,N)=LO(N)
         END DO

         !  PUT 1D ARRAYS OF DATA ONTO 2D LDAS GRID (FOR NEIGHBOR SEARCHING
         !  ALGORITHM BELOW).  AND CHANGE PRECIP FROM MM INTO RATE MM/SEC

         COUNT = 0
         DO Y = 1, NY
            DO X = 1, NX
               if ((lo(x+count)).or.(lohpd(x+count)).or.(locmor(x+count))) then
               lofinalxy(x,y)=.true.
	       endif

               S42DTEMP(X,Y) = (S4OUT(X+COUNT))/3600
               MASKGRADS(X,Y)=0.0
               IF(LOfinalxy(x,y)) MASKGRADS(X,Y)=1.0
               MASKTEMP(X,Y) = LOfinalxy(x,y)

            END DO

            COUNT = COUNT + NX

         END DO

         !  REPLACE MISSING DATA POINTS WITH NEAREST NEIGHBOR IF NEAREST NEIGHBOR
         !  WITHIN ABOUT 240 KM ("BOUND").  ELSE, REPLACE WITH EDAS VALUE.
  rad=20
  deg2rad = 0.01745329
  DO y=1, ldas%nr
     DO x=1, ldas%nc
! out2 new and out is orig array
!        out2(i,j) = out(i,j)
!        IF (out(i, j) .EQ. undef) THEN
        IF ( .NOT. MASKTEMP(X,Y)) THEN
           ifound = 0
! search radius is rad
           DO irad = 1, rad
             numarcs = irad * 8
             p_angle = 360000 / numarcs
             DO c_angle=0,360000-p_angle,p_angle
               angle = c_angle / 1000.00
               j1 = y + INT(COS(angle*deg2rad)*irad)
               i1 = x + INT(SIN(angle*deg2rad)*irad)
               jx = j1
               ix = i1
               IF (jx .LE. 0 ) jx = 1
               IF (jx .GE. ldas%nr ) jx = ldas%nr
               IF (ix .LE. 0 ) ix = 1
               IF (ix .GT. ldas%nc ) ix = ldas%nc
               !write(*,'(f8.2,5i4)') angle,i1,j1,ix,jx,irad
!               IF(out(ix, jx) .NE. undef ) then
               IF(MASKTEMP(ix,jx) ) then

                 ifound = ifound + 1
                 tmp=S42DTEMP(ix,jx)
                 S4FLAG(ii)=1
               ENDIF
! number of points to find
               IF(ifound .GE. 1) EXIT
             ENDDO
             IF (ifound .GE. 1) EXIT
           ENDDO
           IF (ifound .GE. 1) THEN
             S42D(X,Y) = tmp
           ELSE
              S42D(X,Y)=GRID(X,Y)%ETAPRECIP(II)
           ENDIF
        ELSE
        S42D(X,Y)=S42DTEMP(X,Y)
        S4FLAG(ii)=1
        ENDIF !IF MISSING DATA AT POINT

     ENDDO !NC
  ENDDO !NR

!	print *,'finished with replacement'
!         DO Y = 1, NY
!            DO X = 1, NX
!               FOUND = .TRUE.
!               IF (.NOT. MASKTEMP(X,Y)) THEN
!                  FOUND = .FALSE.
!                  IF ((X.LE.(NX-BOUND)).AND.(X.GT.BOUND).AND. &
!                     (Y.LE.(NY-BOUND)).AND.(Y.GT.BOUND)) THEN
!
!                     DIST = 1
!
!                     DO WHILE ((DIST .LE. BOUND).AND.(.NOT. FOUND))
!                        YY = -1*DIST
!                        DO WHILE ((YY .LE. DIST).AND.(.NOT. FOUND))
!                           XX = -1*DIST
!                           DO WHILE ((XX .LE. DIST).AND.(.NOT. FOUND))
!                              IF (MASKTEMP(X+XX,Y+YY)) THEN
!                                 S42D(X,Y) = S42DTEMP(X+XX,Y+YY)
!                                 S4FLAG(ii)=1
!                                 FOUND = .TRUE.
!                              END IF
!                              XX = XX+1
!                           END DO
!                           YY = YY+1
!                        END DO
!                        DIST = DIST+1
!                     END DO
!                  ELSE
!
!                     S42D(X,Y)=GRID(X,Y)%ETAPRECIP(II)
!                     FOUND = .TRUE.
!
!                  END IF
!
!               ELSE
!                  S42D(X,Y)=S42DTEMP(X,Y)
!                  S4FLAG(ii)=1
!               END IF
!
!               IF (.NOT. FOUND) THEN
!                  S42D(X,Y)=GRID(X,Y)%ETAPRECIP(II)
!               END IF

               !! set hourly precip over mexico to NARR Precip for pre-cmorph period
               !       IF (((ldas%yr.eq.2002).and.(ldas%mo.lt.12))
               !     &  .or.(ldas%yr.lt.2002)) then
               !        IF ((grid(x,y)%precipweight.lt.16).and.(y.lt.112)) THEN
               !                S42D(X,Y)=(grid(x,y)%precipweight/16.0)*S42D(X,Y)
               !     &          +((1.0-(grid(x,y)%precipweight/16.0))
               !     &          *GRID(X,Y)%ETAPRECIP(II))
               !        ENDIF
               !        ENDIF
               !
               !! set hourly precip over canada to NARR Precip for whole period
               !        IF((grid(x,y)%precipweight.lt.16).and.(y.ge.112))
               !     &      S42D(X,Y)=(grid(x,y)%precipweight/16.0)*S42D(X,Y)
               !     &      + ((1.0-(grid(x,y)%precipweight/16.0))
               !     &      *GRID(X,Y)%ETAPRECIP(II))

!            END DO
!         END DO

         ! do j=1,224
         ! do i=1,464
         !      if(s42d(i,j).gt.0) print *,'s42d=',s42d(i,j),i,j
         !        enddo
         !        enddo
         !        open (unit=45,file='out.bin',form='unformatted',
         !     &  status='unknown')
         !        write(45) s42d
         !        close(45)
         !        stop


! else stageiv is not available
      ELSE
	count=0

         !  REPLACE MISSING DATA POINTS WITH NEAREST NEIGHBOR IF NEAREST NEIGHBOR
         !  WITHIN ABOUT 240 KM ("BOUND").  ELSE, REPLACE WITH EDAS VALUE.
  rad=20
  deg2rad = 0.01745329
        DO y=1, ldas%nr
        DO x=1, ldas%nc
               S42DTEMPTEMP(X,Y) = S4OUT(X+COUNT)
               if ((lo(x+count)).or.(lohpd(x+count)).or.(locmor(x+count)))  &
                lofinalxy(x,y)=.true.
	enddo	
	count=count+nx
	enddo

	count=0
  DO y=1, ldas%nr
     DO x=1, ldas%nc
! out2 new and out is orig array
!        out2(i,j) = out(i,j)
!        IF (out(i, j) .EQ. undef) THEN
        IF ( .NOT. lofinalxy(x,y)) THEN
           ifound = 0
! search radius is rad
           DO irad = 1, rad
             numarcs = irad * 8
             p_angle = 360000 / numarcs
             DO c_angle=0,360000-p_angle,p_angle
               angle = c_angle / 1000.00
               j1 = y + INT(COS(angle*deg2rad)*irad)
               i1 = x + INT(SIN(angle*deg2rad)*irad)
               jx = j1
               ix = i1
               IF (jx .LE. 0 ) jx = 1
               IF (jx .GE. ldas%nr ) jx = ldas%nr
               IF (ix .LE. 0 ) ix = 1
               IF (ix .GT. ldas%nc ) ix = ldas%nc
               !write(*,'(f8.2,5i4)') angle,i1,j1,ix,jx,irad
!               IF(out(ix, jx) .NE. undef ) then
               IF(lofinalxy(ix,jx) ) then
                 ifound = ifound + 1
                 tmp=S42DTEMPTEMP(ix,jx)
!                 S4FLAG(ii)=1
               ENDIF
! number of points to find
               IF(ifound .GE. 1) EXIT
             ENDDO
             IF (ifound .GE. 1) EXIT
           ENDDO
           IF (ifound .GE. 1) THEN
             S4out(X+count) = tmp
!             s42dtemptemp(x,y)=s4out(x+count)
!             lofinalxy(x,y)=.true.
           ELSE
              S4out(X+count)=GRID(X,Y)%ETAPRECIP(II)
!             s42dtemptemp(x,y)=s4out(x+count)
!             lofinalxy(x,y)=.true.
           ENDIF
        ELSE
!        S42D(X,Y)=S42DTEMPTEMP(X,Y)
!        if ((x.eq.123).and.(y.eq.85)) print *,' in else, set to',S42DTEMP(X,Y)
!        S4FLAG(ii)=1
        ENDIF !IF MISSING DATA AT POINT
                S42D(X,Y) = (S4OUT(X+COUNT))/3600
!               MASKGRADS(X,Y)=0.0
!               IF(LOfinalxy(X,y)) MASKGRADS(X,Y)=1.0
!               MASKTEMP(X,Y) = LOfinalxy(x,y)
     ENDDO !NC
	count=count+nx
  ENDDO !NR

	do y=1,ldas%nr
	do x=1,ldas%nc
	lofinalxy(x,y)=.true.
               MASKGRADS(X,Y)=0.0
               IF(LOfinalxy(X,y)) MASKGRADS(X,Y)=1.0
               MASKTEMP(X,Y) = LOfinalxy(x,y)
	enddo
	enddo


!	count=0
!         DO Y = 1, NY
!            DO X = 1, NX
!               if ((lo(x+count)).or.(lohpd(x+count)).or.(locmor(x+count)))  &
!                lofinal(x+count)=.true.
!                S42D(X,Y) = (S4OUT(X+COUNT))/3600
!               MASKGRADS(X,Y)=0.0
!               IF(LOfinal(X+COUNT)) MASKGRADS(X,Y)=1.0
!               MASKTEMP(X,Y) = LOfinal(X+COUNT)
!            END DO
!            COUNT = COUNT + NX
!         END DO

      END IF

         !  WRITE THE HOLE-FILLED GRID BACK OUT INTO A 1D ARRAY

         COUNT = 0


!	         DO Y = 1, NY
!            DO X = 1, NX
!	briantemp(x,y)=GRID(X,Y)%ETAPRECIP(II)
!	enddo
!	enddo
!                 open (unit=45,file='out.bin',form='unformatted', &
!                status='unknown')
!                 write(45) briantemp
!                 close(45)
!                 stop

         DO Y = 1, NY
            DO X = 1, NX

               ! set hourly precip over mexico to NARR Precip for pre-cmorph period
               IF (((ldas%yr.eq.2002).and.(ldas%mo.lt.12)).or.(ldas%yr.lt.2002)) then
                  IF ((grid(x,y)%precipweight.lt.16).and.(y.lt.112)) THEN
!	print *,'x, y, precipweight,s42d,etaprecip',x,y, &
!       grid(x,y)%precipweight,S42D(X,Y),GRID(X,Y)%ETAPRECIP(II)
                     S42D(X,Y)=(grid(x,y)%precipweight/16.0)*S42D(X,Y) &
                               +((1.0-(grid(x,y)%precipweight/16.0)) &
                               *GRID(X,Y)%ETAPRECIP(II))
                  ENDIF
               ENDIF

               ! set hourly precip over canada to NARR Precip for whole period
               IF((grid(x,y)%precipweight.lt.16).and.(y.ge.112)) &
                   S42D(X,Y)=(grid(x,y)%precipweight/16.0)*S42D(X,Y) &
                   + ((1.0-(grid(x,y)%precipweight/16.0)) &
                   *GRID(X,Y)%ETAPRECIP(II))

               S4(II,X+COUNT) = S42D(X,Y)	

            END DO
            COUNT = COUNT+NX
         END DO


      DO R=1,LDAS%NR
         DO C=1,LDAS%NC
            ! IF(GRID(C,R)%FIMASK.EQ.0) S42D(C,R)=LDAS%UDEF
            ! IF (S4AVAIL(II)) THEN
            S42D24(II,C,R)=S42D(C,R)
            ! ENDIF
         ENDDO
      ENDDO

      !  DO R=1,LDAS%NR
      !     DO C=1,LDAS%NC
      !        IF (.NOT.S4AVAIL(II)) THEN
      !           S42D24(II,C,R)=GRID(C,R)%ETAPRECIP(II)
      !         ENDIF
      !     ENDDO
      !  ENDDO
      !  IF (.NOT.S4AVAIL(II)) print
      !  &  *,'putting model precip into s4 array'
      !  print *,'writing to 88'
      !  WRITE(88) (((S42D24(II,C,R)),C=1,LDAS%NC),R=1,LDAS%NR)
      !  close(88)
      !  print *,'stopping'
      !  stop

      !  IF(GRID(C,R)%FIMASK.EQ.0) S42DTEMP(C,R)=LDAS%UDEF
      !  IF(MASKGRADS(C,R).EQ.0) S42DTEMP(C,R)=LDAS%UDEF
      !  IF(GRID(C,R)%FIMASK.EQ.0) GRID(C,R)%ETAPRECIP(II)=LDAS%UDEF
      !  IF(GRID(C,R)%FIMASK.EQ.0) MASKGRADS(C,R)=LDAS%UDEF
      !  WRITE(88) (((MASKGRADS(C,R)),C=1,LDAS%NC),R=1,LDAS%NR)
      !  WRITE(88) (((GRID(C,R)%ETAPRECIP(II)),C=1,LDAS%NC),R=1,LDAS%NR)
      !  WRITE(88) (((S42DTEMP(C,R)),C=1,LDAS%NC),R=1,LDAS%NR)
      !  WRITE(88) (((S42D(C,R)),C=1,LDAS%NC),R=1,LDAS%NR)
      !  close(88)
      !  stop

      !  SUM HOURLY STAGE IV PRECIPS (WITH NO MISSING DATA) OVER 24H PERIOD

      DO N = 1, NLDAS
         SUM(N) = SUM(N)+S4(II,N)
      END DO

      !  open (unit=45,file='out.bin',form='unformatted',
      !  &  status='unknown')
      !  write(45) sum
      !  close(45)
      !  stop

   END DO

   !  open (unit=89,file='out.bin',form='unformatted',
   !  &  status='unknown')
   !  write(89) ((S4(II,N),n=1,nldas),II=1,24)
   !  close(89)
   !  stop

   !  CLOSE(88)

   !  Read in HIGGINS UNIFIED precip data
   !  Interpolate to 1/8th degree
   !  If no 1/4 degree points to work with, make 1/8th degree precip undefined
   IF (LDAS%UNIFIED.EQ.0) THEN
      UNIFIEDFLAG=0
   ELSE
      UNIFIEDFLAG=0
      OPEN (59, FILE=UNIFIEDNAME,FORM='UNFORMATTED',STATUS='OLD',ERR=76)
      READ (59) HIGGINSTEMP1D
      CLOSE(59)
      UNIFIEDFLAG=1
      76 CONTINUE
     
      PRISMFLAG=0
!  Read in HIGGINS UNIFIED PRISM precip data
      OPEN (59, FILE=PRISMNAME,FORM='UNFORMATTED',STATUS='OLD',ERR=77)
      READ (59) PRISMTEMP1D
      CLOSE(59)
     
! Read in PRISM precip mask
      OPEN (59, FILE='./BCS/N0.125/newmask.txt',STATUS='OLD')
   COUNT = 0
   DO Y = 1, NY
      DO X = 1, NX
         read(59,*) dummyint,dummyint,newmask(x+count)
      END DO
      COUNT = COUNT + NX
   END DO
	close(59)

! using both prism mask and more conservative precipweight mask
! make borders of prism data undefined due to lack of good data
	
   COUNT = 0
   DO Y = 1, NY
      DO X = 1, NX
         if ((grid(x,y)%precipweight.lt.16).or.(newmask(x+count).lt.1)) &
            PRISMTEMP1D(x+count)=LDAS%UDEF
      END DO
      COUNT = COUNT + NX
   END DO



!	print *,'prismtemp1d=',prismtemp1d(89753)
      PRISMFLAG=1
      77 CONTINUE
      
      IF (UNIFIEDFLAG.eq.1) THEN
         ! Interpolate Higgins Unified data

         DO I=1,200
            KGDSUNI(I)=0
         ENDDO

         KGDSUNI(1)=0
         KGDSUNI(2)=720
         KGDSUNI(3)=360
         KGDSUNI(4)=-89750
         KGDSUNI(5)=000
         KGDSUNI(6)=128
         KGDSUNI(7)=89750
         KGDSUNI(8)=360000
         KGDSUNI(9)=500
         KGDSUNI(10)=500
         KGDSUNI(11)=64
         KGDSUNI(12)=0
         KGDSUNI(13)=0
         KGDSUNI(14)=0
         KGDSUNI(15)=0
         KGDSUNI(16)=0
         KGDSUNI(17)=0
         KGDSUNI(18)=0
         KGDSUNI(19)=0
         KGDSUNI(20)=255
         KGDSUNI(21)=0
         KGDSUNI(22)=0
         IPUNI =  2

         DO I=1,20
          IPUNIOPT(I)=0
         ENDDO
         !  IPUNIOPT(1) =  4
         !  IPUNIOPT(2) = -1

         KMUNI = 1
         IBIUNI = 1

         DO I=1,NUNI
         LBUNI(I)=.TRUE.
         ENDDO

         CALL IPOLATES (IPUNI,IPUNIOPT,KGDSUNI,KGDSLDAS,NUNI,NLDAS,KMUNI,IBIUNI, &
                        LBUNI,HIGGINSTEMP1D,NO,RLAT,RLON,IBO,LO,HIGGINSLDAS1D,IRET)

      ENDIF !endif the unified flag = 1


      ! Transfer 1D data to 2D array and convert to mm/day from inches/day
      COUNT = 0
!        print *,'now prismunified=',PRISMUNIFIED(201,194)
!        print *,'unified flag=',UNIFIEDFLAG

      DO Y = 1, NY
         DO X = 1, NX
            if (UNIFIEDFLAG.eq.1) THEN
               IF (HIGGINSLDAS1D(X+COUNT).GE.0) THEN
                  HIGGINSUNIFIED(X,Y) = HIGGINSLDAS1D(X+COUNT)*25.4
               ENDIF
            endif

            if (PRISMFLAG.eq.1) THEN
               IF (PRISMTEMP1D(X+COUNT).GE.0) THEN
                  PRISMUNIFIED(X,Y)=PRISMTEMP1D(X+COUNT)*25.4
               ELSE
!                  PRISMUNIFIED(X,Y)=PRISMTEMP1D(X+COUNT)*25.4
                   PRISMUNIFIED(X,Y)=LDAS%UDEF
                  ! set prism to cpc
                  ! or stage2 or edas if undefined
                  ! PRISMUNIFIED(X,Y)=SUM(X+COUNT)
                  ! PRISMUNIFIED(X,Y)=SUM(X+COUNT)*3600.0
               ENDIF
            endif
            IF (UNIFIEDFLAG.EQ.1) THEN
               !  OVERWRITE HIGGINS UNIFIED WITH HIGGINS PRISM UNIFIED
               !  IF PRISM IS GREATER THAN OR EQUAL TO 0
               IF (PRISMUNIFIED(X,Y).GE.0) THEN
                  HIGGINSUNIFIED(X,Y)=PRISMUNIFIED(X,Y)
               ELSE
                  !  DO NOTHING, don't replace higgins precip with prism precip
               ENDIF

               IF ((GRID(X,Y)%FMASK.GE.1.0).AND.(HIGGINSUNIFIED(X,Y).NE.LDAS%UDEF)) THEN
                  IF (HIGGINSUNIFIED(X,Y).LT.0.0) THEN
                     write (*,'(A22,F11.8,A6,I3,A1,I3,A21)') &
                            'HIGGINS OUT OF RANGE (', &
                            HIGGINSUNIFIED(X,Y),') at (',c,',',r, &
                            ') CORRECTING TO 0'
                     HIGGINSUNIFIED(X,Y)=0.0
                  ENDIF

                  IF (HIGGINSUNIFIED(X,Y).GT.1000.0) THEN
                     write (*,'(A22,F11.8,A6,I3,A1,I3,A20)') &
                            'HIGGINS OUT OF RANGE (', &
                            HIGGINSUNIFIED(X,Y),') at (',c,',',r, &
                            ') CORRECTING TO 1000'
                     HIGGINSUNIFIED(X,Y)=1000.0
                  ENDIF

               ENDIF

            ENDIF !endif the unfied flag =1

         END DO

         COUNT = COUNT + NX

      END DO

   ENDIF  ! UNIFIED FLAG == 0

   !         COUNT = 0
   !         DO Y = 1, NHY
   !           DO X = 1, NHX
   !             HIGGINSTEMP2D(X,Y) = HIGGINSTEMP1D(X+COUNT)*25.4
   !           END DO
   !           COUNT = COUNT + NHX
   !         END DO


   !       do j=1,LDAS%NR
   !                 do i=1,LDAS%NC
   !                        HIGGINSUNIFIED(i,j)=0.0
   !                          ice=0
   !                          do li=1,4
   !                            ia=GRID(i,j)%FHLOCI(li)
   !                            jb=GRID(i,j)%FHLOCJ(li)
   !                             if(GRID(i,j)%FHDIST(li).LT.999) then
   !                              HIGGINSUNIFIED(i,j)=HIGGINSUNIFIED(i,j)+
   !     &                (GRID(i,j)%FHWHT(li)*HIGGINSTEMP2D(ia,jb))
   !                       else
   !                               ice=ice+1
   !                            endif
   !                          enddo
   !                          if(ice.eq.4) HIGGINSUNIFIED(i,j)=0.0
   !                 enddo
   !        enddo

   !  CALL GRIDCPC SUBROUTINE AND RETURN GRID OF CPC DAILY PREC TOTALS
   !    If CPCFLAG is set to 0 in the card (user does not
   !    wish to use CPC data) then set CPCFLAG to 0

   IF (LDAS%CPC.EQ.0) THEN
      CPCFLAG=0
   ELSE
         !  READ IN THE CPC mexico GRID AND INTERPOLATE TO LDAS GRID.
         !  Use quarter degree data for 2002 to present
         !  Use one degree data for 2001 and back
	 !  Use 0.5 degree global CPC precipitation data after 2010 
     IF (((ldas%yr.eq.2001).and.(ldas%mo.eq.12).and.(ldas%da.eq.31)).or.(ldas%yr.ge.2002)) THEN
      
      IF(ldas%yr.ge.2010) THEN       
       CALL GRIDCPCNEW (HNAMEmex,CPC,NLDAS,CPCFLAG)
      ELSE
       CALL GRIDCPCMEX (HNAMEmex,CPC,NLDAS,CPCFLAG)
       ENDIF 
     
      ELSE
       CALL GRIDCPC (HNAMEmex,CPC,NLDAS,CPCFLAG)
      ENDIF
!	print *,'mex gage is',cpc(89753)
!       print *,'cpcflag=',cpcflag
   ENDIF
      print *,'cpcflag=',cpcflag 
   !  Transfer CPC to X,Y grid from 1D grid

   COUNT = 0

   DO Y = 1, NY
      DO X = 1, NX
         CPCTMP(X,Y) = CPC(X+COUNT)
      END DO
      COUNT = COUNT + NX
   END DO

   ! Adjust GAGE Flag
   GAGEFLAG=0

   IF(PRISMFLAG.EQ.1) GAGEFLAG=1
   IF((LDAS%UNIFIED.EQ.2).AND.(UNIFIEDFLAG.NE.0)) GAGEFLAG=1
   IF((LDAS%CPC.EQ.2).AND.(CPCFLAG.NE.0)) GAGEFLAG=2

   IF((LDAS%UNIFIED.EQ.1).AND.(UNIFIEDFLAG.NE.0)) GAGEFLAG=1
   IF((LDAS%CPC.EQ.1).AND.(CPCFLAG.NE.0)) GAGEFLAG=2

   IF(CPCFLAG.EQ.0) GAGEFLAG=0

   !  Initialize Gage array to Undefined
   DO Y = 1, NY
      DO X = 1, NX
         GAGE(X,Y)=LDAS%UDEF
      END DO
   END DO

   !  Transfer CPC or HIGGINS UNIFIED data into Gage array
       !print *,'gageflag,higginsunified,cpc,cpcflag,cpctmp', &
       !gageflag,higginsunified(201,194),ldas%cpc,cpcflag, &
       !cpctmp(201,194)
       
   IF (GAGEFLAG.EQ.1) THEN
      DO Y = 1, NY
         DO X = 1, NX
            GAGE(X,Y)=HIGGINSUNIFIED(X,Y)
            IF ((HIGGINSUNIFIED(X,Y).LT.0).and.(LDAS%CPC.GT.0).and.(CPCFLAG.NE.0)) THEN
               GAGE(X,Y)=CPCTMP(X,Y)
            ENDIF
         END DO
      END DO

   ELSEIF (GAGEFLAG.EQ.2) THEN

      DO Y = 1, NY
         DO X = 1, NX
            GAGE(X,Y)=CPCTMP(X,Y)
         END DO
      END DO

   ELSE

      DO Y = 1, NY
         DO X = 1, NX
            GAGE(X,Y)=LDAS%UDEF
         END DO
      END DO

   ENDIF

   !  OVERWRITE GAGE WITH HIGGINS PRISM UNIFIED
   !  IF PRISM FLAG is 1 and prism IS GREATER THAN OR EQUAL TO 0
!	print *,'prismflag=',prismflag
!	print *,'prismunified=',prismunified(201,194)
!	print *,'before, gage=',gage(201,194)
   IF (PRISMFLAG.EQ.1) THEN
      count=0
      DO Y = 1, NY
         DO X = 1, NX
            IF (PRISMUNIFIED(x,y).GE.0) THEN
               count=count+1
               GAGE(x,y)=PRISMUNIFIED(x,y)
            ELSE
               !  DO NOTHING, don't replace gage precip with prism precip
            ENDIF
         ENDDO
      ENDDO
      if(count.eq.0) GAGEFLAG=0        ! no cpc conus precipitation
   ENDIF

!	print *,'after gage=',gage(201,194)

   !  print *,'opening',NAMES(1)(1:50)//'GAGE.OUT'
   !  OPEN (UNIT=88,FILE=NAMES(1)(1:50)//'GAGE.OUT',FORM=
   !  &  'UNFORMATTED')
   !  WRITE(88) (((GAGE(C,R)),C=1,LDAS%NC),R=1,LDAS%NR)
   !  CLOSE(88)


   DO I = 1, 24

!              open (unit=45,file='out.bin',form='unformatted', &
!             status='unknown')
!              write(45) ((float(GRID(X,Y)%FIMASK),x=1,464),y=1,224)
!              close(45)
!              stop

      !  IF HAVE STAGEIV DATA, BUT DO NOT HAVE 24HR GAGE DATA,
      !  THEN SET FINALPRECIP EQUAL TO THE UNWEIGHTED STAGEIV DATA

      IF ((GAGEFLAG.EQ.0).AND.(S4FLAG(I).EQ.1)) THEN
         !print *,'have s4 but not 24hr data'
         DO Y=1,NY
            DO X=1,NX
               FINALPRECIP(I,X,Y)=S42D24(I,X,Y)
               IF(GRID(X,Y)%FIMASK.EQ.0) FINALPRECIP(I,X,Y)=LDAS%UDEF
            ENDDO
         ENDDO

         CALL SYSTEM ('mkdir -p '//NAMES(I)(1:DIRLEN))

         OPEN(UNIT=68,FILE=NAMES(I),FORM='UNFORMATTED')

         do y =1,ny
            do x=1,nx

               WRITE(68) FINALPRECIP(I,X,Y)
               WRITE(68) S4SAVE2D(I,X,Y)

            enddo
         enddo
         CLOSE(68)

      ENDIF  !ENDIF GAGEFLAG.EQ.0, S4FLAG.EQ.1 loop


      !  IF DO NOT HAVE STAGEIV DATA AND DON'T HAVE 24HR GAGE DATA,
      !  THEN SET FINALPRECIP EQUAL TO THE UNWEIGHTED ETA-BASED DATA
      IF ((GAGEFLAG.EQ.0).AND.(S4FLAG(I).EQ.0)) THEN

         DO Y=1,NY
            DO X=1,NX
               FINALPRECIP(I,X,Y)=GRID(X,Y)%ETAPRECIP(I)
               IF(GRID(X,Y)%FIMASK.EQ.0) FINALPRECIP(I,X,Y)=LDAS%UDEF
            ENDDO
         ENDDO

         CALL SYSTEM ('mkdir -p '//NAMES(I)(1:DIRLEN))

         OPEN(UNIT=68,FILE=NAMES(I),FORM='UNFORMATTED')

         do y =1,ny
            do x=1,nx
               WRITE(68) FINALPRECIP(I,X,Y)
               WRITE(68) S4SAVE2D(I,X,Y)
            enddo
         enddo
         CLOSE(68)

      ENDIF  !ENDIF GAGEFLAG.EQ.0, S4FLAG.EQ.0 loop

   ENDDO ! I=1,24

   IF (GAGEFLAG.NE.0) THEN
      !  COMPUTE WEIGHT FOR EACH HOUR.  WATCH OUT FOR SINGULARITY!  FIND SUM
      !  OF HOURLY WEIGHTS.  SHOULD ALWAYS BE ZERO OR ONE!

      do n=1,nldas
         weightsum(n)=0.0
      enddo

      !        open (unit=45,file='out.bin',form='unformatted',
      !     &  status='unknown')
      !        write(45) sum
      !        close(45)
      !        stop

      DO I = 1, 24
         DO N = 1, NLDAS
            IF (SUM(N).GT.0) THEN
               WEIGHT(I,N) = S4(I,N)/SUM(N)
            ELSE
               WEIGHT(I,N) = 0
            END IF
            WEIGHTSUM(N)=WEIGHTSUM(N)+WEIGHT(I,N)
         END DO
      END DO

      !  PUT WEIGHTS AND WEIGHTSUM ON 2D ARRAYS FOR SEARCHING ALGORITHM BELOW.

      COUNT = 0

      DO Y = 1, NY
         DO X = 1, NX
            WTSUMTMP(X,Y) = WEIGHTSUM(X+COUNT)
            DO I = 1, 24
               WT2DTMP(I,X,Y)=WEIGHT(I,X+COUNT)
            END DO
         END DO
         COUNT = COUNT + NX
      END DO



      ! OPEN (UNIT=88,FILE='HIGGINS.OUT',FORM='UNFORMATTED')
      ! WRITE(88) (((CPCTMP(C,R)),C=1,LDAS%NC),R=1,LDAS%NR)
      ! CLOSE(88)

      ! SEARCH FOR OCCURRENCES WHERE WEIGHTSUM IS ZERO, BUT GAGE PRECIP IS
      ! NON-ZERO.  IN THESE INSTANCES, TAKE HOURLY WEIGHTS FROM NEAREST NEIGHBOR.
      ! IN THIS MANNER, THE SUM OF THE 24 HOURLY PRECIPS TO BE DERIVED BELOW
      ! WILL ALWAYS EQUAL THE GAGE 24H ANALYSIS AT ALL GRID POINTS.

!	print *,'replace missing hourdata'
  rad=20
  deg2rad = 0.01745329
  DO y=1, ldas%nr
     DO x=1, ldas%nc
! out2 new and out is orig array
!        out2(i,j) = out(i,j)
!        IF (out(i, j) .EQ. undef) THEN
        IF ((WTSUMTMP(X,Y).EQ.0).and.(GAGE(X,Y).NE.0)) THEN
!	print *,'x,y,wtsum is 0 but daily gage isnt',x,y
           ifound = 0
! search radius is rad
           DO irad = 1, rad
             numarcs = irad * 8
             p_angle = 360000 / numarcs
             DO c_angle=0,360000-p_angle,p_angle
               angle = c_angle / 1000.00
               j1 = y + INT(COS(angle*deg2rad)*irad)
               i1 = x + INT(SIN(angle*deg2rad)*irad)
               jx = j1
               ix = i1
               IF (jx .LE. 0 ) jx = 1
               IF (jx .GE. ldas%nr ) jx = ldas%nr
               IF (ix .LE. 0 ) ix = 1
               IF (ix .GT. ldas%nc ) ix = ldas%nc
               !write(*,'(f8.2,5i4)') angle,i1,j1,ix,jx,irad
!               IF(out(ix, jx) .NE. undef ) then
               IF(WTSUMTMP(ix, jx) .NE. 0 ) then
                 ifound = ifound + 1
                 do i=1,24
                 tmp24(i)=WT2DTMP(I,ix,jx)
                 enddo
               ENDIF
! number of points to find
               IF(ifound .GE. 1) EXIT
             ENDDO
             IF (ifound .GE. 1) EXIT
           ENDDO
           IF (ifound .GE. 1) THEN
                 do i=1,24
                  WT2D(I,X,Y) = tmp24(i)
	         enddo
           ELSE
              do i=1,24
              WT2D(I,X,Y) = 1.0/24.0
              enddo
           ENDIF
        ELSE
        do i=1,24
        WT2D(I,X,Y)=WT2DTMP(I,X,Y)
        enddo
        ENDIF !IF MISSING DATA AT POINT
     ENDDO !NC
  ENDDO !NR

!	print *,'finished with search'
!      DO Y = 1, NY
!         DO X = 1, NX
!            FOUND = .TRUE.
!            IF ((WTSUMTMP(X,Y).EQ.0).AND.(GAGE(X,Y).NE.0)) THEN
!               FOUND = .FALSE.
!               IF ((X.LE.(NX-BOUND)).AND.(X.GT.BOUND).AND. &
!                  (Y.LE.(NY-BOUND)).AND.(Y.GT.BOUND)) THEN
!
!                  DIST = 1
!                  DO WHILE ((DIST .LE. BOUND).AND.(.NOT. FOUND))
!                     YY = -1*DIST
!                     DO WHILE ((YY .LE. DIST).AND.(.NOT. FOUND))
!                        XX = -1*DIST
!                        DO WHILE ((XX .LE. DIST).AND.(.NOT. FOUND))
!                           IF (WTSUMTMP(X+XX,Y+YY).NE.0) THEN
!                              DO I = 1, 24
!                                 WT2D(I,X,Y) = WT2DTMP(I,X+XX,Y+YY)
!                              END DO
!                              FOUND = .TRUE.
!                           END IF
!                           XX = XX+1
!                        END DO
!                        YY = YY+1
!                     END DO
!                     DIST = DIST+1
!                  END DO
!
!               ELSE
!
!                  DO I = 1, 24
!                     WT2D(I,X,Y)=1.0/24.0
!                  END DO
!
!                  FOUND = .TRUE.
!
!               END IF
!
!            ELSE
!
!               DO I = 1, 24
!                  WT2D(I,X,Y)=WT2DTMP(I,X,Y)
!               END DO
!
!            END IF
!
!            IF (.NOT. FOUND) THEN
!               DO I = 1, 24
!                  WT2D(I,X,Y) = 1.0/24.0
!               END DO
!            END IF
!
!         END DO !X=1,NX
!
!      END DO !Y=1,NY

      !  Set FINALPRECIP (which holds the merged product)
      !  equal to the Gage precip times the appropriate hourly
      !  weight and divided by 3600 seconds to arrive at
      !  a rainfall rate.  This will only happen if GAGE 24hr
      !  product is available....and weights can be
      !  derived from either StageIV hourly or EDAS data

      DO I=1,24

         DO Y=1,NY

            DO X=1,NX

               IF(GRID(X,Y)%FIMASK.NE.0) FINALPRECIP(I,X,Y) = &
                    (GAGE(X,Y)*WT2D(I,X,Y))/3600.0
               IF(GRID(X,Y)%FIMASK.EQ.0) FINALPRECIP(I,X,Y) = LDAS%UDEF
               ! set precip over canada to NARR Precip for whole period
               IF((grid(x,y)%precipweight.lt.16).and.(y.ge.112)) &
!                    FINALPRECIP(I,X,Y)=S42D24(I,X,Y)
		FINALPRECIP(I,X,Y)=(grid(x,y)%precipweight/16.0)* &
                     ((GAGE(X,Y)*WT2D(I,X,Y))/3600.0) &
                   + (1.0-(grid(x,y)%precipweight/16.0)) &
                   *S42D24(I,X,Y)
               if (finalprecip(i,x,y).gt.(125.0/3600)) badprecip(i)=1
!	       if ((x.eq.70).and.(y.eq.80)) print *,'precip=',finalprecip(i,x,y)

!!!	        if (FINALPRECIP(I,X,Y).lt.(0.01/3600.)) FINALPRECIP(I,X,Y)=0.0

            ENDDO ! X=1,NX

         ENDDO ! Y=1,NY

	ENDDO ! enddo the i=1,24

! set finalprecip to narr for whole day if original value was over 125mm per hour

	do i=1,24
         if (badprecip(i).eq.1) flag=1
	enddo

	do i=1,24
	 if (flag.eq.1) then 
           DO Y=1,NY
            DO X=1,NX
	      finalprecip(i,x,y)=GRID(X,Y)%ETAPRECIP(I)
	    ENDDO
           ENDDO
         endif
	enddo




! write out files
      DO I=1,24

         !  open (unit=45,file='out.bin',form='unformatted',
         !  &  status='unknown')
         !  write(45)((grid(x,y)%FIMASK,x=1,464),y=1,224)
         !  close(45)
         !  stop

         !  Make directory structure for Merged precip product
         !  if it doesn't already exist.  Write out 24 merged
         !  precip product files.

         !DO Y=1,NY
         !  DO X=1,NX
         !  print *,'time,x,y,finalprecip=',i,x,y,FINALPRECIP(I,X,Y)
         !  print *,'time,x,y,s4save2d=',i,x,y,S4SAVE2D(I,X,Y)
         !  ENDDO
         !ENDDO

         CALL SYSTEM ('mkdir -p '//NAMES(I)(1:DIRLEN))

         !  OPEN(UNIT=68,FILE=NAMES(I),FORM='UNFORMATTED')
         !  WRITE(68) ((FINALPRECIP(I,X,Y),X=1,NX),Y=1,NY)
         !  WRITE(68) ((S4SAVE2D(I,X,Y),X=1,NX),Y=1,NY)
         !  CLOSE(68)

         !  OPEN(UNIT=68,FILE=NAMES(I),ACCESS='DIRECT',RECL=1)
         !  ire=0
         !  do y =1,ny
         !    do x=1,nx
         !      ire=ire+1
         !      WRITE(68,REC=ire) finalprecip(i,x,y)
         !      WRITE(68,REC=ire+(nx*ny)) finalprecip(i,x,y)
         !    enddo
         !  enddo
         !  CLOSE(68)

         OPEN(UNIT=68,FILE=NAMES(I),FORM='UNFORMATTED')

         do y =1,ny
            do x=1,nx
               WRITE(68) FINALPRECIP(I,X,Y)
               WRITE(68) S4SAVE2D(I,X,Y)
            enddo
         enddo

         CLOSE(68)

         !        OPEN(UNIT=68,FILE=NAMES(I)(1:44),status='unknown')
         !        do y =1,ny
         !      do x=1,nx
         !        WRITE(68,*) 'finalprecip=',i,x,y,FINALPRECIP(I,X,Y)
         !        WRITE(68,*) 's4save2d=',i,x,y,S4SAVE2D(I,X,Y)
         !      enddo
         !      enddo
         !        CLOSE(68)

         !      print*,'after write'
         !          DO Y=1,NY
         !         DO X=1,NX
         !        print *,'time,x,y,finalprecip=',i,x,y,FINALPRECIP(I,X,Y)
         !        print *,'time,x,y,s4save2d=',i,x,y,S4SAVE2D(I,X,Y)
         !        ENDDO
         !        ENDDO
      ENDDO ! enddo the i=1 to 24 loop

   ENDIF  !ENDIF GAGEFLAG.ne.0 loop

   DO I=1,24

      OPEN(UNIT=67,FILE=TRIM(NAMES(I)(1:42))//'.WT',form='unformatted',STATUS='UNKNOWN')
      WRITE(67) ((WT2D(I,X,Y),x=1,464),y=1,224)
      CLOSE(67)

      OPEN(UNIT=68,FILE=TRIM(NAMES(I)(1:42))//'.INFO',STATUS='UNKNOWN')
      WRITE(68,'(A7,I2)') 'UNIFIED    ',UNIFIEDFLAG
      WRITE(68,'(A7,I2)') 'CPC    ',CPCFLAG
      WRITE(68,'(A7,I2)') 'GAGE    ',GAGEFLAG
      WRITE(68,'(A7,I2)') 'STAGE4 ',S4FLAG(I)
      IF (S4FLAG(I).EQ.1) THEN
         WRITE(68,'(A7,A2)')'ETA    ','0'
      ELSEIF (S4FLAG(I).EQ.0) THEN
         WRITE(68,'(A7,A2)')'ETA    ','1'
      ENDIF

      CLOSE(68)

   ENDDO
  
   !    Diagnostic writes that can be deleted in the final
   !    version of the program

   !         OPEN (UNIT=66,FILE='WEIGHT.OUT',FORM='UNFORMATTED')
   !    DO I=1,24
   !         WRITE(66) ((GRID(X,Y)%ETAPRECIP(I),X=1,NX),Y=1,NY)
   !     WRITE(66) ((WT2D(I,X,Y),X=1,NX),Y=1,NY)
   !         WRITE(66) ((FINALPRECIP(I,X,Y),X=1,NX),Y=1,NY)
   !         WRITE(66) ((CPCTMP(X,Y),X=1,NX),Y=1,NY)
   !     WRITE(66) ((((CPCTMP(X,Y)*
   !     &  WT2D(I,X,Y))/3600.0),X=1,NX),Y=1,NY)
   !         WRITE(66) (((CPCTMP(X,Y)*WT2D(I,X,Y)),X=1,NX),Y=1,NY)

   !    ENDDO
   !         CLOSE(66)
   ! End Diagnostic writes

ENDIF  !End of the NEEDPRECIP=1 loop

!  Check to see if it is time (the next hour) to read in new
!  precipitation information from the merged
!  precipitation product files.  If it is time
!  to read in new data, then do so and replace
!  the precip forcing already in the forcing
!  array.  If not time, then overlay the data
!  that is already in the merged precip product
!  array on top of the data in the forcing array.

IF (LDAS%TIME.GT.LDAS%MERGETIME) THEN

   TYR=LDAS%YR
   TMO=LDAS%MO
   TDA=LDAS%DA
   THR=LDAS%HR
   TMN=LDAS%MN
   TSS=LDAS%SS

   !  If current minutes are not equal to zero, then
   !  we can use the next model hour to open up
   !  the required file.  If they equal 0, then we
   !  should use the current hour.

   IF (LDAS%MN.NE.0) THEN

      TTS=60*60
      CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,TTS)

   ENDIF

   LONGYEAR=TYR

   !  Assemble the required merged precip product filename

   OPEN(90,FILE='TEMP',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

   WRITE(90,28,REC=1)'/',TYR,TMO,'/',TMO,TDA,TYR,'.',THR,'.MERGE'
   READ(90,29,REC=1)FDIR2

   DO I=1,25
      IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
   ENDDO

   WRITE(90,26,REC=1) LDAS%MERGEH
   READ(90,22,REC=1) (FBASE1(I),I=1,80)

   C=0

   DO I=1,80

      IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1

   ENDDO

   WRITE(90,21,REC=1)(FBASE1(I),I=1,C), (FDIR2(I),I=1,25)
   READ(90,23,REC=1) TNAME
   CLOSE(90)

   !  Read in merged precip product
   OPEN(UNIT=68,FILE=TNAME,FORM='UNFORMATTED',IOSTAT=STATUS)

   IF(STATUS.EQ.0) THEN

      DO X=1,LDAS%NC
         DO Y=1,LDAS%NR
            GRID(X,Y)%MERGEPRECIP=0.0
            GRID(X,Y)%S4PRECIP=0.0
         ENDDO
      ENDDO

      do y =1,ny
         do x=1,nx
         !  READ(68) FINALPRECIP(I,X,Y)
         READ(68) GRID(X,Y)%MERGEPRECIP
         READ(68) GRID(X,Y)%S4PRECIP

         !  print *,'finalgrid=',x,y,GRID(X,Y)%MERGEPRECIP
         !  print *,'s4save=',x,y,s4save2d(i,x,y)
         enddo
      enddo

!           open (unit=45,file='out.bin',form='unformatted', &
!       status='unknown')
!           write(45)((grid(x,y)%MERGEPRECIP,x=1,464),y=1,224)
!           close(45)
!           stop


      !  READ(68) ((GRID(X,Y)%MERGEPRECIP,X=1,NX),Y=1,NY)
      !  READ(68) ((GRID(X,Y)%S4PRECIP,X=1,NX),Y=1,NY)
      CLOSE(68)

      !  Transfer merged product into the forcing array, overlaying
      !  it on top of the EDAS information (if present already in
      !  the forcing array).

      DO X=1,LDAS%NC
         DO Y=1,LDAS%NR
            !  IF (LDAS%PRECIPMASK.EQ.0) THEN
            if(GRID(X,Y)%forcing(4).ne.0.0.and. &
               GRID(X,Y)%forcing(4).ne.LDAS%UDEF.and. &
               GRID(X,Y)%forcing(5).ne.LDAS%UDEF) then

               xrat=GRID(X,Y)%forcing(5)/GRID(X,Y)%forcing(4)

               if(xrat.gt.1.0) xrat=1.0
               if(xrat.lt.0.0) xrat=0.0
            else
                xrat=0.0
            endif

            GRID(X,Y)%forcing(4)=GRID(X,Y)%MERGEPRECIP


            GRID(X,Y)%forcing(5)=(xrat)* GRID(X,Y)%forcing(4)

            !           ENDIF
            !           IF (LDAS%PRECIPMASK.EQ.1) THEN
            !            if(GRID(X,Y)%forcing(4).ne.0.0.and.
            !     &        GRID(X,Y)%forcing(4).ne.LDAS%UDEF.and.
            !     &        GRID(X,Y)%forcing(5).ne.LDAS%UDEF) then
            !              xrat=GRID(X,Y)%forcing(5)/
            !     &        GRID(X,Y)%forcing(4)
            !              if(xrat.gt.1.0) xrat=1.0
            !              if(xrat.lt.0.0) xrat=0.0
            !            else
            !              xrat=0.0
            !            endif
            !            GRID(X,Y)%forcing(4)=(
            !     &       ((grid(x,y)%precipweight/16.0)*
            !     &       GRID(X,Y)%MERGEPRECIP)
            !     &       +( (1.0-(grid(x,y)%precipweight/16.0))*
            !     &       GRID(X,Y)%forcing(4)) )
            !            GRID(X,Y)%forcing(5)=(xrat)* GRID(X,Y)%forcing(4)
            !           ENDIF

            IF(GRID(X,Y)%FIMASK.EQ.0) GRID(X,Y)%forcing(4)=LDAS%UDEF

         ENDDO ! Y=1,NR

      ENDDO ! X=1,NC
!           open (unit=45,file='out.bin',form='unformatted', &
!       status='unknown')
!           write(45)((grid(x,y)%forcing(4),x=1,464),y=1,224)
!           close(45)
!           stop


      CALL TICK(TTIME,TDOY,TGMT,LONGYEAR,TMO,TDA,THR,ZEROI1,ZEROI2,ZEROI3)

      !  Set the time of the current merged precip product holdings
      !  to the new time

      LDAS%MERGETIME=TTIME

   ENDIF  ! STATUS=0

   !  If it is not time to read in new data, then
   !  simply overlay the old data once again on
   !  top of whatever precip may already be in the forcing
   !  array.

   ELSE

      DO X=1,LDAS%NC
         DO Y=1,LDAS%NR
            !  IF (LDAS%PRECIPMASK.EQ.0) THEN
            if(GRID(X,Y)%forcing(4).ne.0.0.and. &
               GRID(X,Y)%forcing(4).ne.LDAS%UDEF.and. &
               GRID(X,Y)%forcing(5).ne.LDAS%UDEF) then

               xrat=GRID(X,Y)%forcing(5)/GRID(X,Y)%forcing(4)

               if(xrat.gt.1.0) xrat=1.0
               if(xrat.lt.0.0) xrat=0.0

            else
               xrat=0.0
            endif
            GRID(X,Y)%forcing(4)=GRID(X,Y)%MERGEPRECIP
            GRID(X,Y)%forcing(5)=(xrat)* GRID(X,Y)%forcing(4)


            !           ENDIF
            !           IF (LDAS%PRECIPMASK.EQ.1) THEN
            !            if(GRID(X,Y)%forcing(4).ne.0.0.and.
            !     &        GRID(X,Y)%forcing(4).ne.LDAS%UDEF.and.
            !     &        GRID(X,Y)%forcing(5).ne.LDAS%UDEF) then
            !              xrat=GRID(X,Y)%forcing(5)/
            !     &        GRID(X,Y)%forcing(4)
            !              if(xrat.gt.1.0) xrat=1.0
            !              if(xrat.lt.0.0) xrat=0.0
            !            else
            !              xrat=0.0
            !            endif
            !            GRID(X,Y)%forcing(4)=(
            !     &       ((grid(x,y)%precipweight/16.0)*
            !     &       GRID(X,Y)%MERGEPRECIP)
            !     &       +( (1.0-(grid(x,y)%precipweight/16.0))*
            !     &       GRID(X,Y)%forcing(4)) )
            !            GRID(X,Y)%forcing(5)=(xrat)* GRID(X,Y)%forcing(4)
            !           ENDIF

            IF(GRID(X,Y)%FIMASK.EQ.0) GRID(X,Y)%forcing(4)=LDAS%UDEF

         ENDDO
      ENDDO

   ENDIF  ! LDAS%TIME.GT.LDAS%MERGETIME
   !         WRITE(88) (((GRID(c,r)%MERGEPRECIP),C=1,LDAS%NC),R=1,LDAS%NR)
   !         CLOSE(88)

!           open (unit=45,file='out.bin',form='unformatted', &
!       status='unknown')
!           write(45)((grid(x,y)%forcing(4),x=1,464),y=1,224)
!           close(45)
!           stop

!        stop

   RETURN

END


