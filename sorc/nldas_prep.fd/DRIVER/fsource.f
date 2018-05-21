!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale
!  D                                      L  Land Modeling and Data
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! fsource.f:
!
! DESCRIPTION:
!  Outputs forcing status, as follows:
!
!
!  TIME                ETA          NCEP RAD  PRECIP
!                      3A 3F-A 6F-A 1A PK NE ETA SC SS S4
!  01/01/2000 00:00:00 1  0-0  0-0  1  1  0  0   0  0  1     (Ideal Forcing)
!  
!  
!  0=data not available
!  1=data available and used
!  (-1) to (-10)=data from N days in the past used
!
!  ETA-3A:    3hr-EDAS (Analysis)
!  ETA-3F:    3hr-Eta  (Forecast)
!  ETA-3F:    3hr-Eta (Forecast amount)
!  ETA-6F:    6hr-Eta  (Forecast)
!  ETA-6F:    6hr-Eta  (Forecast amount)
!  NCEP-1A:   1hr-NCEP (Official LDAS Forcing)
!  RAD-PK:    1hr Pinker-GOES Radiation
!  RAD-NE:    1hr NESDIS-GOES Radiation
!  PRECIP-ETA: 1hr Uncorrected ETA Precip (No 24hr data available)
!  PRECIP-S4: 1hr Uncorrected Stage-4 Precip (No 24hr data available)
!  PRECIP-ECH: 1hr ETA Precip Corrected With 24hr HGGINS CPC gage-only
!  PRECIP-ECU: 1hr ETA Precip Corrected With 24hr UNIFIED gage-only
!  PRECIP-SCH: 1hr Stage-4 Precip Corrected With 24hr HIGGINSCPC gage-only
!  PRECIP-SCU: 1hr Stage-4 Precip Corrected With 24hr UNIFIED gage-only
!
! REVISION HISTORY:
! 20  Dec 1999:  Paul Houser; Initial code
!  8  Mar 2000: Brian Cosgrove; Major revisions to deal with precip sources
! 27  Feb 2001: Brian Cosgrove; Revisions to add in catchment forcing data
!=========================================================================

      SUBROUTINE FSOURCE(LDAS)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      type (ldasdec) LDAS              

!=== Local variables =====================================================
      INTEGER I,C       !Local counting variable
      INTEGER PRECIP1,PRECIP2,PRECIP3,PRECIP4,PRECIP5
      INTEGER TYR,TMO,TDA,THR,TMN,TSS,TTS,TDOY,STATUS
      REAL *8 TTIME
      REAL TGMT
      CHARACTER *7 DUMCHAR
      CHARACTER TNAME*80
      CHARACTER*1 FBASE1(80),FDIR2(24)

!=== End Variable List ===================================================
C	Initialize precipitation information variables to 0
       LDAS%FSOURCE(9)=0
       LDAS%FSOURCE(10)=0
       LDAS%FSOURCE(11)=0
       LDAS%FSOURCE(12)=0
       LDAS%FSOURCE(13)=0
       LDAS%FSOURCE(14)=0


C       Form Filename of the Preciptation Information File

        TYR=LDAS%YR
        TMO=LDAS%MO
        TDA=LDAS%DA
        THR=LDAS%HR
        TMN=0
        TSS=0
	IF (LDAS%MN.NE.0) THEN
         CALL TICK(TTIME,TDOY,TGMT,TYR,TMO,TDA,THR,TMN,TSS,3600)
	ENDIF

 21   FORMAT(80A1)
 22   FORMAT(40A1)
 23   FORMAT(A80)
 26   FORMAT(A40)
 28   FORMAT(A1,I4,I2,A1,I2,I2,I4,A1,I2,A5)
 29   FORMAT(24A1)
      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
      WRITE(90,28,REC=1)'/',TYR,TMO,'/',TMO,TDA,TYR,
     &  '.',THR,'.INFO'
      READ(90,29,REC=1)FDIR2
      DO I=1,24
       IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO
      WRITE(90,26,REC=1) LDAS%MERGEH
      READ(90,22,REC=1) (FBASE1(I),I=1,40)
      C=0
      DO I=1,40
       IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO
      WRITE(90,21,REC=1)(FBASE1(I),I=1,C), (FDIR2(I),I=1,24)
      READ(90,23,REC=1)TNAME
      CLOSE(90)

      OPEN (UNIT=23,FILE=TNAME,STATUS='OLD',IOSTAT=STATUS)
C	Set fsource precip flags manually if unable to open
C	precipitation information file.  If file open is
C	successful, read in CPC, Stage4 and ETA precipitation
C	status.
      IF ((STATUS.NE.0).OR.(LDAS%PRECSOR(1).EQ.0))THEN
       LDAS%FSOURCE(9)=1
       LDAS%FSOURCE(10)=0
       LDAS%FSOURCE(11)=0
       LDAS%FSOURCE(12)=0
       LDAS%FSOURCE(13)=0
       LDAS%FSOURCE(14)=0
      ELSE
       READ(23,'(A7,I2)') dumchar,precip1
       READ(23,'(A7,I2)') dumchar,precip2
       READ(23,'(A7,I2)') dumchar,precip3
       READ(23,'(A7,I2)') dumchar,precip4
       READ(23,'(A7,I2)') dumchar,precip5

      ENDIF

      IF ((precip3.eq.1).and.(precip4.eq.1)) LDAS%FSOURCE(14)=1
      IF ((precip3.eq.2).and.(precip4.eq.1)) LDAS%FSOURCE(13)=1
      IF ((precip3.eq.0).and.(precip4.eq.1)) LDAS%FSOURCE(12)=1
      IF ((precip3.eq.1).and.(precip4.eq.0).and.(precip5.eq.1))
     &  LDAS%FSOURCE(11)=1
      IF ((precip3.eq.2).and.(precip4.eq.0).and.(precip5.eq.1))
     &  LDAS%FSOURCE(10)=1
      IF ((precip3.eq.0).and.(precip4.eq.0).and.(precip5.eq.1))
     &  LDAS%FSOURCE(9)=1

!      WRITE(*,24)
!     1 LDAS%MO,'/',LDAS%DA,'/',LDAS%YR,LDAS%HR,':',LDAS%MN,':',LDAS%SS,
!     2 (LDAS%FSOURCE(I),I=1,12)

      WRITE(78,24) 
     1 LDAS%MO,'/',LDAS%DA,'/',LDAS%YR,LDAS%HR,':',LDAS%MN,':',LDAS%SS,
     2 (LDAS%FSOURCE(I),I=1,14)
 24   FORMAT(I2,A1,I2,A1,I4,1X,I2,A1,I2,A1,I2,2X,
     &   I3,2(1X,I3,'-',I1),3(I3),(I4),(I4),(I4),(I3),2(I4))

	CLOSE(23)
      RETURN
      END
      
      
