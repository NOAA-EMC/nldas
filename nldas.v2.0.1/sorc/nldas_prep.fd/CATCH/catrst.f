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
! catrst.f: 
!
! DESCRIPTION:
!  This program reads and writes restart files for the catchment model.  
!  This includes all relevant water/energy storages.
!
! REVISION HISTORY:
!  4  Apr 2000: Jeffrey Walker; Initial code
!  19 Jan 2001: Brian Cosgrove; Added CLOSE statement
!=========================================================================
! RESTART FILE FORMAT(fortran sequential binary):
!  YR,MO,DA,HR,MN,SS !Restart time
!  CAT%TC        !Canopy temperature in the 3 catchment fractions [K]
!  CAT%TSURF     !Mean surface temperature of the catchment [K]
!  CAT%TSNOW     !Snow surface temperature of the catchment [K]
!  CAT%QA        !Specific humidity in the canopy air in the 3 fractons [-]
!  CAT%INT       !Canopy interception reservoir storage [mm]
!  CAT%GHT       !Heat content of the 6 soil layers [J.m-2] ??
!  CAT%WESN      !Water equivalent of the 3 snow layers [mm]
!  CAT%HTSN      !Heat content of the 3 snow layers [J.m-2] ??
!  CAT%SNDZ      !Thickness of the 3 snow layers [mm]
!  CAT%CATDEF    !Catchment deficit [mm]
!  CAT%RZEXC     !Root zone excess [mm]
!  CAT%SRFEXC    !Surface excess [mm]
!=========================================================================

      SUBROUTINE CATRST(RW,LDAS,CAT)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE cat_module       ! Catchment model variables
      IMPLICIT NONE      
      type (ldasdec) LDAS              
      type (catdec)  CAT(LDAS%NCATM)

!=== Local Variables =====================================================
      INTEGER :: RW            ! 1=read restart, 2=write restart
      INTEGER :: C,I ! Loop counters

!=== Temporary tile space transfer files (different than in ldas_module)
      INTEGER :: YR,MO,DA,HR,MN,SS  !Time variables
      INTEGER :: NCATM

      REAL :: DZ(6)

      CHARACTER*80 FILEN,MKFYRMO
      CHARACTER*1  FNAME(80),FBASE(80),FSUBS(80),FMKDIR(80)
      CHARACTER*1  FTIME(10),FYRMODIR(80)
      DATA DZ/0.0988,0.1952,0.3859,0.7626,1.5071,10.0/

!=== End Variable Definition =============================================

!=== Read Active Archive File ============================================

      IF ((RW.EQ.1 .AND. LDAS%CAT_IC.EQ.1) .OR.
     &   (RW.EQ.1.AND.LDAS%STARTCODE.EQ.1))THEN

       OPEN (40,FILE=LDAS%CAT_RFILE,FORM='unformatted')
       READ(40) YR,MO,DA,HR,MN,SS,NCATM  !Time

       IF (NCATM .NE. LDAS%NCATM) THEN
        WRITE(*,*) 'Catchment Model Restart File Mismatch'
        STOP
       END IF
       READ(40) CAT%TC(1)
       READ(40) CAT%TC(2)
       READ(40) CAT%TC(3)
       READ(40) CAT%TSURF
       READ(40) CAT%TSNOW
       READ(40) CAT%QA(1)
       READ(40) CAT%QA(2)
       READ(40) CAT%QA(3)
       READ(40) CAT%INT
       READ(40) CAT%GHT(1)
       READ(40) CAT%GHT(2)
       READ(40) CAT%GHT(3)
       READ(40) CAT%GHT(4)
       READ(40) CAT%GHT(5)
       READ(40) CAT%GHT(6)
       READ(40) CAT%WESN(1)
       READ(40) CAT%WESN(2)
       READ(40) CAT%WESN(3)
       READ(40) CAT%HTSN(1)
       READ(40) CAT%HTSN(2)
       READ(40) CAT%HTSN(3)
       READ(40) CAT%SNDZ(1)
       READ(40) CAT%SNDZ(2)
       READ(40) CAT%SNDZ(3)
       READ(40) CAT%CATDEF
       READ(40) CAT%RZEXC
       READ(40) CAT%SRFEXC
       CLOSE(40)

       WRITE(*,*)'Catchment Restart File Read: ',LDAS%CAT_RFILE
       WRITE(*,*)'Catchment Restart File Time: ',YR,MO,DA,HR,MN,SS

! Write to runtime diagnostics file

       WRITE(79,*)'Catchment Restart File Read: ',LDAS%CAT_RFILE
       WRITE(79,*)'Catchment Restart File Time: ',YR,MO,DA,HR,MN,SS

! Establish Model Restart Time  

       IF(LDAS%STARTCODE.EQ.1)THEN
        LDAS%YR=YR
        LDAS%MO=MO 
        LDAS%DA=DA
        LDAS%HR=HR
        LDAS%MN=MN
        LDAS%SS=SS
        CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,YR,MO,DA,
     &                 HR,MN,SS) 
        LDAS%KOS_STIME = LDAS%TIME
        WRITE(*,*)'Catchment Restart File Time Used: ',LDAS%CAT_RFILE
       ENDIF

      ENDIF !RW option 1

!=== In the case of no previously existing restart, the user can specify
!=== the IC's in the card file
      IF(RW.EQ.1 .AND. LDAS%CAT_IC.EQ.3)THEN !Card File Specified Restart
       CAT%TC(1)=LDAS%CAT_IT
       CAT%TC(2)=LDAS%CAT_IT
       CAT%TC(3)=LDAS%CAT_IT
       CAT%TSURF=LDAS%CAT_IT
       CAT%TSNOW=LDAS%CAT_IT
       CAT%QA(1)=0.0002
       CAT%QA(2)=0.0002
       CAT%QA(3)=0.0002
       CAT%INT=0.0
       IF (LDAS%CAT_IT.LT.273.16) THEN ! Soil is frozen
        CAT%GHT(1)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(1)+4.6E5*DZ(1))-
     &              751.5E5*DZ(1)
        CAT%GHT(2)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(2)+4.6E5*DZ(2))-
     &              751.5E5*DZ(2)
        CAT%GHT(3)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(3)+4.6E5*DZ(3))-
     &              751.5E5*DZ(3)
        CAT%GHT(4)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(4)+4.6E5*DZ(4))-
     &              751.5E5*DZ(4)
        CAT%GHT(5)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(5)+4.6E5*DZ(5))-
     &              751.5E5*DZ(5)
        CAT%GHT(6)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(6)+4.6E5*DZ(6))-
     &              751.5E5*DZ(6)
       ELSE ! Soil is thawed
        CAT%GHT(1)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(1)+9.4E5*DZ(1))
        CAT%GHT(2)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(2)+9.4E5*DZ(2))
        CAT%GHT(3)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(3)+9.4E5*DZ(3))
        CAT%GHT(4)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(4)+9.4E5*DZ(4))
        CAT%GHT(5)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(5)+9.4E5*DZ(5))
        CAT%GHT(6)=(LDAS%CAT_IT-273.16)*(13.2E5*DZ(6)+9.4E5*DZ(6))
       END IF

C	Uniform soil temperature initialization of 1 degree C 

c       cat%ght(1)=2.4e6*.55*dz(1)+.5*.45*dz(1)*4.185e6
c       cat%ght(2)=2.4e6*.55*dz(2)+.5*.45*dz(2)*4.185e6
c       cat%ght(3)=2.4e6*.55*dz(3)+.5*.45*dz(3)*4.185e6
c       cat%ght(4)=2.4e6*.55*dz(4)+.5*.45*dz(4)*4.185e6
c       cat%ght(5)=2.4e6*.55*dz(5)+.5*.45*dz(5)*4.185e6
c       cat%ght(6)=2.4e6*.55*dz(6)+.5*.45*dz(6)*4.185e6

       CAT%WESN(1)=0.0
       CAT%WESN(2)=0.0
       CAT%WESN(3)=0.0
       CAT%HTSN(1)=0.0
       CAT%HTSN(2)=0.0
       CAT%HTSN(3)=0.0
       CAT%SNDZ(1)=0.0
       CAT%SNDZ(2)=0.0
       CAT%SNDZ(3)=0.0
       CAT%CATDEF=MAX(CAT%POROS-LDAS%CAT_ISM,0.05)*CAT%ZDPTH(3)
       CAT%RZEXC=0.0
       CAT%SRFEXC=0.0

       WRITE(*,*) 'Catchment Initialized SM/T=',LDAS%CAT_ISM,LDAS%CAT_IT
       WRITE(79,*)'Catchment Initialized SM/T=',LDAS%CAT_ISM,LDAS%CAT_IT
      ENDIF

!=== Set starttime to ldas.crd Stime 
      IF(RW.EQ.1 .AND. LDAS%STARTCODE.EQ.3)THEN 
       LDAS%YR=LDAS%SYR
       LDAS%MO=LDAS%SMO 
       LDAS%DA=LDAS%SDA
       LDAS%HR=LDAS%SHR
       LDAS%MN=LDAS%SMN
       LDAS%SS=LDAS%SSS
       CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,LDAS%YR,
     &                LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS) 
       WRITE(*,*)'Using ldas.crd start time ',LDAS%TIME
       WRITE(79,*)'Using ldas.crd start time ',LDAS%TIME
      ENDIF

!=== Restart Writing (2 file are written - active and archive)
      IF(RW.EQ.2)THEN
       OPEN(40,FILE=LDAS%CAT_RFILE,FORM='unformatted') !Active archive restart

       WRITE(40) LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS,
     1           LDAS%NCATM
       WRITE(40) CAT%TC(1)
       WRITE(40) CAT%TC(2)
       WRITE(40) CAT%TC(3)
       WRITE(40) CAT%TSURF
       WRITE(40) CAT%TSNOW
       WRITE(40) CAT%QA(1)
       WRITE(40) CAT%QA(2)
       WRITE(40) CAT%QA(3)
       WRITE(40) CAT%INT
       WRITE(40) CAT%GHT(1)
       WRITE(40) CAT%GHT(2)
       WRITE(40) CAT%GHT(3)
       WRITE(40) CAT%GHT(4)
       WRITE(40) CAT%GHT(5)
       WRITE(40) CAT%GHT(6)
       WRITE(40) CAT%WESN(1)
       WRITE(40) CAT%WESN(2)
       WRITE(40) CAT%WESN(3)
       WRITE(40) CAT%HTSN(1)
       WRITE(40) CAT%HTSN(2)
       WRITE(40) CAT%HTSN(3)
       WRITE(40) CAT%SNDZ(1)
       WRITE(40) CAT%SNDZ(2)
       WRITE(40) CAT%SNDZ(3)
       WRITE(40) CAT%CATDEF
       WRITE(40) CAT%RZEXC
       WRITE(40) CAT%SRFEXC
       CLOSE(40)
   
       WRITE(*,*)'Catchment Active Restart Written: ',LDAS%CAT_RFILE
       WRITE(79,*)'Catchment Active Restart Written: ',LDAS%CAT_RFILE

!=== Now write the archived restart file
 91    FORMAT(A4,I3,A5,I4,I2,A7,I3,A1)
 92    FORMAT(80A1)
 93    FORMAT(A80)
 94    FORMAT(I4,I2,I2,I2)
 95    FORMAT(10A1)
 96    FORMAT(A40)
 98    FORMAT(A4,I3,A5,I4,I2)
100    FORMAT(A9)
101    FORMAT(A7)
       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
       WRITE(90,94,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
       READ(90,95,REC=1)FTIME
       DO I=1,10
        IF(FTIME(I).EQ.(' '))FTIME(I)='0'
       ENDDO
       WRITE(90,91,REC=1)'/EXP',LDAS%EXPCODE,'/CAT/',LDAS%YR,
     &  LDAS%MO,'/LDAS.E',LDAS%EXPCODE,'.'
       READ(90,92,REC=1) (FNAME(I),I=1,29)
       DO I=1,29
        IF(FNAME(I).EQ.(' '))FNAME(I)='0'
       ENDDO

       WRITE(90,100,REC=1)'mkdir -p '
       READ(90,92,REC=1)(FMKDIR(I),I=1,9)

       WRITE(90,98,REC=1)'/EXP',LDAS%EXPCODE,'/CAT/',
     &   LDAS%YR,LDAS%MO
       READ(90,92,REC=1) (FYRMODIR(I),I=1,18)
       DO I=1,18
        IF(FYRMODIR(I).EQ.(' '))FYRMODIR(I)='0'
       ENDDO

       WRITE(90,101,REC=1)'.CATrst'
       READ(90,92,REC=1) (FSUBS(I),I=1,7)

       WRITE(90,96,REC=1) LDAS%ODIR                       
       READ(90,92,REC=1) (FBASE(I),I=1,80)
       C=0
       DO I=1,80
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO
       WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FNAME(I),I=1,29),
     &                   (FTIME(I),I=1,10),(FSUBS(I),I=1,7 ) 
       READ(90,93,REC=1)FILEN
 
       WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &   (FYRMODIR(I),I=1,18)
       READ(90,93,REC=1)MKFYRMO
       CLOSE(90)
!== Archive File Name Generation Complete
!== Make the directories for the Catchment restart file
       IF(LDAS%GMT.LE.LDAS%WRITEINTCT)THEN             
        CALL SYSTEM(MKFYRMO)
       ENDIF

       OPEN(40,FILE=FILEN,status='unknown',FORM='unformatted')

       WRITE(40) LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS,
     1           LDAS%NCATM
       WRITE(40) CAT%TC(1)
       WRITE(40) CAT%TC(2)
       WRITE(40) CAT%TC(3)
       WRITE(40) CAT%TSURF
       WRITE(40) CAT%TSNOW
       WRITE(40) CAT%QA(1)
       WRITE(40) CAT%QA(2)
       WRITE(40) CAT%QA(3)
       WRITE(40) CAT%INT
       WRITE(40) CAT%GHT(1)
       WRITE(40) CAT%GHT(2)
       WRITE(40) CAT%GHT(3)
       WRITE(40) CAT%GHT(4)
       WRITE(40) CAT%GHT(5)
       WRITE(40) CAT%GHT(6)
       WRITE(40) CAT%WESN(1)
       WRITE(40) CAT%WESN(2)
       WRITE(40) CAT%WESN(3)
       WRITE(40) CAT%HTSN(1)
       WRITE(40) CAT%HTSN(2)
       WRITE(40) CAT%HTSN(3)
       WRITE(40) CAT%SNDZ(1)
       WRITE(40) CAT%SNDZ(2)
       WRITE(40) CAT%SNDZ(3)
       WRITE(40) CAT%CATDEF
       WRITE(40) CAT%RZEXC
       WRITE(40) CAT%SRFEXC
       CLOSE(40) 

       WRITE(*,*)'Catchment Archive Restart Written: ',FILEN
       WRITE(79,*)'Catchment Archive Restart Written: ',FILEN

      ENDIF

      Return
      End






