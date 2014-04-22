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
! readcard.f90:
!
! DESCRIPTION:
!  Reads in LDAS run specifics from ldas.crd
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial code
!  4  Apr 2000: Jeffrey Walker; Added catchment model output interval
!  11 Apr 2000: Brian Cosgrove; Added Elevation correction and Forcing
!               Mask read statements
!  6  Jun 2000: Jon Radakovich; Updated for new version of CLM
!  23 Feb 2001: Urszula Jambor; Added GEOS or GDAS forcing option
!  27 Mar 2001: Jon Gottschalck; Revision of subroutine by implementing namelists
!  05 Sep 2001: Brian Cosgrove; Altered forcing logfile output to include
!               more precip types
!  04 Feb 2002: Jon Gottschalck; Added section to set to Koster tilespace files if necessary
!  15 Apr 2002: Urszula Jambor; Added ECMWF forcing options, also
!               adding 1 & 1/2 degree GLDAS domain options.
!  28 Apr 2002: Kristi Arsenault; Added NOAH LSM code
!  30 Jul 2002: Jon Gottschalck; Added code to handle additional global observed precip sources
!  25 Sep 2002: Jon Gottschalck; Added option to use MODIS UMD vegetation classification
!  01 Oct 2002: Jon Gottschalck; Added code for MODIS LAI data
!  20 Nov 2002: Jon Radakovich; Added startcode=6 (restarting a BC run)
!  11 Dec 2002: Urszula Jambor; Added initialization of AGRMET time variables
!  05 Feb 2003: Jon Gottschalck; Added CLM LSM, v 2.0 variables and program flow
!=========================================================================

      SUBROUTINE READCARD(LDAS)

      USE ldas_module
      IMPLICIT NONE				  
      type (ldasdec) LDAS

!=== Local Variables =====================================================
      INTEGER :: LD                      ! Last day of month
      INTEGER :: I                       ! Loop Counter
      CHARACTER*80 :: NAME               ! File Name to open

      NAMELIST /initial/ldas
      NAMELIST /driver/ldas
      NAMELIST /ldas_run_inputs/ldas
      NAMELIST /gldas1_4/ldas/nldas/ldas/gldas2x2_5/ldas
      NAMELIST /gldas_1/ldas/gldas1_2/ldas
      NAMELIST /mos/ldas/clm1/ldas/cat/ldas/noah/ldas/clm2/ldas
      NAMELIST /gdas/ldas/geos/ldas/eta/ldas/ncep/ldas/nasa/ldas/catch/ldas
      NAMELIST /ecmwf/ldas/recmwf/ldas

!=== End Variable Definition =============================================

      open(10,file='ldas.crd',form='formatted',status='old')


!=== Reading in parameters that need to be initialized to avoid any problems later
      READ(unit=10,NML=initial)

!=== Reading in parameters that are used by all LDAS runs
      READ(unit=10,NML=driver)
      READ(unit=10,NML=ldas_run_inputs)

!=== Open runtime diagnostics file
      CALL OPENFILE(NAME,LDAS%ODIR,LDAS%EXPCODE,LDAS%DFILE)
	print *,'name=',name
      IF(LDAS%STARTCODE.EQ.1.OR.LDAS%STARTCODE.EQ.6)THEN
       OPEN(79,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
      ELSE
       OPEN(79,FILE=NAME,FORM='FORMATTED',STATUS='REPLACE')
      ENDIF
      
!=== Open runtime diagnostics file and write header
      CALL OPENFILE(NAME,LDAS%ODIR,LDAS%EXPCODE,LDAS%FFILE)
      IF(LDAS%STARTCODE.EQ.1.OR.LDAS%STARTCODE.EQ.6)THEN
       OPEN(78,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
      ELSE
       OPEN(78,FILE=NAME,FORM='FORMATTED',STATUS='REPLACE')
       WRITE(78,88)'TIME','ETA','NCEP RAD  PRECIP'
       WRITE(78,89)'   3A  3F-A  6F-A 1A PK NE ETA ECH ECU S4 SCH SCU'
      ENDIF
88    FORMAT(A4,25X,A3,5X,16A)
89    FORMAT(20X,A49)

!=== Open runtime diagnositc file for ETA valid time record
      CALL OPENFILE(NAME,LDAS%ODIR,LDAS%EXPCODE,LDAS%EVTFILE)
      IF(LDAS%STARTCODE.EQ.1.OR.LDAS%STARTCODE.EQ.6)THEN
       OPEN(83,FILE=NAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',POSITION='APPEND')
      ELSE
       OPEN(83,FILE=NAME,STATUS='REPLACE',ACCESS='SEQUENTIAL')	  
      ENDIF

!=== Impose limit on time step
      IF (LDAS%TS.GT.3600) THEN
       PRINT *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       PRINT *,'ERROR, USER TIMESTEP > 3600 SECONDS!!!'
       PRINT *,'RESETTING LDAS%TS TO 3600!!!!'
       PRINT *,'######################################'
       LDAS%TS=3600
      ENDIF      
      IF(LDAS%TS.LT.1) then
       PRINT*,'Timestep can not be less than 1 minute, reset to 15 min.'
       WRITE(79,*)'Timestep can not be less than 1 minute, reset to 15 min.'
       LDAS%TS=15*60
      ENDIF

!=== Set Time
     LDAS%YR=LDAS%SYR
     LDAS%MO=LDAS%SMO
     LDAS%DA=LDAS%SDA
     LDAS%HR=LDAS%SHR
     LDAS%MN=LDAS%SMN
     LDAS%SS=LDAS%SSS
     LDAS%NUMOUTF=0               ! Initialize output time counters
     LDAS%NUMOUTM=0
     LDAS%NUMOUTC1=0
     LDAS%NUMOUTCT=0
     LDAS%NUMOUTC2=0
     CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT, &
&    LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS)

     PRINT*
     PRINT*,'********** GSFC-LDAS Driver **********'
     PRINT*,'Experiment Code: ','-',LDAS%EXPCODE,'-'
     PRINT*,'Starting TIME: ',LDAS%SMO,'/',LDAS%SDA,'/',LDAS%SYR
     PRINT*,'Ending TIME: ',LDAS%EMO,'/',LDAS%EDA,'/',LDAS%EYR
     PRINT* 

     WRITE(79,*)
     WRITE(79,*)'********** GSFC-LDAS Driver **********'
     WRITE(79,*)'Experiment Code: ','-',LDAS%EXPCODE,'-'
     WRITE(79,*)'Starting TIME: ',LDAS%SMO,'/',LDAS%SDA,'/',LDAS%SYR
     WRITE(79,*)'Ending TIME: ',LDAS%EMO,'/',LDAS%EDA,'/',LDAS%EYR
     WRITE(79,*) 

     PRINT*,'DOMAIN details:'
     WRITE(79,*)'DOMAIN details:'
      
!===> Read namelist of parameters depending on the domain
         read(unit=10,NML=nldas)
         PRINT*,'Running NLDAS ','-',LDAS%NC,LDAS%NR,'-'
         WRITE(79,*)'Running NLDAS ','-',LDAS%NC,LDAS%NR,'-'

     PRINT*,'Minimum Tile Area: ',LDAS%MINA
     PRINT*,'Maximum Tiles per grid: ',LDAS%MAXT
     WRITE(79,*)'Minimum Tile Area: ',LDAS%MINA
     WRITE(79,*)'Maximum Tiles per grid: ',LDAS%MAXT
     PRINT*      
     WRITE(79,*)
     
     PRINT*,'LSM details:'
     WRITE(79,*)'LSM details:'
     
!===> Read namelist of parameters depending on land surface model
         read(unit=10,NML=mos)
         PRINT*,'Running MOSAIC LSM:'
         PRINT*,'MOSAIC Active Restart File: ', LDAS%MOS_RFILE
         WRITE(79,*)'Running MOSAIC LSM:'
         WRITE(79,*)'MOSAIC Active Restart File: ', LDAS%MOS_RFILE

     PRINT*
     WRITE(79,*)
     
     PRINT*,'FORCING details:'
     WRITE(79,*)'FORCING details:'
     
!===> Read namelist of parameters depending on type of forcing
         read(unit=10,NML=eta)
         PRINT*,'Using ETA forcing:'
         WRITE(79,*)'Using ETA forcing:'

       IF (LDAS%KOSTER .GT. 0) THEN
        LDAS%NT = LDAS%NKTYPE
        LDAS%MOS_VFILE  = LDAS%MOS_KVFILE
        LDAS%MOS_MVFILE = LDAS%MOS_KMVFILE
        LDAS%MOS_SFILE  = LDAS%MOS_KSFILE
       ENDIF

!=== Setting Satellite LAI variables
       LDAS%LAITIME  = 0.0
       IF (LDAS%LAI .EQ. 2) PRINT*, "Using AVHRR Satellite LAI"
       IF (LDAS%LAI .EQ. 2) WRITE(79,*) "Using AVHRR Satellite LAI"
       IF (LDAS%LAI .EQ. 3) PRINT*, "Using MODIS Satellite LAI"
       IF (LDAS%LAI .EQ. 3) WRITE(79,*) "Using MODIS Satellite LAI"

!=== Set observed radiation times, ensures data is read during first timestep
       LDAS%AGRMTIME1 = 3000.0
       LDAS%AGRMTIME2 = 0.0

!=== Setting global observed precip times to zero to assure that data is read in during first time step
       LDAS%NRLTIME  = 0.0
       LDAS%PERSTIME = 0.0
       LDAS%HUFFTIME = 0.0
       LDAS%CMAPTIME = 0.0

!=== Marking which global observed precip data sources are being used
       IF (LDAS%GPCPSRC(1) .GT. 0) PRINT*, "Using NRL Observed precipitation"       
       IF (LDAS%GPCPSRC(2) .GT. 0) PRINT*, "Using PERSIANN Observed precipitation"
       IF (LDAS%GPCPSRC(3) .GT. 0) PRINT*, "Using HUFFMAN Observed precipitation"
       IF (LDAS%GPCPSRC(4) .GT. 0) PRINT*, "Using CMAP Observed precipitation"
       IF (LDAS%GPCPSRC(1) .GT. 0) WRITE(79,*) "Using NRL Observed precipitation"
       IF (LDAS%GPCPSRC(2) .GT. 0) WRITE(79,*) "Using PERSIANN Observed precipitation"
       IF (LDAS%GPCPSRC(3) .GT. 0) WRITE(79,*) "Using HUFFMAN Observed precipitation"
       IF (LDAS%GPCPSRC(4) .GT. 0) WRITE(79,*) "Using CMAP Observed precipitation"

!=== Determining the highest ranked global observed precip data source
       LDAS%MASTERPCP = 0
       DO I=1,3
        IF (LDAS%GPCPSRC(I) .EQ. 1) LDAS%HIERFLAG = I
       ENDDO

       IF (LDAS%DOMAIN .EQ. 1) THEN
         PRINT*,'Elevation Diff File: ',LDAS%ELEVFILE
         PRINT*,'NARR elev file: ',LDAS%NARRELEV
         PRINT*,'Awips Interp File A: ', LDAS%INTERPA
         PRINT*,'Awips Interp File B: ', LDAS%INTERPAF
         PRINT*,'1/2 Interp File A: ', LDAS%INTERPP
         PRINT*,'1/2 Interp File B: ', LDAS%INTERPPF
         WRITE(79,*)'Elevation Diff File: ',LDAS%ELEVFILE
         WRITE(79,*)'Awips Interp File A: ', LDAS%INTERPA
         WRITE(79,*)'Awips Interp File B: ', LDAS%INTERPAF
         WRITE(79,*)'1/2 Interp File A: ', LDAS%INTERPP
         WRITE(79,*)'1/2 Interp File B: ', LDAS%INTERPPF
       ELSE
	 PRINT*,'Elevation Diff File: ',LDAS%ELEVFILE
	 WRITE(79,*)'Elevation Diff File: ',LDAS%ELEVFILE
       ENDIF
       PRINT*
       WRITE(79,*)

!=== Initialize Statistics files conditional
      LDAS%FOROPEN=0
      LDAS%MOSOPEN=0
      LDAS%CLM1OPEN=0
      LDAS%CLM2OPEN=0
      LDAS%CATOPEN=0
      LDAS%NOAHOPEN=0
     
!=== Initialize Forcing Source Array
      DO I=1,16
       LDAS%FSOURCE(I)=0
      ENDDO

!=== Fix MAXT out of bounds problems
      IF(LDAS%MAXT.GT.LDAS%NT) LDAS%MAXT=LDAS%NT
      IF(LDAS%MAXT.LT.1      ) LDAS%MAXT=1
      
!=== Select which vegetation tile space and mask files

      IF (LDAS%VCLASS .EQ. 2) THEN 
        LDAS%VFILE  = LDAS%V2FILE
	LDAS%MFILE  = LDAS%M2FILE
	LDAS%FMFILE = LDAS%FM2FILE
      ENDIF
      
      PRINT*, 'MISCELLANEOUS details:'
      IF (LDAS%VCLASS .EQ. 1) PRINT*,'AVHRR UMD Vegetation', LDAS%VCLASS
      IF (LDAS%VCLASS .EQ. 2) PRINT*,'MODIS UMD Vegetation', LDAS%VCLASS      
      PRINT*,'Main Precipitation Source: ',LDAS%PRECSOR(1)
      PRINT*,'Precipitation File Regeneration: ',LDAS%REGENERATE
      IF (LDAS%PINKER.EQ.1) PRINT *,'Main Radiation Source: Pinker'
      IF (LDAS%NESDIS.EQ.1) PRINT *,'Main Radiation Source: NESDIS'
      IF ((LDAS%NESDIS.EQ.0).AND.(LDAS%PINKER.EQ.0)) &
&        PRINT *,'Main Radiation Source: ETA/EDAS'
      PRINT*,'TEMP,PRES,HUMID,LWRAD Elev. Adj.',LDAS%TEMPADJ, &
&        LDAS%PRESADJ,LDAS%HUMIDADJ,LDAS%LWRADADJ
      PRINT*,'MASK File: ', LDAS%MFILE
      PRINT*,'Vegetation File: ', LDAS%VFILE
      IF (LDAS%SOIL.EQ.1) PRINT *,'Original vegetation-based soil scheme'
      IF (LDAS%SOIL.EQ.2) PRINT *,'Reynolds soils'
!      PRINT*,'Soil File: ', LDAS%SFILE
      PRINT*

      WRITE(79,*) 'Other RUN details:'
      IF (LDAS%VCLASS .EQ. 1) WRITE(79,*) 'AVHRR UMD Vegetation', LDAS%VCLASS
      IF (LDAS%VCLASS .EQ. 2) WRITE(79,*) 'MODIS UMD Vegetation', LDAS%VCLASS
      WRITE(79,*) 'Main Precipitation Source: ',LDAS%PRECSOR(1)
      WRITE(79,*) 'Precipitation File Regeneration: ',LDAS%REGENERATE
      IF (LDAS%PINKER.EQ.1) WRITE(79,*) 'Main Radiation Source: Pinker'
      IF (LDAS%NESDIS.EQ.1) WRITE(79,*) 'Main Radiation Source: NESDIS'
      IF ((LDAS%NESDIS.EQ.0).AND.(LDAS%PINKER.EQ.0)) &
&        WRITE(79,*) 'Main Radiation Source: ETA/EDAS'
      WRITE(79,*) 'TEMP,PRES,HUMID,LWRAD Elev. Adj.',LDAS%TEMPADJ, &
&         LDAS%PRESADJ,LDAS%HUMIDADJ,LDAS%LWRADADJ
      WRITE(79,*) 'MASK File: ', LDAS%MFILE
      WRITE(79,*) 'Vegetation File: ', LDAS%VFILE
      IF (LDAS%SOIL.EQ.1) WRITE(79,*) 'Original vegtation-based soil scheme'
      IF (LDAS%SOIL.EQ.2) WRITE(79,*) 'Reynolds soils'
!      WRITE(79,*) 'Soil File: ', LDAS%SFILE
      WRITE(79,*)

      CLOSE(10)

      RETURN
      
      END

! Determine the last day of the month
!      SUBROUTINE LDOM(imon,iyr,iday)
!      integer imon,iyr,iday
!      integer dpm(12)  
!      data dpm/31,28,31,30,31,30,31,31,30,31,30,31/
!      iday=dpm(imon)
!      if(imon.eq.2) then
!       itest=mod(iyr,4)
!       if(itest.eq.0) iday=iday+1
!      endif
!      Return
!      End	

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
! forcefile.f: 
!
! DESCRIPTION:
!  This subroutine puts together a filename
!
! REVISION HISTORY:
!  4 Jan 2000: Paul Houser; Original Code
!=========================================================================
      SUBROUTINE OPENFILE(NAME,ODIR,EXPCODE,FFILE)
      IMPLICIT NONE

!=== Local Variables =====================================================
      CHARACTER*80 NAME,MKDIR
      CHARACTER*40 ODIR,FFILE
      INTEGER I,C,D,E,EXPCODE 
      CHARACTER*1 FBASE(80),FNAME(80),FCODE(80),FMKDIR(80)

!=== End Variable Definition =============================================

!=== Put together filename
 92   FORMAT(80A1)
 93   FORMAT(A80)
 94   FORMAT(A4,I3,A1)
 96   FORMAT(A40)
 100  FORMAT(A9)

      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      WRITE(90,96,REC=1) ODIR                       
      READ(90,92,REC=1) (FBASE(I),I=1,80)
      C=0
      DO I=1,80
       IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      WRITE(90,94,REC=1) '/EXP',EXPCODE,'/'                       
      READ(90,92,REC=1) (FCODE(I),I=1,80)
      D=0
      DO I=1,80
       IF(FCODE(I).EQ.(' ').AND.D.EQ.0)D=I-1
      ENDDO

      WRITE(90,96,REC=1) FFILE                       
      READ(90,92,REC=1) (FNAME(I),I=1,80)
      E=0
      DO I=1,80
       IF(FNAME(I).EQ.(' ').AND.E.EQ.0)E=I-1
      ENDDO

      WRITE(90,92,REC=1)(FBASE(I),I=1,C), &
&                       (FCODE(I),I=1,D), (FNAME(I),I=1,E)
      READ(90,93,REC=1)NAME

      WRITE(90,100,REC=1)'mkdir -p '
      READ(90,92,REC=1)(FMKDIR(I),I=1,9)

!=== Make the directories for the output files             
      WRITE(90,92,REC=1) (FMKDIR(I),I=1,9), &
&                        (FBASE(I),I=1,C),(FCODE(I),I=1,D)
      READ(90,93,REC=1)MKDIR
      CALL SYSTEM(MKDIR)
      CLOSE(90)

      RETURN
      END


