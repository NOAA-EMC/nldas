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
! getncep.f: 
!
! DESCRIPTION:
!  Opens, reads, and interpolates NCEP-LDAS forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!  The strategy for missing data is to backwards up to 10 days to get
!  forcing at the same time of day.
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  20 Dec 1999: Paul Houser; Allow for Eta Data Overwrite by NCEP data
!  27 Apr 2000: Brian Cosgrove; Turned zenith angle weighting back on.
!               Changed times supplied to ZTERP from GMT1 and GMT2 to 
!               LDAS%NCEPTIME1 and LDAS%NCEPTIME2
!   4 May 2000: Added 15 minutes to GMT1 and GMT2 to accurately
!               reflect the valid times of the NCEP radiation data
!  18 Aug 2000: Brian Cosgrove; Fixed error in date calculations so that
!               3600 second (1hr) timestep may be used.
!               Added code to make any calculated radiation forcing
!               values undefined if both ncep time 1 and ncep time 2
!               radiation values are undefined
!  27 Feb 2001: Brian Cosgrove; Added CZM into call for ZTERP subroutine
!  07 Mar 2001: Brian Cosgrove; Added code to allow for use of NASA-LDAS data 
!  04 Sep 2001: Brian Cosgrove; Changed tempgmt1,tempgmt2 to real to match
!               tick.f call, changed file name construction.
!  21 Aug 2002: Brian Cosgrove; Removed code that adjusted for 15 minute
!               offset in radiation fields supplied by NOAA or NASA
!               NLDAS realtime or retrospective standard forcing files.
!               This offset no longer exists as it is now dealt with 
!               during forcing file creation through zenith angle correction.
!               The need for this fix was just discovered...unfortunately
!               simulations before the date of this fix using 
!               this fortran subroutine have incorrectly shifted
!               radiation data.
!=========================================================================
      SUBROUTINE GETNCEP(LDAS,GRID)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (griddec) GRID(LDAS%NC,LDAS%NR)   

!=== Local Variables =====================================================
      INTEGER C,R,F,FERROR,TRY,ZDOY,BDOY,BYR,BMO,BDA,BHR,BMN
      INTEGER YR1,MO1,DA1,HR1,MN1,SS1,DOY1,TS1
      INTEGER YR2,MO2,DA2,HR2,MN2,SS2,DOY2,TS2,TEMPBSS,TEMPBTS
      INTEGER TEMPBDOY,TEMPBYR,TEMPBMO,TEMPBDA,TEMPBHR,TEMPBMN
      REAL*8 TIME1,TIME2,BTIME,NEWTIME1,NEWTIME2
      REAL TEMPGMT1,TEMPGMT2
      CHARACTER*80 NAME
      REAL WT1,WT2,ZW1,ZW2,CZB,CZE,CZM,GMT1,GMT2,SWT1,SWT2

      INTEGER FINDTIME1    ! 0=don't get new file for 1st time  (1=get new file)
	INTEGER FINDTIME2    ! 1=Get a new file for 2nd time  (0=don't get new file)
	INTEGER MOVETIME     ! 1=Move Time 2 data into Time 1  

!=== End Variable Definition =============================================
      TRY=-999

!====Assumption will be not to find or move any data
      FINDTIME1=0
	FINDTIME2=0
	MOVETIME=0

!=== Determine Required NCEP Data Times (The previous hour and the future hour)
	IF (((LDAS%TS.EQ.3600).AND.(LDAS%MN.NE.0))
     &  .OR.((LDAS%TS.NE.3600).AND.(LDAS%TSCOUNT.GT.1))) THEN
      YR1=LDAS%YR    !Previous Hour
      MO1=LDAS%MO
      DA1=LDAS%DA
      HR1=LDAS%HR
      MN1=0
      SS1=0
      TS1=0
      CALL TICK(TIME1,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)
      YR2=LDAS%YR    !Next Hour
      MO2=LDAS%MO
      DA2=LDAS%DA
      HR2=LDAS%HR
      MN2=0
      SS2=0
      TS2=60*60
      CALL TICK(TIME2,DOY2,GMT2,YR2,MO2,DA2,HR2,MN2,SS2,TS2)
	ELSEIF ((LDAS%TS.NE.3600).AND.(LDAS%TSCOUNT.EQ.1) 
     &  .AND.(LDAS%MN.NE.0)) THEN
c	print *,'finding normal times'
      YR1=LDAS%YR    !Previous Hour
      MO1=LDAS%MO
      DA1=LDAS%DA
      HR1=LDAS%HR
      MN1=0
      SS1=0
      TS1=0
      CALL TICK(TIME1,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)
      YR2=LDAS%YR    !Next Hour
      MO2=LDAS%MO
      DA2=LDAS%DA
      HR2=LDAS%HR
      MN2=0
      SS2=0
      TS2=60*60
      CALL TICK(TIME2,DOY2,GMT2,YR2,MO2,DA2,HR2,MN2,SS2,TS2)

	ELSE
      YR1=LDAS%YR    !Previous Hour
      MO1=LDAS%MO
      DA1=LDAS%DA
      HR1=LDAS%HR
      MN1=0
      SS1=0
      TS1=-3600
      CALL TICK(TIME1,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)
      YR2=LDAS%YR    !Next Hour
      MO2=LDAS%MO
      DA2=LDAS%DA
      HR2=LDAS%HR
      MN2=0
      SS2=0
      TS2=0
      CALL TICK(TIME2,DOY2,GMT2,YR2,MO2,DA2,HR2,MN2,SS2,TS2)
        ENDIF


	
C
c      IF((TIME1.GT.LDAS%NCEPTIME2).OR.
c     & (TIME1.EQ.LDAS%NCEPTIME2.AND.LDAS%MN.GT.0)) then  
        IF (((LDAS%TS.EQ.3600).AND.(LDAS%MN.NE.0))
     &  .OR.((LDAS%TS.NE.3600).AND.(LDAS%TSCOUNT.GT.1))) THEN
       IF (TIME2.GT.LDAS%NCEPTIME2.AND.LDAS%MN.GT.0) then  
            MOVETIME=1
	    FINDTIME2=1
      ENDIF
        ELSEIF ((LDAS%TS.NE.3600).AND.(LDAS%TSCOUNT.EQ.1)
     &  .AND.(LDAS%MN.NE.0)) THEN
       IF (TIME2.GT.LDAS%NCEPTIME2.AND.LDAS%MN.GT.0) then
            MOVETIME=1
            FINDTIME2=1
      ENDIF

	endif

        if (LDAS%TS.eq.3600) then
       IF (TIME2.GT.LDAS%NCEPTIME2) then
            MOVETIME=1
            FINDTIME2=1
      ENDIF
        endif

	IF(LDAS%TSCOUNT.EQ.1) then    !beginning of the run	
c	print *,'at beginning of run'
	      FINDTIME1=1
		FINDTIME2=1
	      MOVETIME=0
      ENDIF
	
	

      IF(MOVETIME.eq.1) then
c	print *,'moving time2 into time1'
        LDAS%NCEPTIME1=LDAS%NCEPTIME2
        LDAS%SKIPINTP1=LDAS%SKIPINTP2
        DO F=1,LDAS%NF
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
           GRID(C,R)%NCEPDATA1(F)=GRID(C,R)%NCEPDATA2(F)
          ENDDO
         ENDDO
        ENDDO
	 ENDIF    !end of movetime=1
	 
	 IF(FINDTIME1.eq.1) then
!=== The following looks back 10 days, at the same hour to fill data gaps.
        LDAS%SKIPINTP1=0
        FERROR=0
        TRY=0  
        TS1=-60*60*24
        DO WHILE (FERROR.EQ.0)
         TRY=TRY+1
      IF (LDAS%FNCEP.EQ.1) THEN
         CALL NCEPFILE(NAME,LDAS%NCEPDIR,YR1,MO1,DA1,HR1,LDAS)
      ENDIF
      IF (LDAS%FNASA.EQ.1) THEN
         CALL NCEPFILE(NAME,LDAS%NASADIR,YR1,MO1,DA1,HR1,LDAS)
      ENDIF

c	print *,'getting 1',NAME
         CALL RETNCEP(1,LDAS,GRID,NAME,FERROR,1,0)
         IF(FERROR.EQ.1)LDAS%NCEPTIME1=TIME1
         CALL TICK(TIME1,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)
         IF(TRY.GT.11)THEN
          WRITE(*,*)'ERROR: NCEP Data gap exceeds 10 days on file 1'
          WRITE(79,*)'ERROR: NCEP Data gap exceeds 10 days on file 1'
          STOP
         ENDIF
         IF(LDAS%FETA.EQ.1)THEN
          IF(FERROR.eq.0)LDAS%SKIPINTP1=1
          FERROR=1 !Only use current NCEP data if ETA is read
         ENDIF      
        ENDDO
!=== End of Data Search
      ENDIF   !end of findtime=1	   	


      IF(FINDTIME2.eq.1) then  
!=== The following looks back 10 days, at the same hour to fill data gaps.
        LDAS%SKIPINTP2=0
        FERROR=0
        TRY=0  
        TS2=-60*60*24
        DO WHILE (FERROR.EQ.0)
         TRY=TRY+1
      IF (LDAS%FNCEP.EQ.1) THEN
         CALL NCEPFILE(NAME,LDAS%NCEPDIR,YR2,MO2,DA2,HR2,LDAS)
      ENDIF
      IF (LDAS%FNASA.EQ.1) THEN
         CALL NCEPFILE(NAME,LDAS%NASADIR,YR2,MO2,DA2,HR2,LDAS)
      ENDIF
c        print *,'getting 2',NAME

         CALL RETNCEP(2,LDAS,GRID,NAME,FERROR,1,0)
         IF(FERROR.EQ.1)LDAS%NCEPTIME2=TIME2
         CALL TICK(TIME2,DOY2,GMT2,YR2,MO2,DA2,HR2,MN2,SS2,TS2)
         IF(TRY.GT.11)THEN
          WRITE(*,*)'ERROR: NCEP Data gap exceeds 10 days on file 2'
          WRITE(79,*)'ERROR: NCEP Data gap exceeds 10 days on file 2'
          STOP
         ENDIF
         IF(LDAS%FETA.EQ.1)THEN
          IF(FERROR.eq.0)LDAS%SKIPINTP2=1
          FERROR=1 !Only use current NCEP data if ETA is read
         ENDIF           
        ENDDO
!=== End of Data Search

      ENDIF   ! End of findtime2=1

!=== Produce Reporting for fsource.dat file
      IF(LDAS%SKIPINTP1.eq.1.or.LDAS%SKIPINTP2.eq.1)THEN
       LDAS%FSOURCE(6)=0            !Forcing not read
      ELSEIF(TRY.eq.(1))THEN
       LDAS%FSOURCE(6)=1            !Current forcing read 
      ELSEIF(TRY.eq.(-999))THEN
       LDAS%FSOURCE(6)=LDAS%FSOURCE(6)   !Forcing in Memory
      ELSE
       LDAS%FSOURCE(6)=TRY*(-1)+1   !Old forcing read
      ENDIF

        BTIME=LDAS%NCEPTIME1
        CALL TIME2DATE(BTIME,BDOY,GMT1,BYR,BMO,BDA,BHR,BMN)
c	GMT1=GMT1+.25
C	Make new decimal time that reflects added 15 minutes

C       This was commented out so that it no longer shifts by 15 mins
C       since this shift is no longer necessitated by forcing files
C       Forcing files already contain this shift
cSUBROUTINE TICK(TIME,DOY,GMT,YR,MO,DA,HR,MN,SS,TS)
	TEMPBDOY=BDOY
c	TEMPGMT1=GMT1-.25
	TEMPGMT1=GMT1
	TEMPBYR=BYR
	TEMPBMO=BMO
	TEMPBDA=BDA
	TEMPBHR=BHR
	IF (TEMPBHR.EQ.24) TEMPBHR=0
	TEMPBMN=BMN
	TEMPBSS=0
	TEMPBTS=900
	CALL TICK(NEWTIME1,TEMPBDOY,TEMPGMT1,
     &  TEMPBYR,TEMPBMO,TEMPBDA,TEMPBHR,TEMPBMN,
     &  TEMPBSS,TEMPBTS)

        BTIME=LDAS%NCEPTIME2
        CALL TIME2DATE(BTIME,BDOY,GMT2,BYR,BMO,BDA,BHR,BMN)
c	GMT2=GMT2+.25
C       Make new decimal time that reflects added 15 minutes
        TEMPBDOY=BDOY
c        TEMPGMT2=GMT2-.25
        TEMPGMT2=GMT2
        TEMPBYR=BYR
        TEMPBMO=BMO
        TEMPBDA=BDA
        TEMPBHR=BHR
        IF (TEMPBHR.EQ.24) TEMPBHR=0
        TEMPBMN=BMN
        TEMPBSS=0
        TEMPBTS=900
        CALL TICK(NEWTIME2,TEMPBDOY,TEMPGMT2,
     &  TEMPBYR,TEMPBMO,TEMPBDA,TEMPBHR,TEMPBMN,
     &  TEMPBSS,TEMPBTS)


!=== Interpolate Data in time      
      IF(LDAS%SKIPINTP1.ne.1.or.LDAS%SKIPINTP2.ne.1)THEN !Skip interp if getETA routine on, and no NCEP data


       WT1=(LDAS%NCEPTIME2-LDAS%TIME)/(LDAS%NCEPTIME2-LDAS%NCEPTIME1)
       WT2=1.0-WT1
       SWT1=(NEWTIME2-LDAS%TIME)/(NEWTIME2-NEWTIME1)
       SWT2=1.0-SWT1

       DO F=1,LDAS%NF
        IF(F.EQ.3)THEN
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
	   ZDOY=LDAS%DOY

C	Compute and apply zenith angle weights

           CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1       GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,CZM,LDAS,GRID)
           GRID(C,R)%FORCING(F)=GRID(C,R)%NCEPDATA1(F)*ZW1+
     1		                GRID(C,R)%NCEPDATA2(F)*ZW2

c	In cases of small cos(zenith) angles, use linear weighting
C       to avoid overly large weights

           IF((GRID(C,R)%FORCING(F).GT.GRID(C,R)%NCEPDATA1(F).AND.
     1         GRID(C,R)%FORCING(F).GT.GRID(C,R)%NCEPDATA2(F)).AND.
     2         (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
            GRID(C,R)%FORCING(F)=GRID(C,R)%NCEPDATA1(F)*SWT1+
     1                           GRID(C,R)%NCEPDATA2(F)*SWT2
c	print *,'using swt1 and swt2',swt1,swt2
           ENDIF

                 if (GRID(C,R)%FORCING(F).gt.1367) then
                  print *,'warning, SW RADIATION TOO HIGH!!'
                  print *,'it is',GRID(C,R)%FORCING(F)
                  print *,'NCEPDATA1=',GRID(C,R)%NCEPDATA1(F)
                  print *,'NCEPDATA2=',GRID(C,R)%NCEPDATA2(F)
                  print *,'ZW1=',ZW1,'ZW2=',ZW2
		  print *,'SWT1=',SWT1,'SWT2=',SWT2
                  GRID(C,R)%FORCING(F)=GRID(C,R)%NCEPDATA1(F)*SWT1+
     1                           GRID(C,R)%NCEPDATA2(F)*SWT2
                 endif

	IF ( (GRID(C,R)%NCEPDATA1(F).eq.LDAS%UDEF).AND.
     &  (GRID(C,R)%NCEPDATA2(F).eq.LDAS%UDEF) ) then
	GRID(C,R)%FORCING(F)=LDAS%UDEF
	ENDIF
          ENDDO
         ENDDO  	 
        ELSE IF(F.EQ.8.OR.F.EQ.9)THEN	!Do block Precipitation Interpolation
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
           GRID(C,R)%FORCING(F)=GRID(C,R)%NCEPDATA2(F)	
          ENDDO
         ENDDO
        ELSE !Linearly interpolate everything else
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
           GRID(C,R)%FORCING(F)=GRID(C,R)%NCEPDATA1(F)*WT1+
     1                          GRID(C,R)%NCEPDATA2(F)*WT2
          ENDDO
         ENDDO            
        ENDIF    
       ENDDO   

!=== ADJUST PRECIP TO VALUE PER SEC DATA SHOULD COME IN AS VALUE PER 
!=== TIME PERIOD ONE HOUR FOR NCEP, THREE HOUR FOR ETA

       DO C=1,LDAS%NC
        DO R=1,LDAS%NR
         GRID(C,R)%FORCING(8)=GRID(C,R)%FORCING(8)/(60.0*60.0)
         GRID(C,R)%FORCING(9)=GRID(C,R)%FORCING(9)/(60.0*60.0)
        ENDDO
       ENDDO 

      ENDIF !LDAS%SKIPINTP



      RETURN
      END





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
! ncepfile.f: 
!
! DESCRIPTION:
!  This subroutine puts together the ncep data filename
!
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
!=========================================================================
      SUBROUTINE NCEPFILE(NAME,NCEPDIR,YR,MO,DA,HR,LDAS)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      type (ldasdec) LDAS


!=== Local Variables =====================================================
      CHARACTER*80 NAME
      CHARACTER*40 NCEPDIR
      INTEGER YR,MO,DA,HR,I,C 

      CHARACTER*1  FNAME(80),FBASE(80),FSUBS(80)
      CHARACTER*1  FTIME(10),FDIR(15)

!=== End Variable Definition =============================================

!=== Put together filename
 91   FORMAT(A4,I3,A11,I3)
 92   FORMAT(80A1)
 93   FORMAT(A80)
 94   FORMAT(I4,I2,I2,I2)
 95   FORMAT(10A1)
 96   FORMAT(A40)
 97   FORMAT(A12) 
 89   FORMAT(A14)
 98   FORMAT(A1,I4,A1,I4,I2,I2,A1)
 99   FORMAT(15A1)
      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      WRITE(90,98,REC=1)'/',YR,'/',YR,MO,DA,'/'
      READ(90,99,REC=1)FDIR
      DO I=1,15
       IF(FDIR(I).EQ.(' '))FDIR(I)='0'
      ENDDO

      WRITE(90,94,REC=1)YR,MO,DA,HR
      READ(90,95,REC=1)FTIME
      DO I=1,10
       IF(FTIME(I).EQ.(' '))FTIME(I)='0'
      ENDDO

      IF (LDAS%FNCEP.EQ.1) THEN
       WRITE(90,89,REC=1)'.lsmforce_anal'
       READ(90,92,REC=1) (FSUBS(I),I=1,14)
      ENDIF

      IF (LDAS%FNASA.EQ.1) THEN
       WRITE(90,97,REC=1)'.FORCING.GRB'
       READ(90,92,REC=1) (FSUBS(I),I=1,12)
      ENDIF

      WRITE(90,96,REC=1) NCEPDIR                       
      READ(90,92,REC=1) (FBASE(I),I=1,80)

      C=0
      DO I=1,80
       IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO

      IF (LDAS%FNCEP.EQ.1) THEN
        WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FDIR(I),I=1,15),
     1                  (FTIME(I),I=1,10),(FSUBS(I),I=1,14 )
      ENDIF 

      IF (LDAS%FNASA.EQ.1) THEN
        WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FDIR(I),I=1,15),
     1                  (FTIME(I),I=1,10),(FSUBS(I),I=1,12 ) 
	ENDIF
      READ(90,93,REC=1)NAME
      CLOSE(90)

      RETURN
      END












































