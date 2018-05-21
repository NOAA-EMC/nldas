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
! getpar.f: 
!
! DESCRIPTION:
!  Opens, reads, interpolates and overlays radiation forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!
! REVISION HISTORY:
!  19  Dec 2000: Brian Cosgrove; Initial code based on getrad.f
!  05  Feb 2002: Brian Cosgrove; Added year, month into call for retpar
!=========================================================================
      SUBROUTINE GETPAR(LDAS,GRID)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (griddec) GRID(LDAS%NC,LDAS%NR)   

!=== Local Variables =====================================================
      INTEGER C,R,F,TRY,ZDOY,AA,BB,CC,LOOP
      INTEGER YR1,MO1,DA1,HR1,MN1,SS1,DOY1,TS1
      INTEGER YR2,MO2,DA2,HR2,MN2,SS2,DOY2,TS2
      REAL*8 TIME1,TIME2
      CHARACTER*80 NAME1
      REAL WT1,WT2,ZW1,ZW2,CZB,CZE,CZM,GMT1,GMT2,czmean

      
!=== End Variable Definition =============================================

!=== Determine Required Observed Radiation Data Times 
!=== (The previous hour and the future hour)
      YR1=LDAS%YR    !Previous Hour
      MO1=LDAS%MO
      DA1=LDAS%DA
      IF ((LDAS%MN.GE.0).AND.(LDAS%MN.LT.15)) HR1=LDAS%HR-1
      IF (LDAS%MN.GE.15) HR1=LDAS%HR
      MN1=15
      SS1=0
      TS1=0
      CALL TICK(TIME1,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)
      YR2=LDAS%YR    !Next Hour
      MO2=LDAS%MO
      DA2=LDAS%DA
      IF ((LDAS%MN.GE.0).AND.(LDAS%MN.LT.15)) HR2=LDAS%HR-1
      IF (LDAS%MN.GE.15) HR2=LDAS%HR
      MN2=15
      SS2=0
      TS2=60*60
      CALL TICK(TIME2,DOY2,GMT2,YR2,MO2,DA2,HR2,MN2,SS2,TS2)

!=== Check to see if required PAR data in memory, if not, read it
       IF(TIME1.NE.LDAS%PARTIME1)THEN !We need to get new TIME1 data
        IF(TIME1.EQ.LDAS%PARTIME2)THEN !Transfer TIME2 data TIME1
         LDAS%PARTIME1=LDAS%PARTIME2
	 LDAS%PARSTAT1=LDAS%PARSTAT2
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
           GRID(C,R)%PARDATA1=GRID(C,R)%PARDATA2
          ENDDO
         ENDDO
        ELSE  !Get RAD Data for PAR TIME1
!=== The following retreives observed radiation data for Time1
        if ((yr1.ge.1996).and.(yr1.le.2000)) then
         CALL RADFILEPAR(NAME1,LDAS,YR1,MO1,DA1,HR1,1)
         CALL RETPAR(1,LDAS,GRID,NAME1,LDAS%PARSTAT1,1,yr1,mo1)
	else
         CALL RADFILEPAROLD(NAME1,LDAS,YR1,MO1,DA1,HR1,1)
         CALL RETPAR(1,LDAS,GRID,NAME1,LDAS%PARSTAT1,1,yr1,mo1)
         ENDIF
         IF (LDAS%PARSTAT1.NE.0) LDAS%PARTIME1=TIME1
!=== End of Data Search For PAR Time1
        ENDIF
       ENDIF
       IF(TIME2.NE.LDAS%PARTIME2)THEN !Get RAD Data for TIME2
!=== The following retreives observed radiation data for TIME2
        if ((yr2.ge.1996).and.(yr2.le.2000)) then
        CALL RADFILEPAR(NAME1,LDAS,YR2,MO2,DA2,HR2,1)
        CALL RETPAR(2,LDAS,GRID,NAME1,LDAS%PARSTAT2,1,yr2,mo2)
	else
        CALL RADFILEPAROLD(NAME1,LDAS,YR2,MO2,DA2,HR2,1)
        CALL RETPAR(2,LDAS,GRID,NAME1,LDAS%PARSTAT2,1,yr2,mo2)
	ENDIF
        IF(LDAS%PARSTAT2.NE.0)LDAS%PARTIME2=TIME2
!=== End of Data Search For PAR Time2
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%PARDATA2.NE.LDAS%UDEF).AND.
     &      ((GRID(C,R)%PARDATA2.GT.-998.999).AND.
     &      (GRID(C,R)%PARDATA2.LT.-999.01))) THEN
             IF (GRID(C,R)%PARDATA2.LT.0.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'PAR OUT OF RANGE (',
     &          GRID(C,R)%PARDATA2,') at (',c,',',r,
     &              ') CORRECTING TO 0.0'
              GRID(C,R)%PARDATA2=0.0
             ENDIF

             IF (GRID(C,R)%PARDATA2.GT.816) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'PAR OUT OF RANGE (',
     &          GRID(C,R)%PARDATA2,') at (',c,',',r,
     &          ') CORRECTING TO 816'
              GRID(C,R)%PARDATA2=816
             ENDIF
            ENDIF
	ENDDO
	ENDDO

       ENDIF

	if(LDAS%TSCOUNT.EQ.1) then
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%PARDATA1.NE.LDAS%UDEF).AND.
     &      ((GRID(C,R)%PARDATA1.GT.-998.999).AND.
     &      (GRID(C,R)%PARDATA1.LT.-999.01))) THEN
             IF (GRID(C,R)%PARDATA1.LT.0.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'PAR OUT OF RANGE (',
     &          GRID(C,R)%PARDATA1,') at (',c,',',r,
     &              ') CORRECTING TO 0.0'
              GRID(C,R)%PARDATA1=0.0
             ENDIF

             IF (GRID(C,R)%PARDATA1.GT.816) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'PAR OUT OF RANGE (',
     &          GRID(C,R)%PARDATA1,') at (',c,',',r,
     &          ') CORRECTING TO 816'
              GRID(C,R)%PARDATA1=816
             ENDIF
            ENDIF
        ENDDO
        ENDDO
	endif



!=== Print out Status of data holdings

	IF (LDAS%PARSTAT1.EQ.0) WRITE(79,*)'PAR UNAVAILABLE',
     &  MO1,DA1,YR1,HR1
        IF (LDAS%PARSTAT1.EQ.1) WRITE(79,*)'PAR AVAILABLE',
     &  MO1,DA1,YR1,HR1
        IF (LDAS%PARSTAT1.EQ.2) WRITE(79,*)'PAR UNDEFINED',
     &  MO1,DA1,YR1,HR1

        IF (LDAS%PARSTAT2.EQ.0) WRITE(79,*)'PAR UNAVAILABLE',
     &  MO2,DA2,YR2,HR2
        IF (LDAS%PARSTAT2.EQ.1) WRITE(79,*)'PAR AVAILABLE',
     &  MO2,DA2,YR2,HR2
        IF (LDAS%PARSTAT2.EQ.2) WRITE(79,*)'PAR UNDEFINED',
     &  MO2,DA2,YR2,HR2


!=== If either NESDIS or PAR data is available, set
!=== up looping structure for forcing replacement process

          AA=1
          BB=1
	  CC=1

!== Loop through and replace data as possible with NESDIS/PAR data
!== This depends on options specified in ldas.crd as well as actual
!== data holdings.  The primary data set is looped through second, and
!== is placed on top of the secondary data set as well as the original
!== NCEP radiation data.

	 DO LOOP=AA,BB,CC
	  IF (LOOP.EQ.1) THEN
	   IF ((LDAS%PARSTAT1.EQ.1).AND.(LDAS%PARSTAT2.EQ.1)) THEN
!== Compute weights and zenith angle information.  Replace forcing 
!== with PAR data.
            WT1=(LDAS%PARTIME2-LDAS%TIME)/(LDAS%PARTIME2-
     &       LDAS%PARTIME1)
            WT2=1.0-WT1
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
              ZDOY=LDAS%DOY
              CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1         GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2         CZM,LDAS,GRID,czmean)
              GRID(C,R)%PAR=LDAS%UDEF
              IF ((GRID(C,R)%PARDATA1.GT.0.0).AND.
     &          (GRID(C,R)%PARDATA2.GT.0.0)) THEN
                GRID(C,R)%PAR=GRID(C,R)%PARDATA1*ZW1+
     1          GRID(C,R)%PARDATA2*ZW2
               IF((GRID(C,R)%PAR.GT.GRID(C,R)%PARDATA1.AND.
     1          GRID(C,R)%PAR.GT.GRID(C,R)%PARDATA2).AND.
     2          (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%PAR=GRID(C,R)%PARDATA1*WT1+
     1           GRID(C,R)%PARDATA2*WT2
               ENDIF
              ENDIF
              IF ((GRID(C,R)%PARDATA1.GT.0.0).AND.
     &          (GRID(C,R)%PARDATA2.LE.0.0)) THEN
              IF (CZB.GT.0.0) THEN
                GRID(C,R)%PAR=GRID(C,R)%PARDATA1*CZM/CZB
              ELSE
                GRID(C,R)%PAR=GRID(C,R)%PARDATA1*0.0
              ENDIF
               IF((GRID(C,R)%PAR.GT.GRID(C,R)%PARDATA1.AND.
     1             GRID(C,R)%PAR.GT.0.0).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%PAR=GRID(C,R)%PARDATA1*WT1+
     1                             0.0*WT2
               ENDIF
              ENDIF
              IF ((GRID(C,R)%PARDATA1.LE.0.0).AND.
     &          (GRID(C,R)%PARDATA2.GT.0.0)) THEN
              IF (CZE.GT.0.0) THEN 
                GRID(C,R)%PAR= GRID(C,R)%PARDATA2*CZM/CZE
              ELSE
                GRID(C,R)%PAR= GRID(C,R)%PARDATA2*0.0
              ENDIF 
               IF((GRID(C,R)%PAR.GT.0.0.AND.
     1             GRID(C,R)%PAR.GT.GRID(C,R)%PARDATA2).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%PAR=0.0*WT1+
     1                             GRID(C,R)%PARDATA2*WT2
               ENDIF
              ENDIF
              if (GRID(C,R)%PAR.gt.1367) then
              print *,'warning, PAR SW RADIATION TOO HIGH'
               print *,'it is',GRID(C,R)%PAR
	       print *,'PAR1=',GRID(C,R)%PARDATA1
	       print *,'PAR2=',GRID(C,R)%PARDATA2
	       print *,'wt1,wt2,czb,cze,czm,zw1,zw2'
	       print *,wt1,wt2,czb,cze,czm,zw1,zw2
               stop
              endif
             ENDDO
            ENDDO
           ENDIF

          IF ((LDAS%PARSTAT1.EQ.1).AND.(LDAS%PARSTAT2.NE.1)) THEN
!== Compute weights and zenith angle information.  Replace forcing 
!== with zenith extrapolated PAR data
           WT1=(LDAS%PARTIME2-LDAS%TIME)/(LDAS%PARTIME2-
     &     LDAS%PARTIME1)
           WT2=1.0-WT1
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             ZDOY=LDAS%DOY
             CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1             GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2             CZM,LDAS,GRID,czmean)
             GRID(C,R)%PAR=LDAS%UDEF
            IF (GRID(C,R)%PARDATA1.GT.0.0) THEN
            IF (CZB.GT.0.0) THEN
             GRID(C,R)%PAR=GRID(C,R)%PARDATA1*CZM/CZB
            ELSE
             GRID(C,R)%PAR=GRID(C,R)%PARDATA1*0.0
            ENDIF 

              IF((GRID(C,R)%PAR.GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
               GRID(C,R)%PAR=GRID(C,R)%PARDATA1*WT1
              ENDIF

              IF (GRID(C,R)%PAR.gt.1367) then
              print *,'warning, OBSERVED PAR RADIATION TOO HIGH'
               print *,'it is',GRID(C,R)%PAR
               print *,'PAR1=',GRID(C,R)%PARDATA1
               print *,'PAR2=',GRID(C,R)%PARDATA2
               print *,'wt1,wt2,czb,cze,czm,zw1,zw2'
               print *,wt1,wt2,czb,cze,czm,zw1,zw2
               stop
              ENDIF
             ENDIF
            ENDDO
           ENDDO
          ENDIF


          IF ((LDAS%PARSTAT1.NE.1).AND.(LDAS%PARSTAT2.EQ.1)) THEN
!== Compute weights and zenith angle information.  Replace forcing 
!== with zenith extrapolated PAR data
           WT1=(LDAS%PARTIME2-LDAS%TIME)/(LDAS%PARTIME2-
     &     LDAS%PARTIME1)
           WT2=1.0-WT1
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             ZDOY=LDAS%DOY
             CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1       GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2       CZM,LDAS,GRID,czmean)
             GRID(C,R)%PAR=LDAS%UDEF
             IF (GRID(C,R)%PARDATA2.GT.0.0) THEN
             IF (CZE.GT.0.0) THEN 
              GRID(C,R)%PAR=GRID(C,R)%PARDATA2*CZM/CZE
	     ELSE
   	      GRID(C,R)%PAR=GRID(C,R)%PARDATA2*0.0
	     ENDIF
              IF((GRID(C,R)%PAR.GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2        (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
               GRID(C,R)%PAR=0.0*WT1+
     1         GRID(C,R)%PARDATA2*WT2
              ENDIF
              IF (GRID(C,R)%PAR.gt.1367) then
              print *,'warning, OBSERVED PAR RADIATION TOO HIGH'
               print *,'it is',GRID(C,R)%PAR
               print *,'PAR1=',GRID(C,R)%PARDATA1
               print *,'PAR2=',GRID(C,R)%PARDATA2
               print *,'wt1,wt2,czb,cze,czm,zw1,zw2,czm/cze'
               print *,wt1,wt2,czb,cze,czm,zw1,zw2,czm/cze
               stop
              ENDIF
             ENDIF
            ENDDO
           ENDDO
          ENDIF

          IF ((LDAS%PARSTAT1.NE.1).AND.(LDAS%PARSTAT2.NE.1)) THEN
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             GRID(C,R)%PAR=LDAS%UDEF
            ENDDO
           ENDDO
          ENDIF
	 ENDIF   ! END 'IF LOOP.EQ.1' loop

        ENDDO
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
!   GSFC - RAD - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! radfilepar.f: 
!
! DESCRIPTION:
!  This subroutine puts together the radiation data filenames
!
!
! REVISION HISTORY:
!  28  Oct 1999: Brian Cosgrove; Initial code
!=========================================================================
      SUBROUTINE RADFILEPAR(NAME1,LDAS,YR,MO,DA,HR,FLAG)
      USE ldas_module         ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      TYPE (ldasdec) LDAS

!=== Local Variables =====================================================
      CHARACTER*80 NAME1
      INTEGER YR,LOCALYR,MO,DA,HR,I,C,FLAG

      CHARACTER*1  FNAME1(80),FBASE1(80),FSUBS1(80)
      CHARACTER*1  FTIME1(12),FDIR1(8)

      CHARACTER*1  FNAME2(80),FBASE2(80),FSUBS2(80)
      CHARACTER*1  FTIME2(12),FDIR2(8)



!=== End Variable Definition =============================================

!=== Put together filename

	IF (YR.LE.1999) LOCALYR=YR-1900
        IF (YR.GE.2000) LOCALYR=YR-2000

 92   FORMAT(80A1)
 93   FORMAT(A80)
 94   FORMAT(I4,I2,I2,A1,I2)
 95   FORMAT(12A1)
 96   FORMAT(A40)
 97   FORMAT(A6) 
 98   FORMAT(A1,I4,I2,A1)
 99   FORMAT(8A1)

	IF (FLAG.EQ.1) THEN
C	Generate Filename for PAR Radiation
      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      WRITE(90,98,REC=1)'/',YR,MO,'/'
      READ(90,99,REC=1)FDIR1
      DO I=1,8
       IF(FDIR1(I).EQ.(' '))FDIR1(I)='0'
      ENDDO

      WRITE(90,94,REC=1)YR,MO,DA,'.',HR
      READ(90,95,REC=1)FTIME1
      DO I=1,11
       IF(FTIME1(I).EQ.(' '))FTIME1(I)='0'
      ENDDO

      WRITE(90,97,REC=1)'.par.i'
      READ(90,92,REC=1) (FSUBS1(I),I=1,6)

      WRITE(90,96,REC=1) LDAS%PARDIR
      READ(90,92,REC=1) (FBASE1(I),I=1,80)

      C=0
      DO I=1,80
       IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO
      WRITE(90,92,REC=1)(FBASE1(I),I=1,C), (FDIR1(I),I=1,8),
     1                  (FTIME1(I),I=1,11),(FSUBS1(I),I=1,6)
      READ(90,93,REC=1)NAME1
      CLOSE(90)
C       Generate Filename for NESDIS HOURLY Radiation
C       Generate Filename for NESDIS 30 Hour Rotating Radiation

	ENDIF

      RETURN
      END

      SUBROUTINE RADFILEPAROLD(NAME1,LDAS,YR,MO,DA,HR,FLAG)
      USE ldas_module         ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      TYPE (ldasdec) LDAS
 
!=== Local Variables =====================================================
      CHARACTER*80 NAME1
      INTEGER YR,LOCALYR,MO,DA,HR,I,C,FLAG
 
      CHARACTER*1  FNAME1(80),FBASE1(80),FSUBS1(80)
      CHARACTER*1  FTIME1(12),FDIR1(8)
 
      CHARACTER*1  FNAME2(80),FBASE2(80),FSUBS2(80)
      CHARACTER*1  FTIME2(12),FDIR2(8)
 
 
 
!=== End Variable Definition =============================================
 
!=== Put together filename
 
        IF (YR.LE.1999) LOCALYR=YR-1900
        IF (YR.GE.2000) LOCALYR=YR-2000
 
 92   FORMAT(80A1)
 93   FORMAT(A80)
 94   FORMAT(I4,I2,I2,A1,I2)
 95   FORMAT(12A1)
 96   FORMAT(A40)
 97   FORMAT(A6)
 98   FORMAT(A1,I4,I2,A1)
 99   FORMAT(8A1)
 
        IF (FLAG.EQ.1) THEN
C       Generate Filename for PAR Radiation
      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
 
      WRITE(90,98,REC=1)'/',YR,MO,'/'
      READ(90,99,REC=1)FDIR1
      DO I=1,8
       IF(FDIR1(I).EQ.(' '))FDIR1(I)='0'
      ENDDO
 
      WRITE(90,94,REC=1)YR,MO,DA,'.',HR
      READ(90,95,REC=1)FTIME1
      DO I=1,11
       IF(FTIME1(I).EQ.(' '))FTIME1(I)='0'
      ENDDO
 
      WRITE(90,97,REC=1)'.par.i'
      READ(90,92,REC=1) (FSUBS1(I),I=1,6)
 
      WRITE(90,96,REC=1) LDAS%PARDIROLD
      READ(90,92,REC=1) (FBASE1(I),I=1,80)
 
      C=0
      DO I=1,80
       IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO
      WRITE(90,92,REC=1)(FBASE1(I),I=1,C), (FDIR1(I),I=1,8),
     1                  (FTIME1(I),I=1,11),(FSUBS1(I),I=1,6)
      READ(90,93,REC=1)NAME1
      CLOSE(90)
C       Generate Filename for NESDIS HOURLY Radiation
C       Generate Filename for NESDIS 30 Hour Rotating Radiation
 
        ENDIF
 
      RETURN
      END











































