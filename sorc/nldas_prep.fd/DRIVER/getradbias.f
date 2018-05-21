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
! getradbias.f: 
!
! DESCRIPTION:
!  Opens, reads, interpolates and overlays radiation forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!
! REVISION HISTORY:
!  28  Oct 1999: Brian Cosgrove; Initial code
!  27  Apr 2000: Brian Cosgrove' Disabled zenith angle correction cutoff for 
!                cos(zen) less than .2
!  11  May 2000: Brian Cosgrove; Enabled correction cutoffs for cos(zen) less
!                than .1, stop model if computed value is greater than 1367 w/m2
!  08  Jan 2001: Brian Cosgrove; Added check to see if czb or czm is equal
!                to zero before trying to divide by czm or czb.  If it is
!                zero, set radiation value to zero
!  27  Feb 2001: Brian Cosgrove; Added CZM into call for ZTERP subroutine
!  06  Mar 2001: Brian Cosgrove; Changed computation of WT1 and WT2 in cases
!                where the previous hour or the next hour of observed radiation
!                is not available.  Substituted TIME1 for LDAS%PINKTIME1 and
!                LDAS%NESTIME1 and TIME2 for LDAS%PINKTIME2 and LDAS%NESTIME2
!  05 Feb 2002:  Brian Cosgrove; Added year, month into call for retrad,
!                changed checks for undefined so that doesn't look for N.EQ to
!                -999.9, but look to see if GT or LT that values
!  15 May 2007:  Charles Alonge - Bias Correction modification
!=========================================================================
      SUBROUTINE GETRADBC(LDAS,GRID)

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

!      PRINT*,"TIME1,TIME2,RTIME1,RTIME2"
!      PRINT*,TIME1,TIME2,LDAS%RADBCTIME1,LDAS%RADBCTIME2

!=== Check to see if required bias corr data in memory, if not, read it
       IF(TIME1.NE.LDAS%RADBCTIME1)THEN !We need to get new TIME1 data
        IF(TIME1.EQ.LDAS%RADBCTIME2)THEN !Transfer TIME2 data TIME1	
         LDAS%RADBCTIME1=LDAS%RADBCTIME2
	 LDAS%RADBCSTAT1=LDAS%RADBCSTAT2
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
           GRID(C,R)%RADBCDATA1=GRID(C,R)%RADBCDATA2
          ENDDO
         ENDDO
        ELSE  !Get RAD Data for BC TIME1
!=== The following retreives observed radiation data for Time1
         CALL RADBCFILE(NAME1,LDAS,YR1,MO1,DA1,HR1)
         CALL RETRADBC(1,LDAS,GRID,NAME1,LDAS%RADBCSTAT1)
         IF (LDAS%RADBCSTAT1.NE.0) LDAS%RADBCTIME1=TIME1
!=== End of Data Search For Pinker Time1
        ENDIF
       ENDIF

 

       IF(TIME2.NE.LDAS%RADBCTIME2)THEN !Get RAD Data for TIME2
!=== The following retreives observed radiation data for TIME2
        CALL RADBCFILE(NAME1,LDAS,YR2,MO2,DA2,HR2)
        CALL RETRADBC(2,LDAS,GRID,NAME1,LDAS%RADBCSTAT2)

        IF(LDAS%RADBCSTAT2.NE.0)LDAS%RADBCTIME2=TIME2

         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ( (GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%RADBCDATA2.NE.LDAS%UDEF).AND.
     &      (GRID(C,R)%RADBCDATA2.GT.-998.999) ) THEN
             IF (GRID(C,R)%RADBCDATA2.LT.0.0) THEN
              write (*,'(A21,F15.5,A6,I3,A1,I3,A21)')
     &              'RADBC2 OUT OF RANGE (',
     &          GRID(C,R)%RADBCDATA2,') at (',c,',',r,
     &              ') CORRECTING TO 1.0'
              GRID(C,R)%RADBCDATA2=1.0
             ENDIF

             IF (GRID(C,R)%RADBCDATA2.GT.2.2) THEN
              write (*,'(A21,F15.5,A6,I3,A1,I3,A20)')
     &          'RADBC2 OUT OF RANGE (',
     &          GRID(C,R)%RADBCDATA2,') at (',c,',',r,
     &          ') CORRECTING TO 1'
              GRID(C,R)%RADBCDATA2=1.0
             ENDIF
            ENDIF
        ENDDO
        ENDDO

!=== End of Data Search For Pinker Time2
      ENDIF
     
!      PRINT*,"STATS",LDAS%RADBCSTAT1,LDAS%RADBCSTAT2

      IF(LDAS%TSCOUNT.EQ.1) THEN

         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%RADBCDATA1.NE.LDAS%UDEF).AND.
     &      (GRID(C,R)%RADBCDATA1.GT.0)) THEN
             IF (GRID(C,R)%RADBCDATA1.LT.0.0) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A21)')
     &              'RADBC1 OUT OF RANGE (',
     &          GRID(C,R)%RADBCDATA1,') at (',c,',',r,
     &              ') CORRECTING TO 1.0'
              GRID(C,R)%RADBCDATA1=1.0
             ENDIF

             IF (GRID(C,R)%RADBCDATA1.GT.2.2) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A20)')
     &          'RADBC1 OUT OF RANGE (',
     &          GRID(C,R)%RADBCDATA1,') at (',c,',',r,
     &          ') CORRECTING TO 1'
              GRID(C,R)%RADBCDATA1=1.0
             ENDIF
            ENDIF
        ENDDO
        ENDDO
      ENDIF

!== Compute weights and zenith angle information.  Replace forcing
!== with linearly extrapolated BC data
      WT1=(TIME2-LDAS%TIME)/(TIME2-LDAS%RADBCTIME1)
      WT2=1.0-WT1
     
      DO C=1,LDAS%NC
        DO R=1,LDAS%NR
          GRID(C,R)%FORCING(22)=GRID(C,R)%RADBCDATA1*WT1 +
     &                         GRID(C,R)%RADBCDATA2*WT2
       
           if(C.eq.1.and.R.eq.1) then
           write(*,*) GRID(C,R)%FORCING(22)
           endif  
        ENDDO
      ENDDO

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
! radfile.f: 
!
! DESCRIPTION:
!  This subroutine puts together the radiation data filenames
!
!
! REVISION HISTORY:
!  28  Oct 1999: Brian Cosgrove; Initial code
!=========================================================================
      SUBROUTINE RADBCFILE(NAME1,LDAS,YR,MO,DA,HR)
      USE ldas_module         ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      TYPE (ldasdec) LDAS

!=== Local Variables =====================================================
      CHARACTER*80 NAME1
      INTEGER YR,MO,DA,HR,I,C

      CHARACTER*1  FNAME1(80),FBASE1(80),FSUBS1(4)
      CHARACTER*1  FTIME1(4),FDIR1(1)

      CHARACTER*1  FNAME2(80),FBASE2(80),FSUBS2(80)
      CHARACTER*1  FTIME2(12),FDIR2(8)

!=== End Variable Definition =============================================

!=== Put together filename

 92   FORMAT(80A1)
 93   FORMAT(A80)
 94   FORMAT(I2,I2)
 95   FORMAT(4A1)
 96   FORMAT(A80)
 97   FORMAT(A4) 
 98   FORMAT(A1)
 99   FORMAT(A1)

      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      WRITE(90,98,REC=1) '/'
      READ(90,99,REC=1) FDIR1

      WRITE(90,94,REC=1) MO,HR
      READ(90,95,REC=1)FTIME1
      DO I=1,4
       IF(FTIME1(I).EQ.(' '))FTIME1(I)='0'
      ENDDO

      WRITE(90,97,REC=1)'.bin'
      READ(90,95,REC=1) (FSUBS1(I),I=1,4)

      WRITE(90,96,REC=1) LDAS%RADBCDIR
      READ(90,92,REC=1) (FBASE1(I),I=1,80)

      C=0
      DO I=1,80
       IF(FBASE1(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO
      WRITE(90,92,REC=1)(FBASE1(I),I=1,C), (FDIR1(I),I=1,1),
     1                  (FTIME1(I),I=1,4),(FSUBS1(I),I=1,4)
      READ(90,93,REC=1)NAME1
      CLOSE(90)

      RETURN
      END

