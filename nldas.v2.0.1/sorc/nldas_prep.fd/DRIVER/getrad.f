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
! getrad.f: 
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
!=========================================================================
      SUBROUTINE GETRAD(LDAS,GRID)

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

      IF (LDAS%PINKER.NE.0) THEN
!=== Check to see if required Pinker data in memory, if not, read it
       IF(TIME1.NE.LDAS%PINKTIME1)THEN !We need to get new TIME1 data
        IF(TIME1.EQ.LDAS%PINKTIME2)THEN !Transfer TIME2 data TIME1	
         LDAS%PINKTIME1=LDAS%PINKTIME2
	 LDAS%PSTAT1=LDAS%PSTAT2
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
           GRID(C,R)%PINKDATA1=GRID(C,R)%PINKDATA2
          ENDDO
         ENDDO
        ELSE  !Get RAD Data for Pinker TIME1
!=== The following retreives observed radiation data for Time1
	if ((yr1.ge.1996).and.(yr1.le.2000)) then
         CALL RADFILE(NAME1,LDAS,YR1,MO1,DA1,HR1,1)
         CALL RETRAD(1,LDAS,GRID,NAME1,LDAS%PSTAT1,1,yr1,mo1)
	else
         CALL RADFILEOLD(NAME1,LDAS,YR1,MO1,DA1,HR1,1)
         CALL RETRAD(1,LDAS,GRID,NAME1,LDAS%PSTAT1,1,yr1,mo1)
	endif
         IF (LDAS%PSTAT1.NE.0) LDAS%PINKTIME1=TIME1
!=== End of Data Search For Pinker Time1
        ENDIF
       ENDIF
       IF(TIME2.NE.LDAS%PINKTIME2)THEN !Get RAD Data for TIME2
!=== The following retreives observed radiation data for TIME2
        if ((yr2.ge.1996).and.(yr2.le.2000)) then
        CALL RADFILE(NAME1,LDAS,YR2,MO2,DA2,HR2,1)
        CALL RETRAD(2,LDAS,GRID,NAME1,LDAS%PSTAT2,1,yr2,mo2)
	else
        CALL RADFILEOLD(NAME1,LDAS,YR2,MO2,DA2,HR2,1)
        CALL RETRAD(2,LDAS,GRID,NAME1,LDAS%PSTAT2,1,yr2,mo2)
	endif

        IF(LDAS%PSTAT2.NE.0)LDAS%PINKTIME2=TIME2

         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%PINKDATA2.NE.LDAS%UDEF).AND.
     &      ((GRID(C,R)%PINKDATA2.GT.-998.999).AND.
     &      (GRID(C,R)%PINKDATA2.LT.-999.01))) THEN
             IF (GRID(C,R)%PINKDATA2.LT.0.0) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A21)')
     &              'PINKER OUT OF RANGE (',
     &          GRID(C,R)%PINKDATA2,') at (',c,',',r,
     &              ') CORRECTING TO 0.0'
              GRID(C,R)%PINKDATA2=0.0
             ENDIF

             IF (GRID(C,R)%PINKDATA2.GT.1360) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A20)')
     &          'PINKER OUT OF RANGE (',
     &          GRID(C,R)%PINKDATA2,') at (',c,',',r,
     &          ') CORRECTING TO 1360'
              GRID(C,R)%PINKDATA2=1360.0
             ENDIF
            ENDIF
        ENDDO
        ENDDO

!=== End of Data Search For Pinker Time2
       ENDIF
      ENDIF


      IF (LDAS%NESDIS.NE.0) THEN
!=== Check to see if required NESDIS data in memory, if not, read it
       IF(TIME1.NE.LDAS%NESTIME1)THEN !We need to get new TIME1 data
        IF(TIME1.EQ.LDAS%NESTIME2)THEN !Transfer TIME2 data TIME1
         LDAS%NESTIME1=LDAS%NESTIME2
         LDAS%NESSTAT1=LDAS%NESSTAT2
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
           GRID(C,R)%NESDATA1=GRID(C,R)%NESDATA2
          ENDDO
         ENDDO
        ELSE  !Get NESDIS Data for TIME1
!=== The following retreives observed radiation data for NESDIS Time1
         CALL RADFILE(NAME1,LDAS,YR1,MO1,DA1,HR1,2)
         CALL RETRAD(1,LDAS,GRID,NAME1,LDAS%NESSTAT1,2,yr1,mo1)
         IF(LDAS%NESSTAT1.NE.0)LDAS%NESTIME1=TIME1
!=== End of Data Search FOR NESDIS Data Time 1
        ENDIF
       ENDIF
       IF(TIME2.NE.LDAS%NESTIME2)THEN !Get NESDIS Data for TIME2
!=== The following retreives observed radiation data for NESDIS TIME2
        CALL RADFILE(NAME1,LDAS,YR2,MO2,DA2,HR2,2)
        CALL RETRAD(2,LDAS,GRID,NAME1,LDAS%NESSTAT2,2,yr2,mo2)
        IF(LDAS%NESSTAT2.NE.0)LDAS%NESTIME2=TIME2

         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%NESDATA2.NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%NESDATA2.LT.0.0) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A21)')
     &              'NESDIS OUT OF RANGE (',
     &          GRID(C,R)%NESDATA2,') at (',c,',',r,
     &              ') CORRECTING TO 0.0'
              GRID(C,R)%NESDATA2=0.0
             ENDIF

             IF (GRID(C,R)%NESDATA2.GT.1360) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A20)')
     &          'NESDIS OUT OF RANGE (',
     &          GRID(C,R)%NESDATA2,') at (',c,',',r,
     &          ') CORRECTING TO 1360'
              GRID(C,R)%NESDATA2=1360.0
             ENDIF
            ENDIF
        ENDDO
        ENDDO


!=== End of Data Search For NESDIS Data TIME2
       ENDIF
      ENDIF


      IF (LDAS%BRTTMP.NE.0) THEN
!=== Check to see if required BRTTMP data in memory, if not, read it
       IF(TIME1.NE.LDAS%BRTTMPTIME1)THEN !We need to get new TIME1 data
        IF(TIME1.EQ.LDAS%BRTTMPTIME2)THEN !Transfer TIME2 data TIME1
         LDAS%BRTTMPTIME1=LDAS%BRTTMPTIME2
         LDAS%BRTTMPSTAT1=LDAS%BRTTMPSTAT2
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
           GRID(C,R)%BRTTMPDATA1=GRID(C,R)%BRTTMPDATA2
          ENDDO
         ENDDO
        ELSE  !Get Brightness Temp Data for TIME1
!=== The following retreives observed BRTTMP data for BRTTMP Time1
        if ((yr1.ge.1996).and.(yr1.le.2000)) then
         CALL RADFILEB(NAME1,LDAS,YR1,MO1,DA1,HR1,1)
         CALL RETRAD(1,LDAS,GRID,NAME1,LDAS%BRTTMPSTAT1,4,yr1,mo1)
	else
         CALL RADFILEBOLD(NAME1,LDAS,YR1,MO1,DA1,HR1,1)
         CALL RETRAD(1,LDAS,GRID,NAME1,LDAS%BRTTMPSTAT1,4,yr1,mo1)
	endif

         IF(LDAS%BRTTMPSTAT1.NE.0)LDAS%BRTTMPTIME1=TIME1
!=== End of Data Search FOR BRTTMP Data Time 1
        ENDIF
       ENDIF
       IF(TIME2.NE.LDAS%BRTTMPTIME2)THEN !Get BRTTMP Data for TIME2
!=== The following retreives observed radiation data for BRTTMP TIME2
        if ((yr2.ge.1996).and.(yr2.le.2000)) then
        CALL RADFILEB(NAME1,LDAS,YR2,MO2,DA2,HR2,1)
        CALL RETRAD(2,LDAS,GRID,NAME1,LDAS%BRTTMPSTAT2,4,yr2,mo2)
	else
        CALL RADFILEBOLD(NAME1,LDAS,YR2,MO2,DA2,HR2,1) 
        CALL RETRAD(2,LDAS,GRID,NAME1,LDAS%BRTTMPSTAT2,4,yr2,mo2)
	ENDIF
        IF(LDAS%BRTTMPSTAT2.NE.0)LDAS%BRTTMPTIME2=TIME2

         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%BRTTMPDATA2.NE.LDAS%UDEF).AND.
     &  (GRID(C,R)%BRTTMPDATA2.NE.0).and.
     &  ((GRID(C,R)%BRTTMPDATA2.GT.-998.999).AND.
     &  (GRID(C,R)%BRTTMPDATA2.LT.-999.01))) THEN
             IF (GRID(C,R)%BRTTMPDATA2.LT.213.0) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A21)')
     &              'BRTTMP OUT OF RANGE (',
     &          GRID(C,R)%BRTTMPDATA2,') at (',c,',',r,
     &              ') CORRECTING TO 0.0'
              GRID(C,R)%BRTTMPDATA2=213.0
             ENDIF

             IF (GRID(C,R)%BRTTMPDATA2.GT.343) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A20)')
     &          'BRTTMP OUT OF RANGE (',
     &          GRID(C,R)%BRTTMPDATA2,') at (',c,',',r,
     &          ') CORRECTING TO 343'
              GRID(C,R)%BRTTMPDATA2=343.0
             ENDIF
            ENDIF
        ENDDO
        ENDDO

!=== End of Data Search For BRTTMP Data TIME2
       ENDIF
      ENDIF

        if(LDAS%TSCOUNT.EQ.1) then
      IF (LDAS%NESDIS.NE.0) THEN
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%NESDATA1.NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%NESDATA1.LT.0.0) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A21)')
     &              'NESDIS OUT OF RANGE (',
     &          GRID(C,R)%NESDATA1,') at (',c,',',r,
     &              ') CORRECTING TO 0.0'
              GRID(C,R)%NESDATA1=0.0
             ENDIF

             IF (GRID(C,R)%NESDATA1.GT.1360) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A20)')
     &          'NESDIS OUT OF RANGE (',
     &          GRID(C,R)%NESDATA1,') at (',c,',',r,
     &          ') CORRECTING TO 1360'
              GRID(C,R)%NESDATA1=1360.0
             ENDIF
            ENDIF
        ENDDO
        ENDDO
	ENDIF

      IF (LDAS%PINKER.NE.0) THEN
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%PINKDATA1.NE.LDAS%UDEF).AND.
     &      ((GRID(C,R)%PINKDATA1.GT.-998.999).AND.
     &      (GRID(C,R)%PINKDATA1.LT.-999.01))) THEN
             IF (GRID(C,R)%PINKDATA1.LT.0.0) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A21)')
     &              'PINKER OUT OF RANGE (',
     &          GRID(C,R)%PINKDATA1,') at (',c,',',r,
     &              ') CORRECTING TO 0.0'
              GRID(C,R)%PINKDATA1=0.0
             ENDIF

             IF (GRID(C,R)%PINKDATA1.GT.1360) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A20)')
     &          'PINKER OUT OF RANGE (',
     &          GRID(C,R)%PINKDATA1,') at (',c,',',r,
     &          ') CORRECTING TO 1360'
              GRID(C,R)%PINKDATA1=1360.0
             ENDIF
            ENDIF
        ENDDO
        ENDDO
	ENDIF

      IF (LDAS%BRTTMP.NE.0) THEN
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%BRTTMPDATA1.NE.LDAS%UDEF).AND.
     &  (GRID(C,R)%BRTTMPDATA1.NE.0).and.
     &  ((GRID(C,R)%BRTTMPDATA1.GT.-998.999).AND.
     &  (GRID(C,R)%BRTTMPDATA1.LT.-999.01))) THEN

             IF (GRID(C,R)%BRTTMPDATA1.LT.213.0) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A21)')
     &              'BRTTMP OUT OF RANGE (',
     &          GRID(C,R)%BRTTMPDATA1,') at (',c,',',r,
     &              ') CORRECTING TO 0.0'
              GRID(C,R)%BRTTMPDATA1=213.0
             ENDIF

             IF (GRID(C,R)%BRTTMPDATA1.GT.343) THEN
              write (*,'(A21,F10.5,A6,I3,A1,I3,A20)')
     &          'BRTTMP OUT OF RANGE (',
     &          GRID(C,R)%BRTTMPDATA1,') at (',c,',',r,
     &          ') CORRECTING TO 343'
              GRID(C,R)%BRTTMPDATA1=343.0
             ENDIF
            ENDIF
        ENDDO
        ENDDO
	ENDIF
        endif


!=== Print out Status of data holdings

	IF (LDAS%PINKER.NE.0) THEN	
	IF (LDAS%PSTAT1.EQ.0) WRITE(79,*)'PINKER UNAVAILABLE',
     &  MO1,DA1,YR1,HR1
        IF (LDAS%PSTAT1.EQ.1) WRITE(79,*)'PINKER AVAILABLE',
     &  MO1,DA1,YR1,HR1
        IF (LDAS%PSTAT1.EQ.2) WRITE(79,*)'PINKER UNDEFINED',
     &  MO1,DA1,YR1,HR1

        IF (LDAS%PSTAT2.EQ.0) WRITE(79,*)'PINKER UNAVAILABLE',
     &  MO2,DA2,YR2,HR2
        IF (LDAS%PSTAT2.EQ.1) WRITE(79,*)'PINKER AVAILABLE',
     &  MO2,DA2,YR2,HR2
        IF (LDAS%PSTAT2.EQ.2) WRITE(79,*)'PINKER UNDEFINED',
     &  MO2,DA2,YR2,HR2
	ENDIF

	IF (LDAS%NESDIS.NE.0) THEN
        IF (LDAS%NESSTAT1.EQ.0) WRITE(79,*)'NESDIS UNAVAILABLE',
     &  MO1,DA1,YR1,HR1
        IF (LDAS%NESSTAT1.EQ.1) WRITE(79,*)'NESDIS AVAILABLE',
     &  MO1,DA1,YR1,HR1
        IF (LDAS%NESSTAT1.EQ.2) WRITE(79,*)'NESDIS UNDEFINED',
     &  MO1,DA1,YR1,HR1

        IF (LDAS%NESSTAT2.EQ.0) WRITE(79,*)'NESDIS UNAVAILABLE',
     &  MO2,DA2,YR2,HR2
        IF (LDAS%NESSTAT2.EQ.1) WRITE(79,*)'NESDIS AVAILABLE',
     &  MO2,DA2,YR2,HR2
        IF (LDAS%NESSTAT2.EQ.2) WRITE(79,*)'NESDIS UNDEFINED',
     &  MO2,DA2,YR2,HR2
	ENDIF

        IF (LDAS%BRTTMP.NE.0) THEN
        IF (LDAS%BRTTMPSTAT1.EQ.0) WRITE(79,*)'BRTTMP UNAVAILABLE',
     &  MO1,DA1,YR1,HR1
        IF (LDAS%BRTTMPSTAT1.EQ.1) WRITE(79,*)'BRTTMP AVAILABLE',
     &  MO1,DA1,YR1,HR1
        IF (LDAS%BRTTMPSTAT1.EQ.2) WRITE(79,*)'BRTTMP UNDEFINED',
     &  MO1,DA1,YR1,HR1

        IF (LDAS%BRTTMPSTAT2.EQ.0) WRITE(79,*)'BRTTMP UNAVAILABLE',
     &  MO2,DA2,YR2,HR2
        IF (LDAS%BRTTMPSTAT2.EQ.1) WRITE(79,*)'BRTTMP AVAILABLE',
     &  MO2,DA2,YR2,HR2
        IF (LDAS%BRTTMPSTAT2.EQ.2) WRITE(79,*)'BRTTMP UNDEFINED',
     &  MO2,DA2,YR2,HR2
        ENDIF

! Setup fsource array

	LDAS%FSOURCE(7)=0
	LDAS%FSOURCE(8)=0

	IF ((LDAS%NESSTAT1.EQ.1).OR.(LDAS%NESSTAT2.EQ.1)) THEN
          IF (LDAS%NESDIS.EQ.1) LDAS%FSOURCE(8)=1
        ENDIF 

        IF ((LDAS%NESSTAT1.EQ.1).OR.(LDAS%NESSTAT2.EQ.1)) THEN
          IF (LDAS%PINKER.EQ.1) THEN
            IF ((LDAS%PSTAT1.EQ.0).AND.(LDAS%PSTAT2.EQ.0)) THEN
             IF (LDAS%NESDIS.EQ.2) LDAS%FSOURCE(8)=1
            ENDIF
          ENDIF 
        ENDIF 

        IF ((LDAS%PSTAT1.EQ.1).OR.(LDAS%PSTAT2.EQ.1)) THEN
          IF (LDAS%PINKER.EQ.1) LDAS%FSOURCE(7)=1
        ENDIF 

        IF ((LDAS%PSTAT1.EQ.1).OR.(LDAS%PSTAT2.EQ.1)) THEN
          IF (LDAS%NESDIS.EQ.1) THEN
            IF ((LDAS%NESSTAT1.EQ.0).AND.(LDAS%NESSTAT2.EQ.0)) THEN
             IF (LDAS%PINKER.EQ.2) LDAS%FSOURCE(7)=1
            ENDIF
          ENDIF
        ENDIF 


c=== If either BRTTMP data is available, set
c=== up looping structure for forcing replacement process
           IF (LDAS%BRTTMP.EQ.1) THEN
           IF((LDAS%BRTTMPSTAT1.EQ.1).AND.(LDAS%BRTTMPSTAT2.EQ.1))THEN
!== Compute weights and BRTTMP data.
            WT1=(LDAS%BRTTMPTIME2-LDAS%TIME)/(LDAS%BRTTMPTIME2-
     &       LDAS%BRTTMPTIME1)
            WT2=1.0-WT1
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
              GRID(C,R)%OBSBT=LDAS%UDEF
                GRID(C,R)%OBSBT=GRID(C,R)%BRTTMPDATA1*WT1+
     1           GRID(C,R)%BRTTMPDATA2*WT2
              IF ((GRID(C,R)%BRTTMPDATA1.GT.0.0).AND.
     &                  (GRID(C,R)%BRTTMPDATA2.LE.0.0)) THEN
                GRID(C,R)%OBSBT=LDAS%UDEF
c                GRID(C,R)%OBSBT=GRID(C,R)%BRTTMPDATA1*WT1+
c     1                             0.0*WT2
              ENDIF
              IF ((GRID(C,R)%BRTTMPDATA1.LE.0.0).AND.
     &                  (GRID(C,R)%BRTTMPDATA2.GT.0.0)) THEN
                GRID(C,R)%OBSBT=LDAS%UDEF
c                GRID(C,R)%OBSBT=0.0*WT1+
c     1                             GRID(C,R)%BRTTMPDATA2*WT2
              ENDIF

              IF ((GRID(C,R)%BRTTMPDATA1.LE.0.0).AND.
     &                  (GRID(C,R)%BRTTMPDATA2.LE.0.0)) THEN
                GRID(C,R)%OBSBT=LDAS%UDEF
              ENDIF

	if ((GRID(C,R)%OBSBT.ne.LDAS%UDEF).and.
     &  (GRID(C,R)%OBSBT.lt.150)) then	
	print *,'doh, obst,1,2',GRID(C,R)%OBSBT,
     &  GRID(C,R)%BRTTMPDATA1,GRID(C,R)%BRTTMPDATA2
	stop
	endif

             ENDDO
            ENDDO
           ENDIF


          IF((LDAS%BRTTMPSTAT1.EQ.1).AND.(LDAS%BRTTMPSTAT2.NE.1))THEN
!== Compute weights and BRTTMP DATA.
           WT1=(TIME2-LDAS%TIME)/(TIME2-
     &     LDAS%BRTTMPTIME1)
           WT2=1.0-WT1
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             GRID(C,R)%OBSBT=LDAS%UDEF
c            IF (GRID(C,R)%BRTTMPDATA1.GT.0.0) THEN
c               GRID(C,R)%OBSBT=GRID(C,R)%BRTTMPDATA1*WT1
c             ENDIF
            ENDDO
           ENDDO
          ENDIF


          IF((LDAS%BRTTMPSTAT1.NE.1).AND.(LDAS%BRTTMPSTAT2.EQ.1))THEN
!== Compute weights and BRTTMP DATA.
           WT1=(LDAS%BRTTMPTIME2-LDAS%TIME)/(LDAS%BRTTMPTIME2-
     &     TIME1)
           WT2=1.0-WT1
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             GRID(C,R)%OBSBT=LDAS%UDEF
c             IF (GRID(C,R)%BRTTMPDATA2.GT.0.0) THEN
c               GRID(C,R)%OBSBT=0.0*WT1+
c     1         GRID(C,R)%BRTTMPDATA2*WT2
c             ENDIF
            ENDDO
           ENDDO
          ENDIF
          IF((LDAS%BRTTMPSTAT1.NE.1).AND.(LDAS%BRTTMPSTAT2.NE.1))THEN
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             GRID(C,R)%OBSBT=LDAS%UDEF
            ENDDO
           ENDDO
          ENDIF
	ENDIF !endif the brttmp=1

c=== If either NESDIS or PINKER data is available, set
c=== up looping structure for forcing replacement process

	IF ((LDAS%NESDIS.NE.0).OR.(LDAS%PINKER.NE.0)) THEN
	 IF ((LDAS%PINKER.EQ.1).AND.(LDAS%NESDIS.EQ.2)) THEN
	  AA=2
	  BB=1
	  CC=-1
	 ENDIF
         IF ((LDAS%PINKER.EQ.2).AND.(LDAS%NESDIS.EQ.1)) THEN
          AA=1
	  BB=2
          CC=1
         ENDIF
         IF ((LDAS%PINKER.NE.0).AND.(LDAS%NESDIS.EQ.0)) THEN
          AA=1
          BB=1
	  CC=1
         ENDIF
         IF ((LDAS%PINKER.EQ.0).AND.(LDAS%NESDIS.NE.0)) THEN
          AA=2
          BB=2
          CC=1
         ENDIF
!== Loop through and replace data as possible with NESDIS/PINKER data
!== This depends on options specified in ldas.crd as well as actual
!== data holdings.  The primary data set is looped through second, and
!== is placed on top of the secondary data set as well as the original
!== NCEP radiation data.
	 DO LOOP=AA,BB,CC
	  IF (LOOP.EQ.1) THEN
	   IF ((LDAS%PSTAT1.EQ.1).AND.(LDAS%PSTAT2.EQ.1)) THEN
!== Compute weights and zenith angle information.  Replace forcing 
!== with Pinker data.
            WT1=(LDAS%PINKTIME2-LDAS%TIME)/(LDAS%PINKTIME2-
     &       LDAS%PINKTIME1)
            WT2=1.0-WT1
     

            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
              ZDOY=LDAS%DOY
              CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1         GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2         CZM,LDAS,GRID,czmean)
              GRID(C,R)%OBSW=LDAS%UDEF
              IF ((GRID(C,R)%PINKDATA1.GT.0.0).AND.
     &          (GRID(C,R)%PINKDATA2.GT.0.0)) THEN
                GRID(C,R)%OBSW=GRID(C,R)%PINKDATA1*ZW1+
     1          GRID(C,R)%PINKDATA2*ZW2
               IF((GRID(C,R)%OBSW.GT.GRID(C,R)%PINKDATA1.AND.
     1          GRID(C,R)%OBSW.GT.GRID(C,R)%PINKDATA2).AND.
     2          (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%OBSW=GRID(C,R)%PINKDATA1*WT1+
     1           GRID(C,R)%PINKDATA2*WT2
               ENDIF
              ENDIF
              IF ((GRID(C,R)%PINKDATA1.GT.0.0).AND.
     &          (GRID(C,R)%PINKDATA2.LE.0.0)) THEN
                IF (CZB.GT.0.0) THEN
                   GRID(C,R)%OBSW=GRID(C,R)%PINKDATA1*CZM/CZB
                ELSE
                   GRID(C,R)%OBSW=GRID(C,R)%PINKDATA1*0.0
                ENDIF
               IF((GRID(C,R)%OBSW.GT.GRID(C,R)%PINKDATA1.AND.
     1             GRID(C,R)%OBSW.GT.0.0).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%OBSW=GRID(C,R)%PINKDATA1*WT1+
     1                             0.0*WT2
               ENDIF
              ENDIF
              IF ((GRID(C,R)%PINKDATA1.LE.0.0).AND.
     &          (GRID(C,R)%PINKDATA2.GT.0.0)) THEN
                IF (CZE.GT.0.0) THEN
                   GRID(C,R)%OBSW= GRID(C,R)%PINKDATA2*CZM/CZE
	        ELSE
                   GRID(C,R)%OBSW= GRID(C,R)%PINKDATA2*0.0
                ENDIF 
               IF((GRID(C,R)%OBSW.GT.0.0.AND.
     1             GRID(C,R)%OBSW.GT.GRID(C,R)%PINKDATA2).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%OBSW=0.0*WT1+
     1                             GRID(C,R)%PINKDATA2*WT2
               ENDIF
              ENDIF
              IF ((GRID(C,R)%PINKDATA1.GT.0.0).AND.
     &                  (GRID(C,R)%PINKDATA2.GT.0.0)) THEN
               GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA1*ZW1+
     1                            GRID(C,R)%PINKDATA2*ZW2
               IF((GRID(C,R)%FORCING(3).GT.GRID(C,R)%PINKDATA1.AND.
     1             GRID(C,R)%FORCING(3).GT.GRID(C,R)%PINKDATA2).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA1*WT1+
     1                             GRID(C,R)%PINKDATA2*WT2
               ENDIF
              ENDIF
              IF ((GRID(C,R)%PINKDATA1.GT.0.0).AND.
     &                  (GRID(C,R)%PINKDATA2.LE.0.0)) THEN
              IF (CZB.GT.0.0) THEN
                 GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA1*CZM/CZB
              ELSE
                 GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA1*0.0
              ENDIF
               IF((GRID(C,R)%FORCING(3).GT.GRID(C,R)%PINKDATA1.AND.
     1             GRID(C,R)%FORCING(3).GT.0.0).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA1*WT1+
     1                             0.0*WT2
               ENDIF
              ENDIF
              IF ((GRID(C,R)%PINKDATA1.LE.0.0).AND.
     &                  (GRID(C,R)%PINKDATA2.GT.0.0)) THEN
              IF (CZE.GT.0.0) THEN
                 GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA2*CZM/CZE
              ELSE
                 GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA2*0.0
              ENDIF 
               IF((GRID(C,R)%FORCING(3).GT.0.0.AND.
     1             GRID(C,R)%FORCING(3).GT.GRID(C,R)%PINKDATA2).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%FORCING(3)=0.0*WT1+
     1                             GRID(C,R)%PINKDATA2*WT2
               ENDIF
              ENDIF
              if (GRID(C,R)%FORCING(3).gt.1367) then
               print *,'warning, OBSERVED SW RADIATION TOO HIGH'
               print *,'it is',GRID(C,R)%FORCING(3),' at',c,r
	       print *,'pink1=',GRID(C,R)%PINKDATA1
	       print *,'pink2=',GRID(C,R)%PINKDATA2
	       print *,'wt1,wt2,czb,cze,czm,zw1,zw2'
	       print *,wt1,wt2,czb,cze,czm,zw1,zw2
               stop
              endif
             ENDDO
            ENDDO
           ENDIF

          IF ((LDAS%PSTAT1.EQ.1).AND.(LDAS%PSTAT2.NE.1)) THEN
!== Compute weights and zenith angle information.  Replace forcing 
!== with zenith extrapolated PINKER data
           WT1=(TIME2-LDAS%TIME)/(TIME2-
     &     LDAS%PINKTIME1)
           WT2=1.0-WT1
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             ZDOY=LDAS%DOY
             CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1             GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2             CZM,LDAS,GRID,czmean)
             GRID(C,R)%OBSW=LDAS%UDEF
            IF (GRID(C,R)%PINKDATA1.GT.0.0) THEN
            IF (CZB.GT.0.0) THEN 
               GRID(C,R)%OBSW=GRID(C,R)%PINKDATA1*CZM/CZB
            ELSE
               GRID(C,R)%OBSW=GRID(C,R)%PINKDATA1*0.0
            ENDIF 
              IF((GRID(C,R)%OBSW.GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
               GRID(C,R)%OBSW=GRID(C,R)%PINKDATA1*WT1

              ENDIF
              IF (CZB.GT.0.0) THEN
                 GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA1*CZM/CZB
              ELSE
                 GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA1*0.0
              ENDIF
              IF((GRID(C,R)%FORCING(3).GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2            (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
               GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA1*WT1

              ENDIF
              IF (GRID(C,R)%FORCING(3).gt.1367) then
              print *,'warning, OBSERVED SW RADIATION TOO HIGH'
               print *,'it is',GRID(C,R)%FORCING(3)
               print *,'pink1=',GRID(C,R)%PINKDATA1
               print *,'pink2=',GRID(C,R)%PINKDATA2
               print *,'wt1,wt2,czb,cze,czm,zw1,zw2'
               print *,wt1,wt2,czb,cze,czm,zw1,zw2
               stop
              ENDIF
             ENDIF
            ENDDO
           ENDDO
          ENDIF


          IF ((LDAS%PSTAT1.NE.1).AND.(LDAS%PSTAT2.EQ.1)) THEN
!== Compute weights and zenith angle information.  Replace forcing 
!== with zenith extrapolated PINKER data
           WT1=(LDAS%PINKTIME2-LDAS%TIME)/(LDAS%PINKTIME2-
     &     TIME1)
           WT2=1.0-WT1
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             ZDOY=LDAS%DOY
             CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1       GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2       CZM,LDAS,GRID,czmean)
             GRID(C,R)%OBSW=LDAS%UDEF
             IF (GRID(C,R)%PINKDATA2.GT.0.0) THEN
             IF (CZE.GT.0.0) THEN
                GRID(C,R)%OBSW=GRID(C,R)%PINKDATA2*CZM/CZE
             ELSE
                GRID(C,R)%OBSW=GRID(C,R)%PINKDATA2*0.0
             ENDIF
              IF((GRID(C,R)%OBSW.GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2        (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
               GRID(C,R)%OBSW=0.0*WT1+
     1         GRID(C,R)%PINKDATA2*WT2

              ENDIF
              IF (CZE.GT.0.0) THEN
                 GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA2*CZM/CZE
              ELSE
                 GRID(C,R)%FORCING(3)=GRID(C,R)%PINKDATA2*0.0
              ENDIF 
              IF((GRID(C,R)%FORCING(3).GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2        (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
               GRID(C,R)%FORCING(3)=0.0*WT1+
     1         GRID(C,R)%PINKDATA2*WT2

              ENDIF
              IF (GRID(C,R)%FORCING(3).gt.1367) then
              print *,'warning, OBSERVED SW RADIATION TOO HIGH'
               print *,'it is',GRID(C,R)%FORCING(3)
               print *,'pink1=',GRID(C,R)%PINKDATA1
               print *,'pink2=',GRID(C,R)%PINKDATA2
               print *,'wt1,wt2,czb,cze,czm,zw1,zw2,czm/cze'
               print *,wt1,wt2,czb,cze,czm,zw1,zw2,czm/cze
               stop
              ENDIF
             ENDIF
            ENDDO
           ENDDO
          ENDIF

          IF ((LDAS%PSTAT1.NE.1).AND.(LDAS%PSTAT2.NE.1)) THEN
           DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             GRID(C,R)%OBSW=LDAS%UDEF
            ENDDO
           ENDDO
          ENDIF
	 ENDIF   ! END 'IF LOOP.EQ.1' loop

          IF (LOOP.EQ.2) THEN
           IF ((LDAS%NESSTAT1.EQ.1).AND.(LDAS%NESSTAT2.EQ.1)) THEN
!== Compute weights and zenith angle information.  Replace forcing 
!== with NESDIS data.
            WT1=(LDAS%NESTIME2-LDAS%TIME)/(LDAS%NESTIME2-
     &       LDAS%NESTIME1)
            WT2=1.0-WT1
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
              ZDOY=LDAS%DOY
              CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1              GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2              CZM,LDAS,GRID,czmean)
              GRID(C,R)%OBSW=LDAS%UDEF
              IF ((GRID(C,R)%NESDATA1.GT.0.0).AND.
     &          (GRID(C,R)%NESDATA2.GT.0.0)) THEN
                GRID(C,R)%OBSW=GRID(C,R)%NESDATA1*ZW1+
     1          GRID(C,R)%NESDATA2*ZW2
               IF((GRID(C,R)%OBSW.GT.GRID(C,R)%NESDATA1.AND.
     1             GRID(C,R)%OBSW.GT.GRID(C,R)%NESDATA2).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%OBSW=GRID(C,R)%NESDATA1*WT1+
     1                             GRID(C,R)%NESDATA2*WT2
               ENDIF
	      ENDIF

              IF ((GRID(C,R)%NESDATA1.GT.0.0).AND.
     &          (GRID(C,R)%NESDATA2.LE.0.0)) THEN
              IF (CZB.GT.0.0) THEN
                 GRID(C,R)%OBSW=GRID(C,R)%NESDATA1*CZM/CZB
              ELSE
                 GRID(C,R)%OBSW=GRID(C,R)%NESDATA1*0.0
              ENDIF
               IF((GRID(C,R)%OBSW.GT.GRID(C,R)%NESDATA1.AND.
     1             GRID(C,R)%OBSW.GT.0.0).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%OBSW=GRID(C,R)%NESDATA1*WT1+
     1                             0.0*WT2
               ENDIF

              ENDIF
              IF ((GRID(C,R)%NESDATA1.LE.0.0).AND.
     &          (GRID(C,R)%NESDATA2.GT.0.0)) THEN
              IF (CZE.GT.0.0) THEN
                 GRID(C,R)%OBSW= GRID(C,R)%NESDATA2*CZM/CZE
	      ELSE
                 GRID(C,R)%OBSW= GRID(C,R)%NESDATA2*0.0
              ENDIF
               IF((GRID(C,R)%OBSW.GT.0.0.AND.
     1             GRID(C,R)%OBSW.GT.GRID(C,R)%NESDATA2).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1)) THEN
                GRID(C,R)%OBSW=0.0*WT1+
     1                             GRID(C,R)%NESDATA2*WT2
               ENDIF

              ENDIF



              IF ((GRID(C,R)%NESDATA1.GT.0.0).AND.
     &                  (GRID(C,R)%NESDATA2.GT.0.0)) THEN
               GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA1*ZW1+
     1                            GRID(C,R)%NESDATA2*ZW2

               IF((GRID(C,R)%FORCING(3).GT.GRID(C,R)%NESDATA1.AND.
     1             GRID(C,R)%FORCING(3).GT.GRID(C,R)%NESDATA2).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA1*WT1+
     1                             GRID(C,R)%NESDATA2*WT2
               ENDIF
	       ENDIF

              IF ((GRID(C,R)%NESDATA1.GT.0.0).AND.
     &                  (GRID(C,R)%NESDATA2.LE.0.0)) THEN
              IF (CZB.GT.0.0) THEN 
                 GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA1*CZM/CZB
              ELSE
                 GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA1*0.0
              ENDIF

               IF((GRID(C,R)%FORCING(3).GT.GRID(C,R)%NESDATA1.AND.
     1             GRID(C,R)%FORCING(3).GT.0.0).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA1*WT1+
     1                             0.0*WT2
               ENDIF
	       ENDIF
              IF ((GRID(C,R)%NESDATA1.LE.0.0).AND.
     &                  (GRID(C,R)%NESDATA2.GT.0.0)) THEN
              IF (CZE.GT.0.0) THEN
                 GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA2*CZM/CZE
	      ELSE
                 GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA2*0.0
              ENDIF

               IF((GRID(C,R)%FORCING(3).GT.0.0.AND.
     1             GRID(C,R)%FORCING(3).GT.GRID(C,R)%NESDATA2).AND.
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%FORCING(3)=0.0*WT1+
     1                             GRID(C,R)%NESDATA2*WT2
               ENDIF
	       ENDIF


                 if (GRID(C,R)%FORCING(3).gt.1367) then
                  print *,'warning, OBSERVED SW RADIATION TOO HIGH'
                  print *,'it is',GRID(C,R)%FORCING(3)
                  stop
                 endif

             ENDDO
            ENDDO
           ENDIF

          IF ((LDAS%NESSTAT1.EQ.1).AND.(LDAS%NESSTAT2.NE.1)) THEN
            WT1=(TIME2-LDAS%TIME)/(TIME2-
     &       LDAS%NESTIME1)
            WT2=1.0-WT1
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
              ZDOY=LDAS%DOY
              CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1              GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2              CZM,LDAS,GRID,czmean)
              GRID(C,R)%OBSW=LDAS%UDEF
             IF (GRID(C,R)%NESDATA1.GT.0.0) THEN
             IF (CZB.GT.0.0) THEN
                GRID(C,R)%OBSW=GRID(C,R)%NESDATA1*CZM/CZB
             ELSE
	        GRID(C,R)%OBSW=GRID(C,R)%NESDATA1*0.0
             ENDIF
               IF((GRID(C,R)%OBSW.GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%OBSW=GRID(C,R)%NESDATA1*WT1

               ENDIF
               IF (CZB.GT.0.0) THEN
                  GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA1*CZM/CZB
               ELSE
                  GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA1*0.0
               ENDIF 
               IF((GRID(C,R)%FORCING(3).GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA1*WT1

               ENDIF

                 if (GRID(C,R)%FORCING(3).gt.1367) then
                  print *,'warning, OBSERVED SW RADIATION TOO HIGH!!'
                  print *,'it is',GRID(C,R)%FORCING(3)
                  stop
                 endif


                ENDIF

	     ENDDO
	    ENDDO
       	  ENDIF
          IF ((LDAS%NESSTAT1.NE.1).AND.(LDAS%NESSTAT2.EQ.1)) THEN
            WT1=(LDAS%NESTIME2-LDAS%TIME)/(LDAS%NESTIME2-
     &       TIME1)
            WT2=1.0-WT1
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
              ZDOY=LDAS%DOY
              CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     1              GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,
     2              CZM,LDAS,GRID,czmean)
              GRID(C,R)%OBSW=LDAS%UDEF
              IF (GRID(C,R)%NESDATA2.GT.0.0) THEN 
              IF (CZE.GT.0.0) THEN
                 GRID(C,R)%OBSW=GRID(C,R)%NESDATA2*CZM/CZE
	      ELSE
                 GRID(C,R)%OBSW=GRID(C,R)%NESDATA2*0.0
              ENDIF
               IF((GRID(C,R)%OBSW.GT.400.0).AND.
!       Arbitrary cutoff value of 400 W/m2
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%OBSW=0.0*WT1+
     1                             GRID(C,R)%NESDATA2*WT2

               ENDIF

               IF (CZE.GT.0.0) THEN
                  GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA2*CZM/CZE
               ELSE
                  GRID(C,R)%FORCING(3)=GRID(C,R)%NESDATA2*0.0
               ENDIF
               IF((GRID(C,R)%FORCING(3).GT.400.0).AND.
!	Arbitrary cutoff value of 400 W/m2
     2             (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
                GRID(C,R)%FORCING(3)=0.0*WT1+
     1                             GRID(C,R)%NESDATA2*WT2

               ENDIF

                 if (GRID(C,R)%FORCING(3).gt.1367) then
                  print *,'warning, OBSERVED SW RADIATION TOO HIGH!!'
                  print *,'it is',GRID(C,R)%FORCING(3)
                  stop
                 endif


		ENDIF
             ENDDO
            ENDDO
          ENDIF

          IF ((LDAS%NESSTAT1.NE.1).AND.(LDAS%NESSTAT2.NE.1)) THEN
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
              GRID(C,R)%OBSW=LDAS%UDEF
             ENDDO
            ENDDO
          ENDIF
         ENDIF   ! END 'IF LOOP.EQ.2' loop
        ENDDO
       ENDIF
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
! radfile.f: 
!
! DESCRIPTION:
!  This subroutine puts together the radiation data filenames
!
!
! REVISION HISTORY:
!  28  Oct 1999: Brian Cosgrove; Initial code
!=========================================================================
      SUBROUTINE RADFILE(NAME1,LDAS,YR,MO,DA,HR,FLAG)
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
C	Generate Filename for Pinker Radiation
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

      WRITE(90,97,REC=1)'.sda.i'
      READ(90,92,REC=1) (FSUBS1(I),I=1,6)

      WRITE(90,96,REC=1) LDAS%PINKDIR
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

 72   FORMAT(80A1)
 73   FORMAT(A80)
 74   FORMAT(I4,I2,I2,I2)
 75   FORMAT(10A1)
 76   FORMAT(A40)
 77   FORMAT(A5)
 78   FORMAT(A1,I4,I2,A1)
 79   FORMAT(8A1)

        IF (FLAG.EQ.2) THEN
C       Generate Filename for NESDIS Radiation
      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      WRITE(90,78,REC=1)'/',YR,MO,'/'
      READ(90,79,REC=1)FDIR2
      DO I=1,8
       IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO

      WRITE(90,74,REC=1)YR,MO,DA,HR
      READ(90,75,REC=1)FTIME2
      DO I=1,10
       IF(FTIME2(I).EQ.(' '))FTIME2(I)='0'
      ENDDO

      WRITE(90,77,REC=1)'.GOES'
      READ(90,72,REC=1) (FSUBS2(I),I=1,5)

      WRITE(90,76,REC=1) LDAS%NDDIR
      READ(90,72,REC=1) (FBASE2(I),I=1,40)

      C=0
      DO I=1,40
       IF(FBASE2(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO
      WRITE(90,72,REC=1)(FBASE2(I),I=1,C), (FDIR2(I),I=1,8),
     1                  (FTIME2(I),I=1,10),(FSUBS2(I),I=1,5)
      READ(90,73,REC=1)NAME1
      CLOSE(90)

C       Generate Filename for NESDIS HOURLY Radiation
C       Generate Filename for NESDIS 30 Hour Rotating Radiation

        ENDIF

      RETURN
      END

C=====================================================================
C        -----------------RADFILEB subroutine

      SUBROUTINE RADFILEB(NAME1,LDAS,YR,MO,DA,HR,FLAG)
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
C       Generate Filename for Pinker Radiation
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

      WRITE(90,97,REC=1)'.stp.i'
      READ(90,92,REC=1) (FSUBS1(I),I=1,6)

      WRITE(90,96,REC=1) LDAS%BRTTMPDIR
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

 72   FORMAT(80A1)
 73   FORMAT(A80)
 74   FORMAT(I4,I2,I2,I2)
 75   FORMAT(10A1)
 76   FORMAT(A40)
 77   FORMAT(A5)
 78   FORMAT(A1,I4,I2,A1)
 79   FORMAT(8A1)

        IF (FLAG.EQ.2) THEN
C       Generate Filename for NESDIS Radiation
      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

      WRITE(90,78,REC=1)'/',YR,MO,'/'
      READ(90,79,REC=1)FDIR2
      DO I=1,8
       IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO

      WRITE(90,74,REC=1)YR,MO,DA,HR
      READ(90,75,REC=1)FTIME2
      DO I=1,10
       IF(FTIME2(I).EQ.(' '))FTIME2(I)='0'
      ENDDO

      WRITE(90,77,REC=1)'.GOES'
      READ(90,72,REC=1) (FSUBS2(I),I=1,5)

      WRITE(90,76,REC=1) LDAS%NDDIR
      READ(90,72,REC=1) (FBASE2(I),I=1,40)

      C=0
      DO I=1,40
       IF(FBASE2(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO
      WRITE(90,72,REC=1)(FBASE2(I),I=1,C), (FDIR2(I),I=1,8),
     1                  (FTIME2(I),I=1,10),(FSUBS2(I),I=1,5)
      READ(90,73,REC=1)NAME1
      CLOSE(90)

C       Generate Filename for NESDIS HOURLY Radiation
C       Generate Filename for NESDIS 30 Hour Rotating Radiation

        ENDIF

      RETURN
      END



      SUBROUTINE RADFILEOLD(NAME1,LDAS,YR,MO,DA,HR,FLAG)
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
C       Generate Filename for Pinker Radiation
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
 
      WRITE(90,97,REC=1)'.sda.i'
      READ(90,92,REC=1) (FSUBS1(I),I=1,6)
 
      WRITE(90,96,REC=1) LDAS%PINKDIROLD
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
 
 72   FORMAT(80A1)
 73   FORMAT(A80)
 74   FORMAT(I4,I2,I2,I2)
 75   FORMAT(10A1)
 76   FORMAT(A40)
 77   FORMAT(A5)
 78   FORMAT(A1,I4,I2,A1)
 79   FORMAT(8A1)
 
        IF (FLAG.EQ.2) THEN
C       Generate Filename for NESDIS Radiation
      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
 
      WRITE(90,78,REC=1)'/',YR,MO,'/'
      READ(90,79,REC=1)FDIR2
      DO I=1,8
       IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO
 
      WRITE(90,74,REC=1)YR,MO,DA,HR
      READ(90,75,REC=1)FTIME2
      DO I=1,10
       IF(FTIME2(I).EQ.(' '))FTIME2(I)='0'
      ENDDO
 
      WRITE(90,77,REC=1)'.GOES'
      READ(90,72,REC=1) (FSUBS2(I),I=1,5)
 
      WRITE(90,76,REC=1) LDAS%NDDIR
      READ(90,72,REC=1) (FBASE2(I),I=1,40)
 
      C=0
      DO I=1,40
       IF(FBASE2(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO
      WRITE(90,72,REC=1)(FBASE2(I),I=1,C), (FDIR2(I),I=1,8),
     1                  (FTIME2(I),I=1,10),(FSUBS2(I),I=1,5)
      READ(90,73,REC=1)NAME1
      CLOSE(90)
 
C       Generate Filename for NESDIS HOURLY Radiation
C       Generate Filename for NESDIS 30 Hour Rotating Radiation
 
        ENDIF
 
      RETURN
      END
 
C=====================================================================
C        -----------------RADFILEB subroutine
 
      SUBROUTINE RADFILEBOLD(NAME1,LDAS,YR,MO,DA,HR,FLAG)
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
C       Generate Filename for Pinker Radiation
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
 
      WRITE(90,97,REC=1)'.stp.i'
      READ(90,92,REC=1) (FSUBS1(I),I=1,6)
 
      WRITE(90,96,REC=1) LDAS%BRTTMPDIROLD
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
 
 72   FORMAT(80A1)
 73   FORMAT(A80)
 74   FORMAT(I4,I2,I2,I2)
 75   FORMAT(10A1)
 76   FORMAT(A40)
 77   FORMAT(A5)
 78   FORMAT(A1,I4,I2,A1)
 79   FORMAT(8A1)
 
        IF (FLAG.EQ.2) THEN
C       Generate Filename for NESDIS Radiation
      OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
 
      WRITE(90,78,REC=1)'/',YR,MO,'/'
      READ(90,79,REC=1)FDIR2
      DO I=1,8
       IF(FDIR2(I).EQ.(' '))FDIR2(I)='0'
      ENDDO
 
      WRITE(90,74,REC=1)YR,MO,DA,HR
      READ(90,75,REC=1)FTIME2
      DO I=1,10
       IF(FTIME2(I).EQ.(' '))FTIME2(I)='0'
      ENDDO
 
      WRITE(90,77,REC=1)'.GOES'
      READ(90,72,REC=1) (FSUBS2(I),I=1,5)
 
      WRITE(90,76,REC=1) LDAS%NDDIR
      READ(90,72,REC=1) (FBASE2(I),I=1,40)
 
      C=0
      DO I=1,40
       IF(FBASE2(I).EQ.(' ').AND.C.EQ.0)C=I-1
      ENDDO
      WRITE(90,72,REC=1)(FBASE2(I),I=1,C), (FDIR2(I),I=1,8),
     1                  (FTIME2(I),I=1,10),(FSUBS2(I),I=1,5)
      READ(90,73,REC=1)NAME1
      CLOSE(90)
 
C       Generate Filename for NESDIS HOURLY Radiation
C       Generate Filename for NESDIS 30 Hour Rotating Radiation
 
        ENDIF
 
      RETURN
      END






































