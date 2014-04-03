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
! ticktime.f: 
!
! DESCRIPTION:
!  Advance the time 1 timestep.
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Version
!  26 Apr 2000: Brian Cosgrove; Fix for Alpha (Used NINT instead
!               of DINT b/c of rounding/truncating issues)
!  01 May 2000: Brian Cosgrove; Correction in TIME2DATE subroutine,
!               if hour=24, hour set to 0, day, month and year advanced
!=========================================================================

      SUBROUTINE TICKTIME (LDAS)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      type (ldasdec) LDAS              

!=== LOCAL Variables =====================================================
      INTEGER DAYS(13)
      DATA DAYS/31,28,31,30,31,30,31,31,30,31,30,31,31/

!=== End Variable List ===================================================

      LDAS%PDA=LDAS%DA  !Used to determine end-of-day restart writes

      LDAS%SS=LDAS%SS+LDAS%TS

      DO WHILE(LDAS%SS.GT.59)
       LDAS%SS=LDAS%SS-60
       LDAS%MN=LDAS%MN+1
      ENDDO 

      DO WHILE(LDAS%MN.GT.59)
       LDAS%MN=LDAS%MN-60
       LDAS%HR=LDAS%HR+1
      ENDDO 

      DO WHILE(LDAS%HR.GE.24)
       LDAS%HR=LDAS%HR-24
       LDAS%DA=LDAS%DA+1
      ENDDO 
 
      IF((MOD(LDAS%YR,4).EQ.0.AND.MOD(LDAS%YR,100).NE.0)   !correct for leap year
     &    .OR.(MOD(LDAS%YR,400).EQ.0))THEN                 !correct for Y2K
       DAYS(2)=29                  
      ELSE
       DAYS(2)=28
      ENDIF

      DO WHILE(LDAS%DA.GT.DAYS(LDAS%MO))
       LDAS%DA=LDAS%DA-DAYS(LDAS%MO)
       LDAS%MO=LDAS%MO+1
      ENDDO

      DO WHILE(LDAS%MO.GT.12)
       LDAS%MO=LDAS%MO-12
       LDAS%YR=LDAS%YR+1
      ENDDO

!=== Update LDAS current model TIME Variable
      CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,
     1  LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS)
      WRITE(*,24)'GSFC-LDAS Time: ',LDAS%MO,'/',LDAS%DA,'/',
     1  LDAS%YR,LDAS%HR,':',LDAS%MN,':',LDAS%SS

      WRITE(79,24)'GSFC-LDAS Time: ',LDAS%MO,'/',LDAS%DA,'/',
     1  LDAS%YR,LDAS%HR,':',LDAS%MN,':',LDAS%SS

 24   FORMAT(A16,I2,A1,I2,A1,I4,1X,I2,A1,I2,A1,I2)

!=== Set HDF output timing variables
      LDAS%YYYYMMDD=(10000*LDAS%YR)+(100*LDAS%MO)+LDAS%DA
      LDAS%HHMMSS=(10000*LDAS%HR)+(100*LDAS%MN)+LDAS%SS

!=== Check for ENDTIME
      IF(LDAS%ENDCODE.EQ.0)THEN  !End at real-time date (TBD)
       WRITE(*,*)'WARNING: Do not know how to stop in real-time' 
       WRITE(79,*)'WARNING: Do not know how to stop in real-time'
      ENDIF

      IF(LDAS%ENDCODE.EQ.1)THEN  !End on date specified in ldas.crd file
       CALL DATE2TIME(LDAS%ETIME,LDAS%EDOY,LDAS%EGMT,
     1   LDAS%EYR,LDAS%EMO,LDAS%EDA,LDAS%EHR,LDAS%EMN,LDAS%ESS)
       IF(LDAS%TIME.GE.LDAS%ETIME)THEN
        LDAS%ENDTIME=1
        WRITE(*,*) 'GSFC-LDAS Run Completed'
        WRITE(79,*) 'GSFC-LDAS Run Completed'
       ENDIF
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
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! tick.f: 
!
! DESCRIPTION:
!  Advance (or retract) time variables a specified amount 
!  (A nonmodular version of ticktime.f.)
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!=========================================================================

      SUBROUTINE TICK(TIME,DOY,GMT,YR,MO,DA,HR,MN,SS,TS)

      IMPLICIT NONE

!=== LOCAL Variables =====================================================
      INTEGER DAYS(13)
      INTEGER YR,MO,DA,HR,MN,SS,TS,DOY
	INTEGER PRVMO   !previous month
      REAL*8 TIME
      REAL GMT
      DATA DAYS/31,28,31,30,31,30,31,31,30,31,30,31,31/

!=== End Variable List ===================================================
143    format(a1,' YR',i6,' MO',i5,' DY',i5,' HR',i5,
     &   ' MN',i6,' SS',i8,' TS',i8)
        SS=SS+TS
C	write(*,143) 'A',YR,MO,DA,HR,MN,SS,TS
      DO WHILE(SS.GT.59)
       SS=SS-60
       MN=MN+1
      ENDDO 
      DO WHILE(SS.LT.0)
       SS=SS+60
       MN=MN-1
      ENDDO 
C	write(*,143) 'B',YR,MO,DA,HR,MN,SS,TS

      DO WHILE(MN.GT.59)
       MN=MN-60
       HR=HR+1
      ENDDO 

      DO WHILE(MN.LT.0)
       MN=MN+60
       HR=HR-1
      ENDDO 
C	write(*,143) 'C',YR,MO,DA,HR,MN,SS,TS

      DO WHILE(HR.GT.23)
       HR=HR-24
       DA=DA+1
      ENDDO 

      DO WHILE(HR.LT.0)
       HR=HR+24
       DA=DA-1
      ENDDO 
C	write(*,143) 'D',YR,MO,DA,HR,MN,SS,TS

      IF((MOD(YR,4).EQ.0.AND.MOD(YR,100).NE.0)           !correct for leap year
     &    .OR.(MOD(YR,400).EQ.0))THEN               !correct for Y2K
       DAYS(2)=29                  
      ELSE
       DAYS(2)=28
      ENDIF

      DO WHILE(DA.GT.DAYS(MO))
       DA=DA-DAYS(MO)
       MO=MO+1
      ENDDO

      DO WHILE(DA.LT.1)
	
	 PRVMO=MO-1
	 if(MO.eq.1) PRVMO=12
	
       DA=DA+DAYS(PRVMO)
	 
	 if(PRVMO.eq.12) then
	   MO=PRVMO
	   YR=YR-1
	 else
	   MO=PRVMO
	 endif
      ENDDO
C	write(*,143) 'E',YR,MO,DA,HR,MN,SS,TS

      DO WHILE(MO.GT.12)
       MO=MO-12
       YR=YR+1
      ENDDO

      DO WHILE(MO.LT.1)
       MO=MO+12
       YR=YR-1
      ENDDO
C	write(*,143) 'F',YR,MO,DA,HR,MN,SS,TS

!=== Update TIME Variable
      CALL DATE2TIME(TIME,DOY,GMT,YR,MO,DA,HR,MN,SS)
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
! time2date.f & date2time: 
!
! DESCRIPTION:
!  Determines time in years, based on year, month, day hour etc..
!   OR reverse (date2time).  
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Code
!  21 Feb 2002: Brian Cosgrove; Corrected leap year code line.  Days(2)
!               was not being reset to 28 after leaving a leap year,
!               it was staying 29
!=========================================================================

      SUBROUTINE DATE2TIME(TIME,DOY,GMT,YR,MO,DA,HR,MN,SS)

      IMPLICIT NONE
      INTEGER YR,MO,DA,HR,MN,SS,YRDAYS,DOY,DAYS(13),K
      REAL*8 TIME
      REAL GMT
      DATA DAYS /31,28,31,30,31,30,31,31,30,31,30,31,30/

      IF((MOD(YR,4).EQ.0.AND.MOD(YR,100).NE.0)    !correct for leap year
     &    .OR.(MOD(YR,400).EQ.0))THEN             !correct for Y2K
       YRDAYS=366                  
      ELSE
       YRDAYS=365
      ENDIF

      DOY=0
      DO K=1,(MO-1)
       DOY=DOY+DAYS(K)
      ENDDO
      DOY=DOY+DA

      IF(YRDAYS.EQ.366.AND.MO.GT.2)DOY=DOY+1

      TIME=(DFLOAT(YR)+((((((DFLOAT(SS)/60.d0)+DFLOAT(MN))/60.d0)+
     1 DFLOAT(HR))/24.d0)+DFLOAT(DOY-1))/DFLOAT(YRDAYS))

      GMT=( ( (FLOAT(SS)/60.0) +FLOAT(MN)) /60.0)+FLOAT(HR)
      RETURN
      END


      SUBROUTINE TIME2DATE(TIME,DOY,GMT,YR,MO,DA,HR,MN)

      IMPLICIT NONE
      INTEGER YR,MO,DA,HR,MN,SS,YRDAYS,DOY,DAYS(13)
      REAL*8 TIME,TMP
      REAL GMT
      DATA DAYS /31,28,31,30,31,30,31,31,30,31,30,31,30/

      YR  = DINT(TIME)
      TMP =     (TIME) 

      IF((MOD(YR,4).EQ.0.AND.MOD(YR,100).NE.0)    !correct for leap year
     &    .OR.(MOD(YR,400).EQ.0))THEN             !correct for Y2K
       YRDAYS=366                  
      ELSE
       YRDAYS=365
      ENDIF
      IF (YRDAYS.EQ.366) THEN
       DAYS(2)=29
      ELSE
       DAYS(2)=28
      ENDIF

      DOY  = DINT((TMP-YR)*DFLOAT(YRDAYS))+1 
      TMP =      ((TMP-YR)*DFLOAT(YRDAYS))+1 
      HR  = NINT((TMP-DOY)*24.D0) 
      TMP =     ((TMP-DOY)*24.D0) 

      MN  = DINT((TMP-HR)*60.D0) 
      TMP =     ((TMP-HR)*60.D0) 

      SS  = DINT((TMP-MN)*60.D0) 
      MO=1
      DO WHILE (DOY.GT.0)
       DOY=DOY-DAYS(MO)
       MO=MO+1
      ENDDO
      MO=MO-1
      DA=DOY+DAYS(MO)

      GMT=(((FLOAT(SS)/60.0)+FLOAT(MN))/60.0)+FLOAT(HR)

c	print *,'yr,mo,da,gmt'
c	print *,yr,mo,da,gmt
  	IF(GMT.EQ.24) THEN
         GMT=0
         DA=DA+1
         IF (DA.GT.DAYS(MO)) THEN
          DA=1
          MO=MO+1
          IF (MO.GT.12) THEN
           MO=1
           YR=YR+1
          ENDIF
         ENDIF
        ENDIF

c        print *,'yr,mo,da,gmt'
c        print *,yr,mo,da,gmt
  
      RETURN
      END

