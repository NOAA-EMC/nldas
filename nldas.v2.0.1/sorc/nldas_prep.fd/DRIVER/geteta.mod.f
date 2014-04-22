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
! geteta.f: 
!
! DESCRIPTION:
!  Opens, reads, and interpolates EDAS-LDAS forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data 
!
!    The idea is to open a 3HR initialization file
!    if that fails a 3HR forecast file
!    if that fails then a 6HR forecast file
!  If these fail then the time must be rolled back
!  And thus we use "FILETIME" instead of "TIME" 
!    So the next step is a forecast file from 12 hours ago
!    that is forecasting for 15 or 18 hours ahead (i.e. 12+3 or 12+6)
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  25 Oct 1999: Jared Entin; Significant F90 Revision
!  18 Feb 2000: Brian Cosgrove; Changes made for use with precip merging
!  24 May 2001: Brian Cosgrove; Fixed error in calculation of precip
!               per second rate.  Now divides by 6*60*60 if ETA 6hr used
!               instead of only dividing by 3*60*60
!  12 Feb 2002: Brian Cosgrove; Changed code to make use of 84 hour ETA 3 hourly forecasts...
!               previously only used out to 48 hrs
!  31 Jul 2003: Brian Cosgrove; Changed local ftype variable to global ftype variable
!               Previous to this, 3/6? hour eta data was not being used correctly in
!               some situations on sgi machines due to local variable being wiped out
!               each time subroutine exited...does not seem to have affected forcing
!               but this is a safeguard to prevent any possible problems.
!=========================================================================
      SUBROUTINE GETETAMOD(LDAS,GRID,LYR,LMO,LDA,LHR,LMN,try)
	
	USE ldas_module     ! LDAS non-model-specific 1-D variables
	USE grid_module     ! LDAS non-model-specific grid variables
	IMPLICIT NONE
	type (ldasdec) LDAS
	type (griddec) GRID(LDAS%NC,LDAS%NR)
	
	
!==== Local Variables=======================
       INTEGER C,R,F,FERROR,TRY,ZDOY
	 INTEGER YR1,MO1,DA1,HR1,MN1,SS1,DOY1,TS1
	 INTEGER YR2,MO2,DA2,HR2,MN2,SS2,DOY2,TS2
	 INTEGER LYR,LMO,LDA,LHR,LMN
	 INTEGER YR26,MO26,DA26,HR26,MN26,SS26,DOY26,TS26
	 INTEGER TS1H   !to go back twelve hours
	 INTEGER TS2D   !to go forward two days (restoring the role back)
	 INTEGER TS1Day   !to go backwards 1 Day
	 REAL*8 TIME1,TIME2,FAKETIME,TIME26
	 REAL*8 DUMBTIME1,DUMBTIME2,DUMBTIME26
	 REAL*8 TIMENOW
	 REAL*8 FILETIME1,FILETIME2    !these are the time that the file correspond too
	 CHARACTER*80 NAME
	 CHARACTER*80 Prevname  !previous time period's file name
	 INTEGER Precflag       !do you need Prevname to adjust precip
	 REAL WT1,WT2,ZW1,ZW2,CZB,CZE,GMT1,GMT2,GMT26
	 INTEGER LFB      !this is the loop for looking backward or forward in time
	 INTEGER FCA      !Forecast Amount, 1=12 hrs, 2=24hrs, etc
	 INTEGER dataflag   !1=edas/eta, 0=ncep(hourly)
c	 INTEGER ftype      !file type 1=EDAS, 2=ETA3hr, 3=ETA6hr
	 INTEGER ORDER    !1=time before, 2=time after
	 
	 INTEGER FINDTIME1     ! 0=don't get new file for 1st time (or 2nd)
	 INTEGER FINDTIME2     ! 1=Get a new file 1st time (or 2nd)
	 INTEGER MOVETIME      ! 1=Move Time 2 data into Time 1
	 
	 dataflag=1   !for Eta data of some sort (i.e. EDAS, 3,6hr forecasts)
	 
C   Assumption will be not to find or move any data
         FINDTIME1=0
	 FINDTIME2=0
	 MOVETIME=0	 
	 
	 
!=== End Variable Definition =======================	 
	
!=== Determine Required EDAS Data Times (The previous hour and the future hour)
!=== The Adjustment of the Hour and the direction will be done
!=== in the subroutines that generate the names b/c it's different
!=== for the three or six hour time steps 
      YR1=LYR    !Previous Hour
      MO1=LMO
      DA1=LDA
      HR1=LHR
      MN1=LMN
      SS1=0
      TS1=0
      CALL TICK(TIMENOW,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)
	YR1=LYR    !Previous Hour
      MO1=LMO
      DA1=LDA
      HR1=3*((LHR)/3)
      MN1=0
      SS1=0
      TS1=0
      CALL TICK(TIME1,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)	
      YR2=LYR    !Next Hour
      MO2=LMO
      DA2=LDA
      HR2=3*((LHR)/3)
	MN2=0
	SS2=0
	TS2=3*60*60
	CALL TICK(TIME2,DOY2,GMT2,YR2,MO2,DA2,HR2,MN2,SS2,TS2)	 
      YR26=LYR
	MO26=LMO
	DA26=LDA
	HR26=6*((LHR)/6)
	MN26=0
	SS26=0
	TS26=6*60*60
	CALL TICK(TIME26,DOY26,GMT26,
     &	    YR26,MO26,DA26,HR26,MN26,SS26,TS26)
      DUMBTIME1=TIME1
	DUMBTIME2=TIME2
	DUMBTIME26=TIME26
	
      if(TIMENOW.GT.LDAS%PRECIPTIME2) then
           MOVETIME=1
           FINDTIME2=1
        endif
        if(TIMENOW.LT.LDAS%PRECIPTIME1) then   !beginning of the run
           FINDTIME1=1
           FINDTIME2=1
        ENDIF
	if((LDAS%PRECIPTIME1.EQ.0).and.(LDAS%PRECIPTIME2.EQ.0)) then
           FINDTIME1=1
           FINDTIME2=1
        endif


20    format('MOVE',i2,2x,'FIND1',i2,2x,'FIND2',i2)
	 
!===  Check to see if required data in memory, if not, read it
	  IF(FINDTIME1.eq.1) THEN !We need to get new TIME1 data
	      ORDER=1        !says that it's the "before" time
            FERROR=0
		LFB=1
		TS1=-60*60*24
		TS1H=-60*60*12
		DO WHILE (FERROR.eq.0)
		  TRY=TRY+1
		  CALL ETAEDASPRECIPFILE(NAME,LDAS%EDASDIR,YR1,MO1,DA1,HR1,LFB)
		  Precflag=0   !becuase EDAS only accumulates every 3 hours not longer
		  ldas%ftype=1   !says that the data will be EDAS
		  CALL RETETAMOD(ORDER,LDAS,GRID,NAME,FERROR,
     &  ldas%ftype,
     &           dataflag,prevname,precflag,try)

 		  FCA=0

C  EDAS failed go to ETA 3HR or 6HR Forecast files		  
		  IF(FERROR.eq.0) then
		    do while (FERROR.eq.0.and.FCA.le.6)
		     CALL ETA3HRFILE(NAME,LDAS%ETA3HRDIR,YR1,MO1,DA1,HR1,
     &		     LFB,FCA,Prevname,Precflag)
		  ldas%ftype=2   !says that the data will be ETA 3 hr
		     CALL RETETAMOD(ORDER,LDAS,GRID,NAME,FERROR,
     &  ldas%ftype,
     &           dataflag,prevname,precflag,try)
		     if(FERROR.eq.0) then
		       CALL ETA6HRFILE(NAME,LDAS%ETA6HRDIR,
     &		 YR1,MO1,DA1,HR1,LFB,FCA,Prevname,Precflag)
		  ldas%ftype=3   !says that the data will be 6 hr ETA
c		  if(Precflag.eq.1) print*,'brian1',Prevname
			 CALL RETETAMOD(ORDER,LDAS,GRID,NAME,FERROR,
     &  ldas%ftype,
     &           dataflag,prevname,precflag,try)
		     endif
		     if(FERROR.eq.0) then
		      FCA=FCA+1
			CALL Tick(DUMBTIME1,DOY1,GMT1,
     &		       YR1,MO1,DA1,HR1,MN1,SS1,TS1H)		     
		     endif
 		     
		    enddo  !the do while 2 conditions
		    TS2D=TS1H*(-1)*FCA
		    Call Tick(DUMBTIME1,DOY1,GMT1,
     &                 YR1,MO1,DA1,HR1,MN1,SS1,TS2D)
		  Endif  !the if Ferror=0
			
		  IF(FERROR.eq.1) then
		    LDAS%PRECIPTIME1=TIME1
		  ELSE
		    CALL TICK(DUMBTIME1,DOY1,GMT1,
     &                 YR1,MO1,DA1,HR1,MN1,SS1,TS1)
		    FCA=0
		  ENDIF
c  Limit the number of back looking to ten days		  
		  IF(TRY.GT.100) THEN
		    WRITE(*,*) 'ERROR: ETA Data gap exceeds 10 days on file 1'
		    STOP
		  ENDIF
c		  
		ENDDO   !the while Ferror=0
		
		
C  You need to make sure that the second loop is going
C  to be set off here too		
		
	   ENDIF    !time1.ne.LDASPRECIPTIME1
C
!  REPEAT FOR TIME 2

	   IF(MOVETIME.eq.1) THEN !Transfer TIME2 data to TIME1
	     LDAS%PRECIPTIME1=LDAS%PRECIPTIME2	
	     FINDTIME2=1
	     do F=1,LDAS%NF
	      do C=1,LDAS%NC
		 do R=1,LDAS%NR
		   GRID(C,R)%PRECIPDATA1(1)=GRID(C,R)%PRECIPDATA2(1)
		 enddo
		enddo
	     enddo
C
         ENDIF  ! if MOVETIME=1



C          Beginning point
	  IF(FINDTIME2.eq.1) THEN !We need to get new TIME2 data
           ORDER=2
	     FERROR=0
		LFB=1
		TS1=-60*60*24
		TS1H=-60*60*12
		DO WHILE (FERROR.eq.0)
		  TRY=try+1
		  CALL ETAEDASPRECIPFILE(NAME,LDAS%EDASDIR,YR2,MO2,DA2,HR2,LFB)
		  Precflag=0   !becuase EDAS only accumulates every 3 hours not longer
		  ldas%ftype=1   !data will be EDAS
		  CALL RETETAMOD(ORDER,LDAS,GRID,NAME,FERROR,
     &  ldas%ftype,
     &           dataflag,prevname,precflag,try)
		  if(FERROR.eq.1) LDAS%PRECIPTIME2=TIME2
 		  FCA=0

C  EDAS failed, go to ETA 3HR or 6HR Forecast files		  
		  IF(FERROR.eq.0) then
		    do while (FERROR.eq.0.and.FCA.le.6)
		     CALL ETA3HRFILE(NAME,LDAS%ETA3HRDIR,YR2,MO2,DA2,HR2,
     &		     LFB,FCA,Prevname,Precflag)
		  ldas%ftype=2   !data will be ETA 3hr
c		  if(Precflag.eq.1) print*,Prevname
		     CALL RETETAMOD(ORDER,LDAS,GRID,NAME,FERROR,
     &  ldas%ftype,
     &           dataflag,prevname,precflag,try)
		     if(FERROR.eq.1) LDAS%PRECIPTIME2=TIME2
		     if(FERROR.eq.0) then
		       CALL ETA6HRFILE(NAME,LDAS%ETA6HRDIR,
     &		 YR26,MO26,DA26,HR26,LFB,FCA,Prevname,Precflag)
		       ldas%ftype=3    !data will be ETA 6HR
                     CALL RETETAMOD(ORDER,LDAS,GRID,NAME,FERROR,
     &  ldas%ftype,
     &            dataflag,prevname,precflag,try)
		       if(FERROR.eq.1) LDAS%PRECIPTIME2=TIME26
     
		      endif  ! if FERROR=0
		        if(FERROR.eq.0) then   !Still haven't gotten good day, role back forecast
		        FCA=FCA+1
C  TS1H=-60*60*12
			  CALL Tick(DUMBTIME2,DOY2,GMT2,
     &			 YR2,MO2,DA2,HR2,MN2,SS2,TS1H)
			  CALL Tick(DUMBTIME26,DOY26,GMT26,
     &                   YR26,MO26,DA26,HR26,MN26,SS26,TS1H)					     
		       endif
		    enddo  !the do while 2 conditions
		    TS2D=TS1H*(-1)*FCA          ! Roll back 12HRs * FCA
C  This is to restore the hours subtracted for forecasts
C  Back to the original time so that we can then roll back a day		
		    Call Tick(DUMBTIME2,DOY2,GMT2,
     &                 YR2,MO2,DA2,HR2,MN2,SS2,TS2D)
		    CALL Tick(DUMBTIME26,DOY26,GMT26,
     &                 YR26,MO26,DA26,HR26,MN26,SS26,TS2D)	
		  Endif  !the if Ferror=0
			
		  IF(FERROR.eq.0) then
		    TS1Day=-24*60*60   ! minus one Day
		    CALL TICK(DUMBTIME2,DOY2,GMT2,
     &		     YR2,MO2,DA2,HR2,MN2,SS2,TS1Day)
		    CALL TICK(DUMBTIME26,DOY26,GMT26,
     &                 YR26,MO26,DA26,HR26,MN26,SS26,TS1Day)		    
		    FCA=0
		  ENDIF
c  Limit the number of back looking to ten days		  
		  IF(TRY.GT.100) THEN
		    WRITE(*,*) 'ERROR: ETA Data gap exceeds 10 days on file 2'
		    STOP
		  ENDIF
c		  
		ENDDO   !the while Ferror=0
		
	   ENDIF    ! IF FINDTIME=1

	   
	   
!===  Interpolate Data in Time
	   DO C=1,LDAS%NC
	    DO R=1,LDAS%NR
	     GRID(C,R)%PRECIP(1)=GRID(C,R)%PRECIPDATA2(1)
	    ENDDO
	   ENDDO

	   
!=== ADJUST PRECIP TO VALUE PER SEC DATA SHOULD COME IN AS VALUE PER
!=== TIME PERIOD THREE HOUR FOR ETA	 
!--- POSSIBLE ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR  ERROR ERROR ERROR ERROR  ERROR ERROR ERROR ERROR 

         IF((ldas%ftype.eq.1).or.(ldas%ftype.eq.2)) then !We got ETA 3hr or EDAS data
         DO C=1,LDAS%NC
            DO R=1,LDAS%NR
            GRID(C,R)%PRECIP(1)=GRID(C,R)%PRECIP(1)/(3.0*60.0*60.0)
            ENDDO
           ENDDO
         ENDIF

         IF (ldas%ftype.eq.3) then !We got ETA 6hr data
         DO C=1,LDAS%NC
            DO R=1,LDAS%NR
            GRID(C,R)%PRECIP(1)=GRID(C,R)%PRECIP(1)/(6.0*60.0*60.0)
            ENDDO
           ENDDO
         ENDIF

	   RETURN 
	   END
	   

      
       
