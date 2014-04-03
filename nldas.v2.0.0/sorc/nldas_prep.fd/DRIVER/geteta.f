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
!  11 Apr 2000: Brian Cosgrove; Fixed name construction error 
!               in Subroutine ETA6HRFILE 
!  27 Apr 2000: Brian Cosgrove; Added correction for use of old shortwave
!               data with opposite sign convention from recent shortwave data.
!               Added capability to use time averaged shortwave and longwave data.
!               Altered times which are passed into ZTERP--used to be GMT1 and GMT2,
!               now they are LDAS%ETATIME1 and LDAS%ETATIME2
!  11 May 2000: Brian Cosgrove; Added checks for SW values that are too high
!               due to zenith angle weighting...if too high, use linear weighting.
!               Also, if cos(zen) less than .1, then use linear weighting to
!               avoid computed values of SW that are too high.
!  18 May 2000: Brian Cosgrove; Corrected line of code in ETAEDASNAME which
!               assigned wrong year directory variable when constructing
!               EDAS filename
!  26 May 2000: Jared Entin; Changed numerical bound of the TRY variable
!               to fix a rollback problem.
!  5 June 2000: Brian Cosgrove; Fixed a problem with the correction of the negative
!               radiation sign convention. Prior to fix, was not correcting negative
!               values of -999.9...now it changes all negative values to positive ones.
!  18 Aug 2000: Brian Cosgrove; Fixed undefined value problem in check for negative
!               radiation values over land points.
!  27 Feb 2001: Brian Cosgrove; Added CZM into call for ZTERP subroutine
!  24 May 2001: Brian Cosgrove; Fixed error in calculation of precip
!               per second rate.  Now divides by 6*60*60 if ETA 6hr used
!               instead of only dividing by 3*60*60
!   4 Sep 2001: Brian Cosgrove; Added range checks for output variables.  Added
!               fix for specific humidity...if using EDAS forcing, then correct
!               specific humidity using 0-30mb specific humidity and 2 meter
!               pressure and temperature .  Correct variables to ALMA/LDAS
!               standards.  Keep ETA SW as seperate variable.
!   5 Feb 2002: Brian Cosgrove; Changed local ftype variable to global ftype variable
!               Previous to this, 3/6? hour eta data was not being used correctly in
!               some situations 
!  12 Feb 2002: Brian Cosgrove; Changed code to make use of 84 hour ETA 3 hourly forecasts...
!               previously only used out to 48 hrs
!  1 Dec 2006:  Brian Cosgrove; Changed code to use NARR instead of EDAS
!  11 Dec 2006: Brian Cosgrove; Altered code so that LW and SW data in second memory position
!               is used every three hours since NARR radiation fields are stored
!               backwards compared to EDAS....that is, a 00Z NARR files contains
!               radiation valid from 00Z to 03Z.
!=========================================================================
      SUBROUTINE GETETA(LDAS,GRID)
	
	use ldas_module     ! LDAS non-model-specific 1-D variables
	use grid_module     ! LDAS non-model-specific grid variables
	implicit none
	type (ldasdec) ldas
	type (griddec) grid(ldas%nc,ldas%nr)
	
	
!==== Local Variables=======================
       INTEGER  ii   
       INTEGER C,R,F,FERROR,TRY,ZDOY,cc,rr
	 INTEGER YR1,MO1,DA1,HR1,MN1,SS1,DOY1,TS1,BDOY,BYR,BMO
	 INTEGER YR2,MO2,DA2,HR2,MN2,SS2,DOY2,TS2,BDA,BHR,BMN
	 INTEGER YR26,MO26,DA26,HR26,MN26,SS26,DOY26,TS26,check
	 INTEGER zzHR1,zzMN1  !for diagnostics
	 INTEGER TS1H   !to go back twelve hours
	 INTEGER TS2D   !to go forward two days (restoring the role back)
	 INTEGER TS1Day   !to go backwards 1 Day
	 REAL*8 TIME1,TIME2,FAKETIME,TIME26
	 REAL*8 DUMBTIME1,DUMBTIME2,DUMBTIME26
	 REAL*8 TIMENOW,BTIME
	 REAL*8 FILETIME1,FILETIME2    !these are the time that the file correspond too
	 CHARACTER*80 NAME
	 CHARACTER*80 Prevname  !previous time period's file name
	 INTEGER Precflag       !do you need Prevname to adjust precip
	 REAL WT1,WT2,ZW1,ZW2,CZB,CZE,GMT1,GMT2,GMT26,CZM,czmean
	 REAL CRATIO,outradmask(464,224)
	 INTEGER LFB      !this is the loop for looking backward or forward in time
	 INTEGER FCA      !Forecast Amount, 1=12 hrs, 2=24hrs, etc
	 INTEGER dataflag   !1=EDAS/eta, 0=ncep(hourly)
c	 INTEGER ftype      !file type 1=EDAS, 2=ETA3hr, 3=ETA6hr
         INTEGER shftype1 !spec humid file type 1=EDAS, 2=ETA3 3=ETA6 for past data
	 INTEGER shftype2 !spec humid file type 1=EDAS, 2=ETA3, 3=ETA6 for future data 
	 INTEGER ORDER    !1=time before, 2=time after
	 
	 INTEGER FINDTIME1     ! 0=don't get new file for 1st time (or 2nd)
	 INTEGER FINDTIME2     ! 1=Get a new file 1st time (or 2nd)
	 INTEGER MOVETIME      ! 1=Move Time 2 data into Time 1
	 INTEGER SUAVE
	 real tempgmt1,tempdumb2,tempgmt2,tempyr2,tempmo2,tempda2
         real tempmn2,tempss2,forwardtime,tempdoy2,temphr2
	 real rho,newtemp,vap,RH,W,Ws,E,Es
	 SUAVE=0
	 
	 dataflag=1   !for Eta data of some sort (i.e. EDAS, 3,6hr forecasts)
	 
C   Assumption will be not to find or move any data
       FINDTIME1=0
	 FINDTIME2=0
	 MOVETIME=0	 

	check=0	 
!=== End Variable Definition =======================	 

	
!=== Determine Required EDAS Data Times (The previous hour and the future hour)
!=== The Adjustment of the Hour and the direction will be done
!=== in the subroutines that generate the names b/c it's different
!=== for the three or six hour time steps 
      YR1=LDAS%YR    !Previous Hour
      MO1=LDAS%MO
      DA1=LDAS%DA
      HR1=LDAS%HR
	zzHR1=HR1
      MN1=LDAS%MN
	zzMN1=MN1
      SS1=0
      TS1=0
      CALL TICK(TIMENOW,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)
	YR1=LDAS%YR    !Previous Hour
      MO1=LDAS%MO
      DA1=LDAS%DA
      HR1=3*((LDAS%HR)/3)
      MN1=0
      SS1=0
      TS1=0
      CALL TICK(TIME1,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)


      YR2=LDAS%YR    !Next Hour
      MO2=LDAS%MO
      DA2=LDAS%DA
      HR2=3*((LDAS%HR)/3)
	MN2=0
	SS2=0
	TS2=3*60*60
	CALL TICK(TIME2,DOY2,GMT2,YR2,MO2,DA2,HR2,MN2,SS2,TS2)	 



      YR26=LDAS%YR
	MO26=LDAS%MO
	DA26=LDAS%DA
	HR26=6*((LDAS%HR)/6)
	MN26=0
	SS26=0
	TS26=6*60*60
	CALL TICK(TIME26,DOY26,GMT26,
     &	    YR26,MO26,DA26,HR26,MN26,SS26,TS26)
     
      DUMBTIME1=TIME1
	DUMBTIME2=TIME2
	DUMBTIME26=TIME26

        if(TIMENOW.GT.LDAS%ETATIME2) then
	   MOVETIME=1
	   FINDTIME2=1
	endif

	if(LDAS%TSCOUNT.EQ.1) then     !beginning of the run
	   FINDTIME1=1
	   FINDTIME2=1
           MOVETIME=0
	endif
      
20    format('MOVE',i2,2x,'FIND1',i2,2x,'FIND2',i2)
    

!==== If FINDTIME2=1 then you want initialize the FSOURCE info
         if(FINDTIME2.eq.1) then
	     LDAS%FSOURCE(1)=0     !EDAS Status
	     LDAS%FSOURCE(2)=0     !ETA 3hr day status
	     LDAS%FSOURCE(3)=0     !ETA 3hr forecast amount
	     LDAS%FSOURCE(4)=0     !ETA 6hr day status
	     LDAS%FSOURCE(5)=0     !ETA 6hr forecast amount
	   endif

	 
!===  Check to see if required data in memory, if not, read it
	  IF(FINDTIME1.eq.1) THEN !We need to get new TIME1 data
	      ORDER=1        !says that it's the "before" time
            FERROR=0
		TRY=0
		LFB=1
		TS1=-60*60*24
		TS1H=-60*60*12
		DO WHILE (FERROR.eq.0)
		  TRY=TRY+1
		  CALL ETAEDASFILE(NAME,LDAS%EDASDIR,YR1,MO1,DA1,HR1,
     &            LFB)
		  Precflag=0   !becuase EDAS only accumulates every 3 hours not longer
		  ldas%ftype=1   !says that the data will be EDAS

		  CALL RETETA(ORDER,LDAS,GRID,NAME,FERROR,ldas%ftype,
     &           dataflag,prevname,precflag,try)
c	print* ,'FIRST EDAS ferror= ',ferror
 		  FCA=0

C  EDAS failed go to ETA 3HR or 6HR Forecast files		  
		  IF(FERROR.eq.0) then
		    do while (FERROR.eq.0.and.FCA.le.6)
		     
		     CALL ETA3HRFILE(NAME,LDAS%ETA3HRDIR,YR1,MO1,
     &               DA1,HR1,
     &		     LFB,FCA,Prevname,Precflag)
c		  print*,'PRINTING ETA 3hr NAME, next line'
C		  print*,NAME
		  ldas%ftype=2   !says that the data will be ETA 3 hr
c		  if(Precflag.eq.1) print*,Prevname
c		  print*,'Done with ETA 3hr Name'
		     CALL RETETA(ORDER,LDAS,GRID,NAME,FERROR,ldas%ftype,
     &           dataflag,prevname,precflag,try)
c		             print* ,'eta3 ferror= ',ferror
	
		     if(FERROR.eq.0) then
		       CALL ETA6HRFILE(NAME,LDAS%ETA6HRDIR,
     &		 YR1,MO1,DA1,HR1,LFB,FCA,Prevname,Precflag)
c		  print*,'PRINTING ETA 6hr NAME, next line'
C		  print*,NAME
		  ldas%ftype=3   !says that the data will be 6 hr ETA
c		  if(Precflag.eq.1) print*,'brian 1',Prevname
c		  print*,'Done with ETA 6hr Name'     
			 CALL RETETA(ORDER,LDAS,GRID,NAME,FERROR,
     &                   ldas%ftype,
     &           dataflag,prevname,precflag,try)
c        print* ,'eta6 ferror= ',ferror

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
		    LDAS%ETATIME1=TIME1
		  ELSE
c		    print*,'Roll back one day and start again'
		    CALL TICK(DUMBTIME1,DOY1,GMT1,
     &                 YR1,MO1,DA1,HR1,MN1,SS1,TS1)
c		    print*,'Forecast amount was',FCA
		    FCA=0
		  ENDIF
c  Limit the number of back looking to ten days		  
		  IF(TRY.GT.10) THEN    !1) THEN
             WRITE(*,*) 'ERROR: ETA Data gap exceeds 
     &  10 days on file 1'
		    STOP
		  ENDIF
c		  
		ENDDO   !the while Ferror=0

!c	print *,'brian first NOW, ldas%ftype=',ldas%ftype
!	IF (ldas%ftype.EQ.3) THEN
!!c	print *,'yr=',ldas%yr,'month=',ldas%mo
!        IF ((LDAS%YR.EQ.1998).AND.(LDAS%MO.GT.1).AND.
!     &  (LDAS%MO.LT.12)) THEN
!c	print *,'CALLING ETA6SPFH FIRST TIME',NAME,PREVNAME
!	  CALL RETETASPFH(ORDER,LDAS,GRID,NAME,FERROR,
!     &                   ldas%ftype,
!     &           dataflag,prevname,precflag,try)
!	ENDIF
!	ENDIF

		
		
C  You need to make sure that the second loop is going
C  to be set off here too		
		
	   ENDIF    !time1.ne.LDASETATIME1
	shftype1=ldas%ftype
C
!  REPEAT FOR TIME 2

	   IF(MOVETIME.eq.1) THEN !Transfer TIME2 data to TIME1
	     LDAS%ETATIME1=LDAS%ETATIME2	
	     FINDTIME2=1 
	     do F=1,LDAS%NF
	      do C=1,LDAS%NC
		 do R=1,LDAS%NR
		   GRID(C,R)%ETADATA1(F)=GRID(C,R)%ETADATA2(F)
		 enddo
		enddo
	     enddo
C	      print*,'done with loop'
           LDAS%EVT1=LDAS%EVT2
C
         ENDIF  ! if MOVETIME=1



C          Beginning point
	  IF(FINDTIME2.eq.1) THEN !We need to get new TIME2 data
c	     print*,'beginning Finding Time 2'
           ORDER=2
	     FERROR=0
C            print*,'try',TRY
		TRY=0
C		print*,'TRY',TRY,'  LFB',LFB
		LFB=1
		TS1=-60*60*24
		TS1H=-60*60*12
C            print*,'TS1',TS1,'TS1H',TS1H,'  LFB',LFB
C		print*,'FERROR',FERROR
		
		DO WHILE (FERROR.eq.0)
c		  print*,'TRY is ',TRY
		  TRY=TRY+1
		  
		  if(TRY.eq.1) then
		    LDAS%FSOURCE(1)=1
		    LDAS%FSOURCE(2)=1
		    LDAS%FSOURCE(4)=1		  
		  else
		    LDAS%FSOURCE(1)=(-1*TRY)+1
		    LDAS%FSOURCE(2)=(-1*TRY)+1
		    LDAS%FSOURCE(4)=(-1*TRY)+1
		  endif
		  
c		  print*,'Calling EDAS Namer'
		  CALL ETAEDASFILE(NAME,LDAS%EDASDIR,YR2,MO2,DA2,
     &            HR2,LFB)
		  Precflag=0   !becuase EDAS only accumulates every 3 hours not longer
c		  print*,'PRINTING EDAS NAME, next line'
c		  print*,NAME
		  ldas%ftype=1   !data will be EDAS
c		  print*,'Done with EDAS Name, Calling Retrive'
		  CALL RETETA(ORDER,LDAS,GRID,NAME,FERROR,ldas%ftype,
     &           dataflag,prevname,precflag,try)
c              print*,'Done with Retrieve of EDAS, FERROR=',FERROR
		  if(FERROR.eq.1) LDAS%ETATIME2=TIME2
c		  if(FERROR.eq.1) print*,'EDAS FILE READ in and APPROVED'
 		  FCA=0

C  EDAS failed, go to ETA 3HR or 6HR Forecast files		  
		  IF(FERROR.eq.0) then
		    do while (FERROR.eq.0.and.FCA.le.6)

		     CALL ETA3HRFILE(NAME,LDAS%ETA3HRDIR,YR2,
     &               MO2,DA2,HR2,
     &		     LFB,FCA,Prevname,Precflag)

		  ldas%ftype=2   !data will be ETA 3hr
c		  if(Precflag.eq.1) print*,'brian',Prevname
c
		     CALL RETETA(ORDER,LDAS,GRID,NAME,FERROR,ldas%ftype,
     &           dataflag,prevname,precflag,try)
c  
		     if(FERROR.eq.1) LDAS%ETATIME2=TIME2    !ETA 3hr READ IN and Approved
c	
		     if(FERROR.eq.0) then
		       CALL ETA6HRFILE(NAME,LDAS%ETA6HRDIR,
     &		 YR26,MO26,DA26,HR26,LFB,FCA,Prevname,Precflag)
c	
c		
		       ldas%ftype=3    !data will be ETA 6HR
c		       if(Precflag.eq.1) print*,Prevname
                     CALL RETETA(ORDER,LDAS,GRID,NAME,FERROR,ldas%ftype,
     &            dataflag,prevname,precflag,try)
		       if(FERROR.eq.1) LDAS%ETATIME2=TIME26
     
		      endif  ! if FERROR=0
c	
		        if(FERROR.eq.0) then   !Still haven't gotten good day, role back forecast
		        FCA=FCA+1
C  TS1H=-60*60*12
                    LDAS%FSOURCE(3)=FCA
			  LDAS%FSOURCE(5)=FCA
			  CALL Tick(DUMBTIME2,DOY2,GMT2,
     &			 YR2,MO2,DA2,HR2,MN2,SS2,TS1H)
			  CALL Tick(DUMBTIME26,DOY26,GMT26,
     &                   YR26,MO26,DA26,HR26,MN26,SS26,TS1H)					     
		       endif
c 		      print*,'END While BETA, FCA wuz',FCA
		    enddo  !the do while 2 conditions
		    
		    TS2D=TS1H*(-1)*FCA          ! Roll back 12HRs * FCA
C  This is to restore the hours subtracted for forecasts
C  Back to the original time so that we can then roll back a day		
c                print*,'BTW TS2D is',TS2D    
		    Call Tick(DUMBTIME2,DOY2,GMT2,
     &                 YR2,MO2,DA2,HR2,MN2,SS2,TS2D)
		    CALL Tick(DUMBTIME26,DOY26,GMT26,
     &                 YR26,MO26,DA26,HR26,MN26,SS26,TS2D)	
		  Endif  !the if Ferror=0
			
		  IF(FERROR.eq.0) then
c		    print*,'Roll back one day and start again'
		    TS1Day=-24*60*60   ! minus one Day
		    CALL TICK(DUMBTIME2,DOY2,GMT2,
     &		     YR2,MO2,DA2,HR2,MN2,SS2,TS1Day)
		    CALL TICK(DUMBTIME26,DOY26,GMT26,
     &                 YR26,MO26,DA26,HR26,MN26,SS26,TS1Day)		    
c		    print*,'Forecast amount was',FCA
		    FCA=0
		  ENDIF
c  Limit the number of back looking to ten days		  
		  IF(TRY.GT.10) THEN
	     WRITE(*,*) 'ERROR: ETA Data gap exceeds 
     &  10 days on file 2'
		    STOP
		  ENDIF
c		  
		ENDDO   !the while Ferror=0

        shftype2=ldas%ftype
!        IF (ldas%ftype.EQ.3) THEN
!c        print *,'yr=',ldas%yr,'month=',ldas%mo
!        IF ((LDAS%YR.EQ.1998).AND.(LDAS%MO.GT.1).AND.
!     &  (LDAS%MO.LT.12)) THEN
!          CALL RETETASPFH(ORDER,LDAS,GRID,NAME,FERROR,ldas%ftype,
!     &            dataflag,prevname,precflag,try)
!        ENDIF
!	ENDIF
       
		if(ldas%ftype.eq.1) then !We got EDAS data
		 do ii=2,5
		  LDAS%FSOURCE(ii)=0
		 enddo
		endif
		if(ldas%ftype.eq.2) then !We got ETA 3hr data
		 LDAS%FSOURCE(1)=0
		 LDAS%FSOURCE(4)=0
		 LDAS%FSOURCE(5)=0
		endif
		if(ldas%ftype.eq.3) then !We got ETA 6hr data
		 LDAS%FSOURCE(1)=0
		 LDAS%FSOURCE(2)=0
		 LDAS%FSOURCE(3)=0
		endif
                
        if(LDAS%TSCOUNT.EQ.1) then     !beginning of the run

        DO F=1,LDAS%NF
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
!UGRD VGRD PRESN APCPN ACPCP CAPE TMP PRES HGT SPFH UGRD VGRD HPBL TMP SPFH FRICV SFEXC CRAIN PEVAP DSWRF DLWRF
!      data nineoctet/33,34,134,202,63,157,11,1,7,51,33,34,
!     &  221,11,51,253,208,140,228,204,205/

           IF ((F.EQ.7).or.(F.eq.14)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF ((GRID(C,R)%ETADATA1(F).GT.330.0).OR.
     &          (GRID(C,R)%ETADATA1(F).LT.213.0)) THEN
              newtemp=SQRT(SQRT(ABS(GRID(C,R)%ETADATA1(21))/5.67E-08))
              write (*,'(A19,F10.5,A6,I3,A1,I3,A16,F10.5)')
     &              'TEMP OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &           ') CORRECTING TO ',newtemp
              GRID(C,R)%ETADATA1(F)=newtemp
             ENDIF
            ENDIF
           ENDIF
           
           IF ((F.EQ.10)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).LE.0.0000001) THEN
              vap=(1.0/100.0) * 6.112 *
     &         10**((7.5*(GRID(C,R)%ETADATA1(7)-273.15))/
     &         (237.7+(GRID(C,R)%ETADATA1(7)-273.15)))
              rho=0.62197 * ( vap / ((GRID(C,R)%ETADATA1(8)/100.0)
     &         + vap*(.062197-1) ) )
              write (*,'(A19,F11.8,A6,I3,A1,I3,A16,F11.8)')
     &              'SPFH OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &              ') CORRECTING TO ',rho
              GRID(C,R)%ETADATA1(F)=rho
             ENDIF
             IF (GRID(C,R)%ETADATA1(F).GT.0.03) THEN
              write (*,'(A19,F11.8,A6,I3,A1,I3,A20)')
     &          'SPFH OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 0.03'
              GRID(C,R)%ETADATA1(F)=0.03
             ENDIF
            ENDIF
           ENDIF

           IF ((F.EQ.15)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).LE.0.0000001) THEN
              vap=(1.0/100.0) * 6.112 *
     &         10**((7.5*(GRID(C,R)%ETADATA1(14)-273.15))/
     &         (237.7+(GRID(C,R)%ETADATA1(14)-273.15)))
              rho=0.62197 * ( vap / ((GRID(C,R)%ETADATA1(3)/100.0)
     &         + vap*(.062197-1) ) )
              write (*,'(A19,F11.8,A6,I3,A1,I3,A16,F11.8)')
     &              'SPFH OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &              ') CORRECTING TO ',rho
              GRID(C,R)%ETADATA1(F)=rho
             ENDIF
             IF (GRID(C,R)%ETADATA1(F).GT.0.03) THEN
              write (*,'(A19,F11.8,A6,I3,A1,I3,A20)')
     &          'SPFH OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 0.03'
              GRID(C,R)%ETADATA1(F)=0.03
             ENDIF
            ENDIF
           ENDIF

           IF (F.EQ.20) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).GT.1367) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'EDSWRF OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 1367'
              GRID(C,R)%ETADATA1(F)=1367
             ENDIF
            ENDIF
           ENDIF

           IF (F.EQ.21) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).GT.750) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'EDLWRF OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 750'
              GRID(C,R)%ETADATA1(F)=750
             ENDIF
            ENDIF
           ENDIF
           
           IF ((F.EQ.1).or.(F.eq.11)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).LT.-75) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'UGRD OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &              ') CORRECTING TO -75'
              GRID(C,R)%ETADATA1(F)=-75
             ENDIF

             IF (GRID(C,R)%ETADATA1(F).GT.75) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'UGRD OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 75'
              GRID(C,R)%ETADATA1(F)=75
             ENDIF

            ENDIF
           ENDIF
           IF ((F.EQ.2).or.(F.eq.12)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).LT.-75) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'VGRD OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &              ') CORRECTING TO -75'
              GRID(C,R)%ETADATA1(F)=-75
             ENDIF

             IF (GRID(C,R)%ETADATA1(F).GT.75) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'VGRD OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 75'
              GRID(C,R)%ETADATA1(F)=75
             ENDIF

            ENDIF
           ENDIF

           IF ((F.EQ.8).or.(f.eq.3)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).LT.5000) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'PRES OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &              ') CORRECTING TO 5000'
              GRID(C,R)%ETADATA1(F)=5000
             ENDIF

             IF (GRID(C,R)%ETADATA1(F).GT.110000) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'PRES OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 110000'
              GRID(C,R)%ETADATA1(F)=110000
             ENDIF

            ENDIF
           ENDIF

           IF (F.EQ.4) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).LT.0.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'EAPCP OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &              ') CORRECTING TO 0'
              GRID(C,R)%ETADATA1(F)=0.0
             ENDIF

             IF (GRID(C,R)%ETADATA1(F).GT.140.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'EAPCP OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 140'
              GRID(C,R)%ETADATA1(F)=140.0
             ENDIF

            ENDIF
           ENDIF

           IF (F.EQ.5) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%ETADATA1(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%ETADATA1(F).LT.0.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'EACPCP OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &              ') CORRECTING TO 0'
              GRID(C,R)%ETADATA1(F)=0.0
             ENDIF

             IF (GRID(C,R)%ETADATA1(F).GT.140.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'EACPCP OUT OF RANGE (',
     &          GRID(C,R)%ETADATA1(F),') at (',c,',',r,
     &          ') CORRECTING TO 140'
              GRID(C,R)%ETADATA1(F)=140.0
             ENDIF

            ENDIF
           ENDIF
        ENDDO
        ENDDO
        ENDDO  
	ENDIF

        DO F=1,LDAS%NF
         DO C=1,LDAS%NC
          DO R=1,LDAS%NR
!UGRD VGRD PRESN APCPN ACPCP CAPE TMP PRES HGT SPFH UGRD VGRD HPBL TMP SPFH FRICV SFEXC CRAIN PEVAP DSWRF DLWRF
!      data nineoctet/33,34,134,202,63,157,11,1,7,51,33,34,
!     &  221,11,51,253,208,140,228,204,205/

           IF ((F.EQ.7).or.(F.eq.14)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF ((GRID(C,R)%etadata2(F).GT.330.0).OR.
     &          (GRID(C,R)%etadata2(F).LT.213.0)) THEN
              newtemp=SQRT(SQRT(ABS(GRID(C,R)%etadata2(21))/5.67E-08))
              write (*,'(A19,F10.5,A6,I3,A1,I3,A16,F10.5)')
     &              'TEMP OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &           ') CORRECTING TO ',newtemp
              GRID(C,R)%etadata2(F)=newtemp
             ENDIF
            ENDIF
           ENDIF

           IF ((F.EQ.10)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).LE.0.0000001) THEN
              vap=(1.0/100.0) * 6.112 *
     &         10**((7.5*(GRID(C,R)%etadata2(7)-273.15))/
     &         (237.7+(GRID(C,R)%etadata2(7)-273.15)))
              rho=0.62197 * ( vap / ((GRID(C,R)%etadata2(8)/100.0)
     &         + vap*(.062197-1) ) )
              write (*,'(A19,F11.8,A6,I3,A1,I3,A16,F11.8)')
     &              'SPFH OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &              ') CORRECTING TO ',rho
              GRID(C,R)%etadata2(F)=rho
             ENDIF
             IF (GRID(C,R)%etadata2(F).GT.0.03) THEN
              write (*,'(A19,F11.8,A6,I3,A1,I3,A20)')
     &          'SPFH OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 0.03'
              GRID(C,R)%etadata2(F)=0.03
             ENDIF
            ENDIF
           ENDIF

           IF ((F.EQ.15)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).LE.0.0000001) THEN
              vap=(1.0/100.0) * 6.112 *
     &         10**((7.5*(GRID(C,R)%etadata2(14)-273.15))/
     &         (237.7+(GRID(C,R)%etadata2(14)-273.15)))
              rho=0.62197 * ( vap / ((GRID(C,R)%etadata2(3)/100.0)
     &         + vap*(.062197-1) ) )
              write (*,'(A19,F11.8,A6,I3,A1,I3,A16,F11.8)')
     &              'SPFH OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &              ') CORRECTING TO ',rho
              GRID(C,R)%etadata2(F)=rho
             ENDIF
             IF (GRID(C,R)%etadata2(F).GT.0.03) THEN
              write (*,'(A19,F11.8,A6,I3,A1,I3,A20)')
     &          'SPFH OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 0.03'
              GRID(C,R)%etadata2(F)=0.03
             ENDIF
            ENDIF
           ENDIF

           IF (F.EQ.20) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).GT.1367) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'EDSWRF OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 1367'
              GRID(C,R)%etadata2(F)=1367
             ENDIF
            ENDIF
           ENDIF

           IF (F.EQ.21) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).GT.750) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'EDLWRF OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 750'
              GRID(C,R)%etadata2(F)=750
             ENDIF
            ENDIF
           ENDIF

           IF ((F.EQ.1).or.(F.eq.11)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).LT.-75) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'UGRD OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &              ') CORRECTING TO -75'
              GRID(C,R)%etadata2(F)=-75
             ENDIF

             IF (GRID(C,R)%etadata2(F).GT.75) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'UGRD OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 75'
              GRID(C,R)%etadata2(F)=75
             ENDIF

            ENDIF
           ENDIF
           IF ((F.EQ.2).or.(F.eq.12)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).LT.-75) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'VGRD OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &              ') CORRECTING TO -75'
              GRID(C,R)%etadata2(F)=-75
             ENDIF

             IF (GRID(C,R)%etadata2(F).GT.75) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'VGRD OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 75'
              GRID(C,R)%etadata2(F)=75
             ENDIF

            ENDIF
           ENDIF

           IF ((F.EQ.8).or.(f.eq.3)) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).LT.5000) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'PRES OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &              ') CORRECTING TO 5000'
              GRID(C,R)%etadata2(F)=5000
             ENDIF

             IF (GRID(C,R)%etadata2(F).GT.110000) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'PRES OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 110000'
              GRID(C,R)%etadata2(F)=110000
             ENDIF

            ENDIF
           ENDIF

           IF (F.EQ.4) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).LT.0.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'EAPCP OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &              ') CORRECTING TO 0'
              GRID(C,R)%etadata2(F)=0.0
             ENDIF

             IF (GRID(C,R)%etadata2(F).GT.140.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'EAPCP OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 140'
              GRID(C,R)%etadata2(F)=140.0
             ENDIF

            ENDIF
           ENDIF

           IF (F.EQ.5) THEN
            IF ((GRID(C,R)%FMASK.GE.1.0).AND.
     &      (GRID(C,R)%etadata2(F).NE.LDAS%UDEF)) THEN
             IF (GRID(C,R)%etadata2(F).LT.0.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A21)')
     &              'EACPCP OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &              ') CORRECTING TO 0'
              GRID(C,R)%etadata2(F)=0.0
             ENDIF

             IF (GRID(C,R)%etadata2(F).GT.140.0) THEN
              write (*,'(A21,F11.8,A6,I3,A1,I3,A20)')
     &          'EACPCP OUT OF RANGE (',
     &          GRID(C,R)%etadata2(F),') at (',c,',',r,
     &          ') CORRECTING TO 140'
              GRID(C,R)%etadata2(F)=140.0
             ENDIF

            ENDIF
           ENDIF
        ENDDO
        ENDDO
        ENDDO

	   ENDIF    ! IF FINDTIME2=1

	BTIME=LDAS%ETATIME1
	CALL TIME2DATE(BTIME,BDOY,GMT1,BYR,BMO,BDA,BHR,BMN)
        BTIME=LDAS%ETATIME2
        CALL TIME2DATE(BTIME,BDOY,GMT2,BYR,BMO,BDA,BHR,BMN)

!===  Interpolate Data in Time
         WT1=(LDAS%ETATIME2-LDAS%TIME)/(LDAS%ETATIME2-LDAS%ETATIME1)
	   WT2=1.0-WT1
c	print *,'wt1=',wt1
c	print *,'wt2=',wt2
	   DO F=1,LDAS%NF


	    IF(F.eq.20) then     !Shortwave
c	print *,'LDAS%GMT,GMT1,GMT2'
c	print *,LDAS%GMT,GMT1,GMT2
	     IF (LDAS%SHORTFLAG.EQ.1) THEN    !Got Instantaneous SW
	      DO C=1,LDAS%NC
		 DO R=1,LDAS%NR
		   ZDOY=LDAS%DOY

                   CALL ZTERP(1,GRID(C,R)%LAT,GRID(C,R)%LON,
     &                GMT1,GMT2,LDAS%GMT,ZDOY,
     &                ZW1,ZW2,CZB,CZE,CZM,LDAS,GRID,czmean)
                   GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)*ZW1+
     &                                  GRID(C,R)%ETADATA2(F)*ZW2
c       In cases of small cos(zenith) angles, use linear weighting
C       to avoid overly large weights
                   
           IF ( (GRID(C,R)%FORCING(F).gt.GRID(C,R)%ETADATA1(F))
     &     .AND. (GRID(C,R)%FORCING(F).gt.GRID(C,R)%ETADATA2(F))
     &     .AND. (CZB.LT.0.1.OR.CZE.LT.0.1))THEN
		   GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)*WT1+
     &		               GRID(C,R)%ETADATA2(F)*WT2
	   ENDIF
                 if (GRID(C,R)%FORCING(F).gt.1367) then
                  print *,'warning, Instantaneous SW RADIATION TOO HIGH'
	          print *,'using linear weighting'
                  print *,'it is',GRID(C,R)%FORCING(F)
                  print *,'ETADATA1=',GRID(C,R)%ETADATA1(F)
                  print *,'ETADATA2=',GRID(C,R)%ETADATA2(F)
                  print *,'ZW1=',ZW1,'ZW2=',ZW2
                  print *,'WT1=',WT1,'WT2=',WT2
                   GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)*WT1+
     &                         GRID(C,R)%ETADATA2(F)*WT2
                 endif

                IF(GRID(C,R)%FIMASK.EQ.0) GRID(C,R)%FORCING(F)=LDAS%UDEF

                 if ((GRID(C,R)%FORCING(F).ne.LDAS%UDEF).and.
     &               (GRID(C,R)%FORCING(F).lt.0) ) then
                    print *,'1 warning!!!!  SW radiation is negative!!'
                    print *,'sw=',GRID(C,R)%FORCING(F),'...negative'
                    print *,'eta1 equaled',GRID(C,R)%ETADATA1(F)
                    print *,'eta2 equaled',GRID(C,R)%ETADATA2(F)
		    print *,'forcing mask=',GRID(C,R)%FIMASK
		    print *,'i,j,f=',c,r,f
                    stop
                 endif
                   GRID(C,R)%ETASW=GRID(C,R)%FORCING(F)
		 ENDDO
		ENDDO
	    ENDIF
            
             IF (LDAS%SHORTFLAG.EQ.2) THEN    !Got Time Averaged SW
	gmt2=gmt2-0.05
	if (gmt2.lt.0) gmt2=24+gmt2
	outradmask=0.0
              DO C=1,LDAS%NC
                 DO R=1,LDAS%NR
                   ZDOY=LDAS%DOY
        if ((mod(ldas%hr,3).ne.0).or.(ldas%mn.ne.0)) then
!       not on 3 hour file dividing line, so use first file in memory
                   CALL ZTERP(0,GRID(C,R)%LAT,GRID(C,R)%LON,
     &                GMT1,GMT2,LDAS%GMT,ZDOY,
     &                ZW1,ZW2,CZB,CZE,CZM,LDAS,GRID,czmean)

         GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)*ZW1
!   QC FOR SMALL ZENITH ANGLES
         IF(GRID(C,R)%FORCING(F) .GT. 1367.0*CZM) THEN
         CRATIO = GRID(C,R)%ETADATA1(F)/(1367.0*CZMEAN)
         GRID(C,R)%FORCING(F) = CRATIO * (1367.0 * CZM)
       ENDIF

	else
	if((c.eq.1).and.(r.eq.1))  then
         gmt2=gmt2+0.05
         if (gmt2.ge.24) gmt2=gmt2-24
         tempgmt1=gmt2
         tempgmt2=hr2
         tempgmt2=tempgmt2-0.05
         if (tempgmt2.lt.0) tempgmt2=24+tempgmt2
	endif

!       on 3 hour file dividing line, so use second file in memory

                   CALL ZTERP(0,GRID(C,R)%LAT,GRID(C,R)%LON,
     &                tempGMT1,tempgmt2,LDAS%GMT,ZDOY,
     &                ZW1,ZW2,CZB,CZE,CZM,LDAS,GRID,czmean)
	 GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA2(F)*ZW1
!QC FOR SMALL ZENITH ANGLES

         IF(GRID(C,R)%FORCING(F) .GT. 1367.0*CZM) THEN
	print *,'input was gmt1,gmt2,model gmt',
     &  tempGMT1,tempgmt2,LDAS%GMT
	print *,'value too high, its',GRID(C,R)%FORCING(F), 
     &  'max is',1367.0*CZM, 'etadata2=',GRID(C,R)%ETADATA2(F)
         CRATIO = GRID(C,R)%ETADATA2(F)/(1367.0*CZMEAN)
!         GRID(C,R)%FORCING(F) = CRATIO*(1367.0*CZM)
!	outradmask(c,r)=(1367.0*CZM)-GRID(C,R)%FORCING(F)
!	check=1
         GRID(C,R)%FORCING(F) = 1367.0*CZM
	ENDIF
	endif

                IF(GRID(C,R)%FIMASK.EQ.0) GRID(C,R)%FORCING(F)=LDAS%UDEF

	         if ((GRID(C,R)%FORCING(F).ne.LDAS%UDEF).and.
     &               (GRID(C,R)%FORCING(F).lt.0) ) then
                    print *,'2 warning!!!!  SW radiation is negative!!'
		    print *,'sw=',GRID(C,R)%FORCING(F),'...negative'
                    print *,'eta2 equaled',GRID(C,R)%ETADATA2(F)
                    print *,'forcing mask=',GRID(C,R)%FIMASK

                    stop
	         endif

	         if (GRID(C,R)%FORCING(F).gt.1367) then
                  print *,'warning, SW AVG RADIATION TOO HIGH!!'
                  print *,'it is',GRID(C,R)%FORCING(F),'at',c,r,
     &  'czm,zw1',czm,zw1
                  print *,'etadata1,etadata2=',
     &              GRID(C,R)%ETADATA1(F),GRID(C,R)%ETADATA2(F)
        if ((mod(ldas%hr,3).ne.0).or.(ldas%mn.ne.0)) then
!       not on 3 hour file dividing line, so use first file in memory
                  print *,'setting equal to etadata1',
     &            GRID(C,R)%ETADATA1(F)
                  GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)
	else
!       on 3 hour file dividing line, so use second file in memory
                  print *,'setting equal to etadata2',
     &            GRID(C,R)%ETADATA2(F)
                  GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA2(F)
		 endif
	         endif
                   GRID(C,R)%ETASW=GRID(C,R)%FORCING(F)
                 ENDDO
                ENDDO

            ENDIF

!	if (check.eq.1) then
!	open (unit=45,file='radchange.bin',form='unformatted',
!     &  status='unknown')
!        write(45) ((GRID(C,R)%FORCING(20),c=1,464),r=1,224)
!	write(45) outradmask
!	close(45)
!	stop	
!	endif
	    ELSE IF(F.eq.4.or.F.eq.5.or.F.eq.18.or.F.eq.19) then    ! A precip variable Do Block Interpolation
	      DO C=1,LDAS%NC
		 DO R=1,LDAS%NR
		  GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)
		 ENDDO
		ENDDO

	    ELSE IF (F.eq.21) then     !Longwave
		 IF (LDAS%LONGFLAG.EQ.1) THEN    !Got Instantaneous LW
		  DO C=1,LDAS%NC
		   DO R=1,LDAS%NR
		     GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)*WT1+
     &		                          GRID(C,R)%ETADATA2(F)*WT2 
                ENDDO
		   ENDDO
		 ENDIF    

             IF (LDAS%LONGFLAG.EQ.2) THEN    !Got Time Averaged LW
              DO C=1,LDAS%NC
               DO R=1,LDAS%NR
	if ((mod(ldas%hr,3).ne.0).or.(ldas%mn.ne.0)) then
!	not on 3 hour file dividing line, so use first file in memory	
c	if((c.eq.1).and.(r.eq.1)) 
c     &  print *,'not on 3 hour file dividing line, so use first LW'
		GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)
	else
! 	on 3 hour file dividing line, so use second file in memory
!        if((c.eq.1).and.(r.eq.1)) 
!     &	print *,'on 3 hour file dividing line, so use second LW'
                GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA2(F)
	endif
               ENDDO
              ENDDO
	     ENDIF
	      
	    ELSE	   !Linearly interpolate everything else	
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
               GRID(C,R)%FORCING(F)=GRID(C,R)%ETADATA1(F)*WT1+
     &                           GRID(C,R)%ETADATA2(F)*WT2
             ENDDO
            ENDDO  	    
	     
	    

	    ENDIF
	   ENDDO   !the F loop
	  
!=== ADJUST PRECIP TO VALUE PER SEC DATA SHOULD COME IN AS VALUE PER
!=== TIME PERIOD THREE HOUR FOR ETA	 
!--- POSSIBLE ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR  ERROR ERROR ERROR ERROR  ERROR ERROR ERROR ERROR 


         if((ldas%ftype.eq.1).or.(ldas%ftype.eq.2)) then !We got ETA 3hr or EDAS data
c	print *,'ldas%ftype is  or 2'
         DO C=1,LDAS%NC
	    DO R=1,LDAS%NR
	     GRID(C,R)%FORCING(4)=GRID(C,R)%FORCING(4)/(3.0*60.0*60.0)
	     GRID(C,R)%FORCING(5)=GRID(C,R)%FORCING(5)/(3.0*60.0*60.0) 
	    ENDDO
	   ENDDO
	 endif

         if (ldas%ftype.eq.3) then !We got ETA 6hr or EDAS data
c	print *,'ldas%ftype is three, dividing by 6 for six hour precip'
         DO C=1,LDAS%NC
            DO R=1,LDAS%NR
             GRID(C,R)%FORCING(4)=GRID(C,R)%FORCING(4)/(6.0*60.0*60.0)
             GRID(C,R)%FORCING(5)=GRID(C,R)%FORCING(5)/(6.0*60.0*60.0)
            ENDDO
           ENDDO
         endif
 

C
84       format('NOW',i4,4i3,2x,'PVT ',a22,' NXT ',a22)
          if(FINDTIME2.eq.1) then
	      write(83,*) 'NXT-Name: ',NAME
            write(83,84) YR1,MO1,DA1,zzHR1,zzMN1,
     &     LDAS%EVT1,LDAS%EVT2
           endif
c
	 
c	print *,'name=',name 
	   RETURN 
	   END
	   

C
C
!!!!!SSSSS  SUBROUTINES    SUBROUTINES    SUBROUTINES   SSSSS
C
C
!=======================================================
!
!  DESCRIPTION:
!   This subroutine puts together EDAS file name
!
!=======================================================
	   
	   SUBROUTINE ETAEDASFILE(NAME,EDASDIR,YR,MO,DA,HR,LFB)
	   
	   IMPLICIT  NONE
	   
!=== Local Variables
          CHARACTER*80 NAME
	    CHARACTER*80 EDASDIR
	    INTEGER YR,MO,DA,HR,I,C
	    
	    INTEGER LFB   ! 1 = look backward , 2 = look forward in time
	    
	    INTEGER UYR,UMO,UDA,UHR,UMN,USS,TA
	    INTEGER SUYR  !two digit year
	    REAL*8 TIME1,DUMTIME
	    REAL GMT1,DUMGMT
	    INTEGER DOY1,DUMDOY
	    CHARACTER*1 FNAME(80),FBASE(80),FSUBS(80),UHRCHAR(2)
	    CHARACTER*1 FTIME(8),FDIR(15)
	    CHARACTER*2 CHINSERT
	    CHARACTER*2  SPCODE(4)   !Special Code
	    CHARACTER*2  HRCODE
	    INTEGER PNTSPC  !points to correct Special Code
	    INTEGER IBREAK
	    
	    DATA SPCODE/'12','09','06','03'/
	    
!==== End Variable Definition
     
    
     
!=== Put together filename
64       FORMAT(I4)
65       FORMAT(4a1)
91       FORMAT(A4,i3,A11,I3)
92       FORMAT(80A1)
93       FORMAT(A80)
94       FORMAT(I4,I2,I2)
95       FORMAT(8a1)
96       FORMAT(a80)
97       FORMAT(8A1,2A1,A9)
98       FORMAT(a1,i4,a1,i4,i2,i2,a1)
99       FORMAT(15a1)
C
!==== Make variables for the time used to create the file
!====  We don't want these variables being passed out
          UYR=YR
	    UMO=MO
	    UDA=DA

          UMN=0
	    USS=0	 
C  The hour needs to be change to
C  a multiple of 3 hour,
C  	    
	   UHR=(HR/3)*3

C        print *,'before uyr,umo,uda,uhr,umn,uss',uyr,umo,uda,uhr,umn,uss
         CALL TICK(DUMTIME,DUMDOY,DUMGMT,UYR,UMO,UDA,UHR,UMN,USS,
     &    -0)
c	print *,'after uyr,umo,uda,uhr,umn,uss',uyr,umo,uda,uhr,umn,uss 
C 
C
C  If the time is 12 or later the file is time stamped
C  with the next day.  So check for that first
 
c narr-a_221_20030502_1500_000.grb
c         if(UHR.ge.12) then
c	     TA=60*60*12         
c	     CALL TICK(TIME1,DOY1,GMT1,UYR,UMO,UDA,UHR,UMN,USS,TA)
c	     HRCODE='00'
c	   else
c	     HRCODE='12'
c         endif
	     PNTSPC=(UHR/3)+1
	   
         OPEN(90,FILE='temp',FORM='FORMATTED',
     &	   ACCESS='DIRECT',RECL=80)

         WRITE(90,96,REC=1) EDASDIR
         READ(90,92,REC=1) (FBASE(I),I=1,80)

         if(UYR.ge.2000) SUYR=UYR-2000
	   if(UYR.le.1999) SUYR=UYR-1900
 
         WRITE(90,98,rec=1)'/',UYR,'/',UYR,UMO,UDA,'/'
         READ(90,99,REC=1) FDIR
         DO I=1,15
          IF(FDIR(I).EQ.(' ')) FDIR(I)='0'
         ENDDO

         WRITE(90,'(I2)',rec=1)UHR
         READ(90,'(2A1)',REC=1) UHRCHAR
         DO I=1,2
          IF(UHRCHAR(I).EQ.(' ')) UHRCHAR(I)='0'
         ENDDO

	  
C         IBREAK=10     !all files have 4 digit years now
c	  	print *,'uyr,umo,uda',uyr,umo,uda 
	     WRITE(90,94,REC=1) UYR,UMO,UDA
           READ(90,95,REC=1) FTIME
           DO I=1,8
             IF(FTIME(I).EQ.(' ')) FTIME(I)='0'	     
           ENDDO	

c	print *,'ftime=',ftime   
c      	print *,'uhr=',uhr 
c narr-a_221_20030502_1500_000.grb 
c 97       FORMAT(A11,8A1,A1,2A1,A10)

c1999113000.NARR.grb
         WRITE(90,97,REC=1)FTIME,UHRCHAR,'.NARR.grb'
         READ(90,92,REC=1) (FSUBS(I),I=1,19)
        C=0
        DO I=1,80
          IF(FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1 
        ENDDO
	
	   WRITE(90,92,REC=1) (FBASE(I),I=1,C),(FDIR(I),I=1,15),
     &  (FSUBS(I),I=1,19)

       READ(90,93,REC=1) NAME
c	print *,'name=',name
	 CLOSE(90)
	 RETURN
	 END     
         	    	   
	  
!=======================================================
!
!  DESCRIPTION:
!   This subroutine puts together EDAS file name
!
!=======================================================

           SUBROUTINE ETAEDASPRECIPFILE(NAME,EDASDIR,YR,MO,DA,HR,LFB)

           IMPLICIT  NONE

!=== Local Variables
          CHARACTER*80 NAME
            CHARACTER*80 EDASDIR
            INTEGER YR,MO,DA,HR,I,C

            INTEGER LFB   ! 1 = look backward , 2 = look forward in time

            INTEGER UYR,UMO,UDA,UHR,UMN,USS,TA
            INTEGER SUYR  !two digit year
            REAL*8 TIME1,DUMTIME
            REAL GMT1,DUMGMT
            INTEGER DOY1,DUMDOY
            CHARACTER*1 FNAME(80),FBASE(80),FSUBS(80),UHRCHAR(2)
            CHARACTER*1 FTIME(8),FDIR(15)
            CHARACTER*2 CHINSERT
            CHARACTER*2  SPCODE(4)   !Special Code
            CHARACTER*2  HRCODE
            INTEGER PNTSPC  !points to correct Special Code
            INTEGER IBREAK

            DATA SPCODE/'12','09','06','03'/

!==== End Variable Definition



!=== Put together filename
64       FORMAT(I4)
65       FORMAT(4a1)
91       FORMAT(A4,i3,A11,I3)
92       FORMAT(80A1)
93       FORMAT(A80)
94       FORMAT(I4,I2,I2)
95       FORMAT(8a1)
96       FORMAT(a80)
97       FORMAT(8A1,2A1,A9)
98       FORMAT(a1,i4,a1,i4,i2,i2,a1)
99       FORMAT(15a1)
C
!==== Make variables for the time used to create the file
!====  We don't want these variables being passed out
          UYR=YR
            UMO=MO
            UDA=DA

          UMN=0
            USS=0
C  The hour needs to be change to
C  a multiple of 3 hour,
C
           UHR=(HR/3)*3

C        print *,'before uyr,umo,uda,uhr,umn,uss',uyr,umo,uda,uhr,umn,uss
         CALL TICK(DUMTIME,DUMDOY,DUMGMT,UYR,UMO,UDA,UHR,UMN,USS,
     &    -10800)
c       print *,'after uyr,umo,uda,uhr,umn,uss',uyr,umo,uda,uhr,umn,uss
C
C
C  If the time is 12 or later the file is time stamped
C  with the next day.  So check for that first

c narr-a_221_20030502_1500_000.grb
c         if(UHR.ge.12) then
c            TA=60*60*12
c            CALL TICK(TIME1,DOY1,GMT1,UYR,UMO,UDA,UHR,UMN,USS,TA)
c            HRCODE='00'
c          else
c            HRCODE='12'
c         endif
             PNTSPC=(UHR/3)+1

         OPEN(90,FILE='temp',FORM='FORMATTED',
     &     ACCESS='DIRECT',RECL=80)

         WRITE(90,96,REC=1) EDASDIR
         READ(90,92,REC=1) (FBASE(I),I=1,80)

         if(UYR.ge.2000) SUYR=UYR-2000
           if(UYR.le.1999) SUYR=UYR-1900

         WRITE(90,98,rec=1)'/',UYR,'/',UYR,UMO,UDA,'/'
         READ(90,99,REC=1) FDIR
         DO I=1,15
          IF(FDIR(I).EQ.(' ')) FDIR(I)='0'
         ENDDO

         WRITE(90,'(I2)',rec=1)UHR
         READ(90,'(2A1)',REC=1) UHRCHAR
         DO I=1,2
          IF(UHRCHAR(I).EQ.(' ')) UHRCHAR(I)='0'
         ENDDO


C         IBREAK=10     !all files have 4 digit years now
c               print *,'uyr,umo,uda',uyr,umo,uda
             WRITE(90,94,REC=1) UYR,UMO,UDA
           READ(90,95,REC=1) FTIME
           DO I=1,8
             IF(FTIME(I).EQ.(' ')) FTIME(I)='0'
           ENDDO

c       print *,'ftime=',ftime
c       print *,'uhr=',uhr
c narr-a_221_20030502_1500_000.grb
c 97       FORMAT(A11,8A1,A1,2A1,A10)

c1999113000.NARR.grb
         WRITE(90,97,REC=1)FTIME,UHRCHAR,'.NARR.grb'
         READ(90,92,REC=1) (FSUBS(I),I=1,19)
        C=0
        DO I=1,80
          IF(FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1
        ENDDO

           WRITE(90,92,REC=1) (FBASE(I),I=1,C),(FDIR(I),I=1,15),
     &  (FSUBS(I),I=1,19)

       READ(90,93,REC=1) NAME
c       print *,'name=',name
         CLOSE(90)
         RETURN
         END
 
!=======================================================
!
!  DESCRIPTION:
!   This subroutine puts toghether ETA 3hr file name
!
!=======================================================
	
	   
	   SUBROUTINE ETA3HRFILE(NAME,ETA3HRDIR,YR,MO,DA,HR,
     &	   LFB,FCA,Prevname,Precflag)

	   IMPLICIT  NONE
	   
!=== Local Variables
          CHARACTER*80 NAME
	    CHARACTER*80 Prevname
	    CHARACTER*40 ETA3HRDIR
	    INTEGER YR,MO,DA,HR,I,C
	    
	    INTEGER LFB       ! 1 = look back in time, 2 = look forward in time
	    INTEGER FCA       ! 1 = roll back 12 hrs, 2 = 24hrs, 3= 36hrs, 4= 48hrs
	    INTEGER Precflag  ! 0=no,1=yes  subtract out previous precip amount
	    INTEGER UYR,UMO,UDA,UHR,UMN,USS,TA
	    INTEGER UDOYR   !local day of year
	    REAL UGMT
	    INTEGER SUYR  !two digit year
	    REAL*8 TIME1
	    REAL GMT1
	    INTEGER DOY1
	    CHARACTER*1 FNAME(80),FBASE(80),FSUBS(80)
	    CHARACTER*1 FTIME(10),FDIR(8)
	    CHARACTER*2 CHINSERT
            CHARACTER*2  SPCODE(28)   !Special Code
            INTEGER PNTSPC  !points to correct Special Code
            DATA SPCODE/'03','06','09','12','15','18','21',
     &     '24','27','30','33','36','39','42','45','48',
     &     '51','54','57','60','63','66','69','72','75',
     &     '78','81','84' /
C	    write(*,*) 'In Make 3 hr name'
!==== End Variable Definition
     
!=== Put together filename
91       FORMAT(A4,i3,A11,I3)
92       FORMAT(80A1)
93       FORMAT(A80)
94       FORMAT(I4,I2,I2,a2)
95       FORMAT(10a1)
96       FORMAT(a40)
97       FORMAT(a1,a4,a1,a7,a2,a8)
98       FORMAT(a1,i4,i2,a1)
99       FORMAT(8a1)
C
!==== Make variables for the time used to create the file
!====  We don't want these variables being passed out
          UYR=YR
	    UMO=MO
	    UDA=DA
	    UMN=0
	    USS=0     
C  The hour needs to be change to
C  a multiple of 3 hour
C  	    
	     UHR=(HR/3)*3    
C
C  If the hour is 00 then we need the previous day of the year
C  If the hour is 03,06,09,12 then today using '00'
C  If the hour is 15,18,21 then today using '00'
 
         if(UHR.eq.0) then
C	      print*,'inside'
	      TA=(-30)*(60)
C		print*,'TIME1',TIME1
C		print*,'UYR',UYR,' UMO',UMO,' UDA', UDA
C		print*,'UMN',UMN,'  USS',USS,' TA',TA
	      CALL TICK(TIME1,UDOYR,UGMT,UYR,UMO,UDA,UHR,UMN,USS,TA)
		CHINSERT='12'
		PNTSPC=4+(FCA*4)
         else if(UHR.ge.3.and.UHR.le.12) then
	      CHINSERT='00'
		PNTSPC=(UHR/3)+(FCA*4)
         else              !UHR equals 15,18,21                   
	      CHINSERT='12'
		PNTSPC=((UHR-12)/3)+(FCA*4)
         endif
C
C
 
         OPEN(90,FILE='temp',FORM='FORMATTED',
     &	   ACCESS='DIRECT',RECL=80)

         WRITE(90,96,REC=1) ETA3HRDIR
         READ(90,92,REC=1) (FBASE(I),I=1,80)
 	
    
         WRITE(90,98,rec=1)'/',UYR,UMO,'/'
         READ(90,99,REC=1) FDIR
         DO I=1,8
          IF(FDIR(I).EQ.(' '))FDIR(I)='0'
         ENDDO
         if(UYR.ge.2000) SUYR=UYR-2000
	   if(UYR.le.1999) SUYR=UYR-1900
         WRITE(90,94,REC=1) UYR,UMO,UDA,CHINSERT
         READ(90,95,REC=1)FTIME
         DO I=1,10
           IF(FTIME(I).EQ.(' '))FTIME(I)='0'
         ENDDO
	   
         WRITE(90,97,REC=1) '.AWIPSF',SPCODE(PNTSPC),'.tm00.sg'
         READ(90,92,REC=1) (FSUBS(I),I=1,17)
    
    
         C=0
         DO I=1,80
           IF(FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1 
         ENDDO
         WRITE(90,92,REC=1) (FBASE(I),I=1,C),(FDIR(I),I=1,8),
     &  (FTIME(I),I=1,10),(FSUBS(I),I=1,17)
         READ(90,93,REC=1) NAME
	   CLOSE(90)

C   If the hour is NOT 3 or 15 you'll need the previous
C   forecast period to subtract out the precip
         if(UHR.ne.03.and.UHR.ne.15) then
           precflag=1
           PNTSPC=PNTSPC-1
C	     
           OPEN(90,FILE='temp',FORM='FORMATTED',
     &	   ACCESS='DIRECT',RECL=80)   
           WRITE(90,97,REC=1) '.AWIPSF',SPCODE(PNTSPC),'.tm00.sg'
           READ(90,92,REC=1) (FSUBS(I),I=1,17)
           C=0
           DO I=1,80
             IF(FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1 
           ENDDO
           WRITE(90,92,REC=1) (FBASE(I),I=1,C),(FDIR(I),I=1,8),
     &    (FTIME(I),I=1,10),(FSUBS(I),I=1,17)
           READ(90,93,REC=1) Prevname
	     CLOSE(90)
         else
	     precflag=0
	   endif
C         print*,'done with 3 hr construction'
	   
	   RETURN
	   END
!=======================================================
!
!  DESCRIPTION:
!   This subroutine puts toghether ETA 6hr file name
!
!=======================================================
	   
	   SUBROUTINE ETA6HRFILE(NAME,ETA6HRDIR,YR,MO,DA,HR,
     &	   LFB,FCA,Prevname,Precflag)

	   IMPLICIT  NONE
	   
!=== Local Variables
          CHARACTER*80 NAME
	    CHARACTER*80 Prevname
	    CHARACTER*40 ETA6HRDIR
	    INTEGER YR,MO,DA,HR,I,C
	    
	    INTEGER LFB       ! 1 = look back in time, 2 = look forward in time
	    INTEGER FCA       ! Forecast amount ( 1= 12 hours back, 2 = 24 hrs back, etc)
	    INTEGER Precflag  ! (0=no 1=yes) subtract out previous precip amount
	    INTEGER UYR,UMO,UDA,UHR,UMN,USS,TA
	    INTEGER UDOYR  !local day of year
	    REAL UGMT
	    INTEGER SUYR   !two digit year
	    INTEGER IBREAK   
	    REAL*8 TIME1
	    REAL GMT1
	    INTEGER DOY1
	    CHARACTER*1 FNAME(80),FBASE(80),FSUBS(80)
	    CHARACTER*1 FTIME(10),FDIR(8)
	    CHARACTER*2  CHINSERT
            CHARACTER*2  FCACODE(14)   !Special Code
            INTEGER PNTSPC  !points to correct Special Code

            DATA FCACODE/'06','12','18','24','30','36','42','48',
     &                   '54','60','66','72','78','84'/
 
!==== End Variable Definition
     
!=== Put together filename
91       FORMAT(A4,i3,A11,I3)
92       FORMAT(80A1)
93       FORMAT(A80)
94      FORMAT(I4,I2,I2,a2)
95      FORMAT(10a1)
96       FORMAT(a40)
97       FORMAT(a7,a2,a8)
98       FORMAT(a1,i4,i2,a1)
99       FORMAT(8a1)
C
!==== Make variables for the time used to create the file
!====  We don't want these variables being passed out
          UYR=YR
	    UMO=MO
	    UDA=DA
	    UMN=0
	    USS=0     
C  The hour needs to be change to
C  a multiple of 6 hour
C  	    
	    UHR=(HR/6)*6
	   
C
C    
C  If the hour is 00 then we need the previous day of the year
C  If the hour is 03,06,09,12 then today using '00'
C  If the hour is 15,18,21 then today using '00'
 
         if(UHR.eq.0) then
	      TA=(-30)*(60)
	      CALL TICK(TIME1,UDOYR,UGMT,UYR,UMO,UDA,UHR,UMN,USS,TA)
		CHINSERT='12'
		PNTSPC=2+(FCA*2)
         else if(UHR.eq.6) then 
	      CHINSERT='00'
		PNTSPC=1+(FCA*2)
	   else if(UHR.eq.12) then
	      CHINSERT='00'
		PNTSPC=2+(FCA*2)
         else  !equals 18
	      CHINSERT='12'
		PNTSPC=1+(FCA*2)
         endif
C
C
 
         OPEN(90,FILE='temp',FORM='FORMATTED',
     &	   ACCESS='DIRECT',RECL=80)

         WRITE(90,96,REC=1) ETA6HRDIR
         READ(90,92,REC=1) (FBASE(I),I=1,80)
 	
         if(UYR.ge.2000) SUYR=UYR-2000
	   if(UYR.le.1999) SUYR=UYR-1900
	   
	   
         WRITE(90,98,rec=1)'/',UYR,UMO,'/'
         READ(90,99,REC=1) FDIR
         DO I=1,8
          IF(FDIR(I).EQ.(' '))FDIR(I)='0'
         ENDDO

         IBREAK=10
	   

	     WRITE(90,94,REC=1) UYR,UMO,UDA,CHINSERT
	     READ(90,95,REC=1) FTIME
	     DO I=1,IBREAK
	       IF(FTIME(I).EQ.(' ')) FTIME(I)='0'
	     ENDDO

	   	   
         WRITE(90,97,REC=1) '.AWIPSF',FCACODE(PNTSPC),'.tm00.sg'
         READ(90,92,REC=1) (FSUBS(I),I=1,17)
    
    
         C=0
         DO I=1,80
           IF(FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1 
         ENDDO
	   
	    WRITE(90,92,REC=1) (FBASE(I),I=1,C),(FDIR(I),I=1,8),
     &  (FTIME(I),I=1,10),(FSUBS(I),I=1,17)
    
         READ(90,93,REC=1) NAME
	   CLOSE(90)

C   Now we need to determine if a name should be made
C   for the previous time period so the precip is for a 6 hour period
         if(UHR.eq.12.or.UHR.eq.00) then
C
	     PNTSPC=PNTSPC-1
	     Precflag=1
C	     
           OPEN(90,FILE='temp',FORM='FORMATTED',
     &	   ACCESS='DIRECT',RECL=80)	 
           WRITE(90,97,REC=1) '.AWIPSF',FCACODE(PNTSPC),'.tm00.sg'
           READ(90,92,REC=1) (FSUBS(I),I=1,17)    
           C=0
           DO I=1,80
             IF(FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1 
           ENDDO
           WRITE(90,92,REC=1) (FBASE(I),I=1,C),(FDIR(I),I=1,8),
     &    (FTIME(I),I=1,10),(FSUBS(I),I=1,17)
           READ(90,93,REC=1) Prevname
	     CLOSE(90)
	   else
	     Precflag=0  
	   endif	   
	            
C         print*,'Done with 6 hr construction'



           RETURN
           END

