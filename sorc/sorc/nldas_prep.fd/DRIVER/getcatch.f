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
! getcatch.f: 
!
! DESCRIPTION:
!  Opens, reads, and interpolates CATCHMENT-LDAS forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data 
!
!
! REVISION HISTORY:
!  20 Feb 2001: Brian Cosgrove; Initial code based heavily on geteta.f
!                               Added CZM into call for ZTERP subroutine
!=========================================================================
      SUBROUTINE GETcatch(LDAS,GRID)
	
	USE ldas_module     ! LDAS non-model-specific 1-D variables
	USE grid_module     ! LDAS non-model-specific grid variables
	IMPLICIT NONE
	type (ldasdec) LDAS
	type (griddec) GRID(LDAS%NC,LDAS%NR)
	
	
!==== Local Variables=======================
       INTEGER  ii,jj
       INTEGER C,R,F,FERROR,TRY,ZDOY,cc,rr
	 INTEGER YR1,MO1,DA1,HR1,MN1,SS1,DOY1,TS1,BDOY,BYR,BMO
	 INTEGER YR2,MO2,DA2,HR2,MN2,SS2,DOY2,TS2,BDA,BHR,BMN
	 INTEGER YR26,MO26,DA26,HR26,MN26,SS26,DOY26,TS26
	 INTEGER zzHR1,zzMN1  !for diagnostics
	 INTEGER TS1H   !to go back twelve hours
	 INTEGER TS2D   !to go forward two days (restoring the role back)
	 INTEGER TS1Day   !to go backwards 1 Day
	 REAL*8 TIME1,TIME2,FAKETIME,TIME26
	 REAL*8 DUMBTIME1,DUMBTIME2,DUMBTIME26
	 REAL*8 TIMENOW,BTIME
	 REAL*8 FILETIME1,FILETIME2    !these are the time that the file correspond too
	 CHARACTER*80 NAME(9)
	 CHARACTER*80 Prevname  !previous time period's file name
	 REAL WT1,WT2,ZW1,ZW2,CZB,CZE,GMT1,GMT2,GMT26,CZM
         REAL vp(ldas%nc,ldas%nr)
	 INTEGER LFB      !this is the loop for looking backward or forward in time
	 INTEGER FCA      !Forecast Amount, 1=12 hrs, 2=24hrs, etc
	 INTEGER ftype      !file type 1=catch, 2=catch3hr, 3=catch6hr
	 INTEGER ORDER    !1=time before, 2=time after
	 
	 INTEGER FINDTIME1     ! 0=don't get new file for 1st time (or 2nd)
	 INTEGER FINDTIME2     ! 1=Get a new file 1st time (or 2nd)
	 INTEGER MOVETIME      ! 1=Move Time 2 data into Time 1
	 INTEGER SUAVE
	 INTEGER COUNT
	 CHARACTER*40 var
c	print *,'in getcatch'	 
	 SUAVE=0
	 
	 
	 
C   Assumption will be not to find or move any data
       FINDTIME1=0
	 FINDTIME2=0
	 MOVETIME=0	 
	 
!=== End Variable Definition =======================	 
	
!=== Determine Required catch Data Times (The previous hour and the future hour)
!=== The Adjustment of the Hour and the direction will be done
!=== in the subroutines that generate the names b/c it's different
!=== for the three or six hour time steps 
      YR1=LDAS%YR    !Current Hour
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
      HR1=6*((LDAS%HR)/6)
      MN1=0
      SS1=0
      TS1=0
      CALL TICK(TIME1,DOY1,GMT1,YR1,MO1,DA1,HR1,MN1,SS1,TS1)

      YR2=LDAS%YR    !Next Hour
      MO2=LDAS%MO
      DA2=LDAS%DA
      HR2=6*((LDAS%HR)/6)
	MN2=0
	SS2=0
	TS2=6*60*60
        CALL TICK(TIME2,DOY2,GMT2,YR2,MO2,DA2,HR2,MN2,SS2,TS2)

     
      DUMBTIME1=TIME1
	DUMBTIME2=TIME2

        if(TIMENOW.GT.LDAS%catchTIME2) then
	   MOVETIME=1
	   FINDTIME2=1
	endif

	if(LDAS%TSCOUNT.EQ.1) then     !beginning of the run
c	print *,'at beginning of run#############'
	   FINDTIME1=1
	   FINDTIME2=1
           MOVETIME=0
	endif
      
20    format('MOVE',i2,2x,'FIND1',i2,2x,'FIND2',i2)
    

!==== If FINDTIME2=1 then you want initialize the FSOURCE info
         if(FINDTIME2.eq.1) then
	     LDAS%FSOURCE(16)=0     !catch Status
	   endif

	 
!===  Check to see if required data in memory, if not, read it
	  IF(FINDTIME1.eq.1) THEN !We need to get new TIME1 data
c	print *,'need to get time 1 data#############'
	      ORDER=1        !says that it's the "before" time
            FERROR=0
		TRY=0
		LFB=1
		TS1=-60*60*24
		TS1H=-60*60*12
		  TRY=TRY+1
c	print *,'calling catchfile for',YR1,MO1,DA1,HR1
                VAR='/TEMP-AIR'
                CALL catchFILE(NAME(1),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)
                VAR='/TEMP-DEW'
                CALL catchFILE(NAME(2),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)
                VAR='/RAD-SW'
                CALL catchFILE(NAME(3),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)
                VAR='/RAD-LW'
                CALL catchFILE(NAME(4),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)
                VAR='/WIND'
                CALL catchFILE(NAME(5),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)
                VAR='/WIND'
                CALL catchFILE(NAME(6),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)
                VAR='/PRES-SRF'
                CALL catchFILE(NAME(7),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)
                VAR='/PREC-TOTL'
                CALL catchFILE(NAME(8),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)
                VAR='/PREC-CONV'
                CALL catchFILE(NAME(9),LDAS%catchDIR,YR1,MO1,DA1,HR1,
     &          LFB,count,var)


c        print *,'count=',count

c	print*,'bc called catchfile'
		  CALL RETcatch(ORDER,LDAS,GRID,NAME,FERROR,ftype,
     &           prevname,count)

                    LDAS%CATCHTIME1=TIME1
	   ENDIF    !time1.ne.LDAScatchTIME1
C
!  REPEAT FOR TIME 2

	   IF(MOVETIME.eq.1) THEN !Transfer TIME2 data to TIME1
c	    print*,'Transfer time 2 to time 1 data'
	     LDAS%catchTIME1=LDAS%catchTIME2	
	     FINDTIME2=1 
	     do F=1,9
	      do C=1,LDAS%NC
		 do R=1,LDAS%NR
		   GRID(C,R)%catchDATA1(F)=GRID(C,R)%catchDATA2(F)
		 enddo
		enddo
	     enddo
C	      print*,'done with loop'
C
         ENDIF  ! if MOVETIME=1



C          Beginning point
	  IF(FINDTIME2.eq.1) THEN !We need to get new TIME2 data
c        print *,'need to get time 2 data#############'
c	     print*,'beginning Finding Time 2'
           ORDER=2
		TRY=0
C		print*,'TRY',TRY,'  LFB',LFB
		LFB=1
		TS1=-60*60*24
		TS1H=-60*60*12
C            print*,'TS1',TS1,'TS1H',TS1H,'  LFB',LFB
C		print*,'FERROR',FERROR
		
c		  print*,'TRY is ',TRY
		  TRY=TRY+1
		  
c        print *,'calling catchfile for',YR2,MO2,DA2,HR2
                VAR='/TEMP-AIR'
                CALL catchFILE(NAME(1),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
                VAR='/TEMP-DEW'
                CALL catchFILE(NAME(2),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
                VAR='/RAD-SW'
                CALL catchFILE(NAME(3),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
                VAR='/RAD-LW'
                CALL catchFILE(NAME(4),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
                VAR='/WIND'
                CALL catchFILE(NAME(5),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
                VAR='/WIND'
                CALL catchFILE(NAME(6),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
                VAR='/PRES-SRF'
                CALL catchFILE(NAME(7),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
                VAR='/PREC-TOTL'
                CALL catchFILE(NAME(8),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
                VAR='/PREC-CONV'
                CALL catchFILE(NAME(9),LDAS%catchDIR,YR2,MO2,DA2,HR2,
     &          LFB,count,var)
c	print *,'count=',count
 
c		  print*,'PRINTING catch NAME, next line'
c		  print*,NAME
		  ftype=1   !data will be catch
c		  print*,'Done with catch Name, Calling Retrive'
		  CALL RETcatch(ORDER,LDAS,GRID,NAME,FERROR,ftype,
     &           prevname,count)
c              print*,'Done with Retrieve of catch, FERROR=',FERROR
		  if(FERROR.eq.1) LDAS%catchTIME2=TIME2
c		  if(FERROR.eq.1) print*,'catch FILE READ in and APPROVED'
 		  FCA=0

c		 
c	print*,'got catch data, setting fsource 16 to 1' 
		  LDAS%FSOURCE(16)=1
		
		
	   ENDIF    ! IF FINDTIME2=1



	BTIME=LDAS%catchTIME1
	CALL TIME2DATE(BTIME,BDOY,GMT1,BYR,BMO,BDA,BHR,BMN)
        BTIME=LDAS%catchTIME2
        CALL TIME2DATE(BTIME,BDOY,GMT2,BYR,BMO,BDA,BHR,BMN)

!===  Interpolate Data in Time


         WT1=(LDAS%catchTIME2-LDAS%TIME)/(LDAS%catchTIME2-
     &  LDAS%catchTIME1)
	   WT2=1.0-WT1
c	print *,'wt1=',wt1
c	print *,'wt2=',wt2
	   DO F=1,9
	    IF(F.eq.3) then     !Shortwave
c	print *,'LDAS%GMT,GMT1,GMT2'
c	print *,LDAS%GMT,GMT1,GMT2

              DO C=1,LDAS%NC
                 DO R=1,LDAS%NR
                   ZDOY=LDAS%DOY

c	print *,'calling zterp for',c,r
c	print *,'GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,CZM'
                   CALL ZTERP(0,GRID(C,R)%LAT,GRID(C,R)%LON,
     &                GMT1,GMT2,LDAS%GMT,ZDOY,
     &                ZW1,ZW2,CZB,CZE,CZM,LDAS,GRID)
c	print *,'called zterp for',c,r
c        print *,'GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,CZM'
c        print *,GMT1,GMT2,LDAS%GMT,ZDOY,ZW1,ZW2,CZB,CZE,CZM



                   GRID(C,R)%FORCING(F)=GRID(C,R)%catchDATA1(F)*ZW1

                IF(GRID(C,R)%FIMASK.EQ.0) GRID(C,R)%FORCING(F)=LDAS%UDEF

	         if ((GRID(C,R)%FORCING(F).ne.LDAS%UDEF).and.
     &               (GRID(C,R)%FORCING(F).lt.0) ) then
                    print *,'2 warning!!!!  SW radiation is negative!!'
		    print *,'sw=',GRID(C,R)%FORCING(F),'...negative'
                    print *,'catch1 equaled',GRID(C,R)%catchDATA1(F)
                    print *,'forcing mask=',GRID(C,R)%FIMASK

                    stop
	         endif

	         if (GRID(C,R)%FORCING(F).gt.1367) then
		  print *,'warning, SW AVG RADIATION TOO HIGH!!'
                  print *,'it is',GRID(C,R)%FORCING(F),'at c,r',c,r
	          print *,'set = catchdata1',GRID(C,R)%catchDATA1(F)
                  print *,'catchdata2 was',GRID(C,R)%catchDATA2(F)
		  GRID(C,R)%FORCING(F)=GRID(C,R)%catchDATA1(F)
		 endif
                   GRID(C,R)%catchSW=GRID(C,R)%FORCING(F)
                 ENDDO
                ENDDO
		
	    ELSE IF(F.eq.8.or.F.eq.9) then    ! A precip variable Do Block Interpolation
	      DO C=1,LDAS%NC
		 DO R=1,LDAS%NR
		  GRID(C,R)%FORCING(F)=GRID(C,R)%catchDATA1(F)
		 ENDDO
		ENDDO
	    ELSE IF (F.eq.4) then     !Longwave
              DO C=1,LDAS%NC
               DO R=1,LDAS%NR
                GRID(C,R)%FORCING(F)=GRID(C,R)%catchDATA1(F)
               ENDDO
              ENDDO
	                ELSE IF(F.eq.5.or.F.eq.6) then    !Winds, u,v are equal 
                DO C=1,LDAS%NC
                 DO R=1,LDAS%NR
              GRID(C,R)%catchDATA1(F)=SQRT((GRID(C,R)%catchDATA1(F)*
     &                                GRID(C,R)%catchDATA1(F))/2.0)
              GRID(C,R)%catchDATA2(F)=SQRT((GRID(C,R)%catchDATA2(F)*
     &                                GRID(C,R)%catchDATA2(F))/2.0)
               GRID(C,R)%FORCING(F)=GRID(C,R)%catchDATA1(F)*WT1+
     &                           GRID(C,R)%catchDATA2(F)*WT2
                 ENDDO
                ENDDO

	
	    ELSE	   !Linearly interpolate everything else	
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
c	       print *,'grid',f,'equals',GRID(C,R)%catchDATA1(F)*WT1,
c     &  '+',GRID(C,R)%catchDATA2(F)*WT2
               GRID(C,R)%FORCING(F)=GRID(C,R)%catchDATA1(F)*WT1+
     &                           GRID(C,R)%catchDATA2(F)*WT2
             ENDDO
            ENDDO  	    
	     
	    ENDIF

	   ENDDO   !the F loop
	  
            DO C=1,LDAS%NC
             DO R=1,LDAS%NR
	IF (GRID(C,R)%FORCING(1).NE.LDAS%UDEF) THEN
c	print *,c,r,GRID(C,R)%FORCING(1)-273.15
               vp(c,r)=6.112*(10**((7.5*
     &         (GRID(C,R)%FORCING(1)-273.15))/(237.7+
     &            (GRID(C,R)%FORCING(1)-273.15))))
c	print*,c,r,GRID(C,R)%FORCING(2),GRID(C,R)%FORCING(7)/100.0,
c     &         vp(c,r)
               GRID(C,R)%FORCING(2)=(0.622*vp(c,r))/
     &         ((GRID(C,R)%FORCING(7)/100.0)-(0.378*vp(c,r)))

c	print *,c,r,GRID(C,R)%FORCING(2)
	ENDIF
             ENDDO
            ENDDO
!=== ADJUST PRECIP TO VALUE PER SEC DATA SHOULD COME IN AS VALUE PER
!=== TIME PERIOD THREE HOUR FOR catch	 
!--- POSSIBLE ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR  ERROR ERROR ERROR ERROR  ERROR ERROR ERROR ERROR 
         DO C=1,LDAS%NC
	    DO R=1,LDAS%NR
	     GRID(C,R)%FORCING(8)=GRID(C,R)%FORCING(8)/(6.0*60.0*60.0)
	     GRID(C,R)%FORCING(9)=GRID(C,R)%FORCING(9)/(6.0*60.0*60.0) 
	    ENDDO
	   ENDDO

C
84       format('NOW',i4,4i3,2x,'PVT ',a22,' NXT ',a22)
          if(FINDTIME2.eq.1) then
	      write(83,*) 'NXT-Name: ',NAME
            write(83,84) YR1,MO1,DA1,zzHR1,zzMN1,
     &     LDAS%EVT1,LDAS%EVT2
           endif
c
	  
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
!   This subroutine puts together catch file name
!
!=======================================================
	   
	   SUBROUTINE catchFILE(NAME,catchDIR,YR,MO,DA,HR,LFB,
     &  count,var)
	   
	   IMPLICIT  NONE
	   
!=== Local Variables
          CHARACTER*80 NAME
	    CHARACTER*40 catchDIR
	    INTEGER YR,MO,DA,HR,I,C,count,CC
	    
	    INTEGER LFB   ! 1 = look backward , 2 = look forward in time
	    
	    INTEGER UYR,UMO,UDA,UHR,UMN,USS,TA
	    INTEGER SUYR  !two digit year
	    REAL*8 TIME1,DUMTIME
	    REAL GMT1,DUMGMT
	    INTEGER DOY1,DUMDOY
	    CHARACTER*1 FNAME(80),FBASE(40),FSUBS(80),FBASEVAR(40)
	    CHARACTER*1 FTIME(6),FDIR(6)
	    CHARACTER*2 CHINSERT
	    CHARACTER*2  SPCODE(4)   !Special Code
	    CHARACTER*2  HRCODE
            CHARACTER*40 VAR
	    INTEGER PNTSPC  !points to correct Special Code
	    INTEGER IBREAK
	    
	    DATA SPCODE/'12','09','06','03'/
	    
!==== End Variable Definition
     
	IF ((DA-1).NE.0) THEN
c	print *,'day greater than 1'
	  COUNT=( (4*(DA-1)) + (HR/6) +1 )
	ELSE
c	print *,'day equals one',hr
          COUNT=(HR/6)+1
        ENDIF 
c	print *,'in filesub, count=',count
     
!=== Put together filename
90       FORMAT(80A1)
91       FORMAT(A4,i3,A11,I3)
92       FORMAT(40A1)
93       FORMAT(A80)
94       FORMAT(I4,I2)
95       FORMAT(6a1)
96       FORMAT(a40)
97       FORMAT(a12,a2,a3)
98       FORMAT(a1,I4,a1)
99       FORMAT(6a1)
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
c	   UHR=(HR/6)*6

C	Because of catch naming convention, subtract off three hours from
c	time in order to open the correct file.
C         print*,'Jared, HR=',HR,'  UHR=',UHR
C        print *,'before uyr,umo,uda,uhr,umn,uss',uyr,umo,uda,uhr,umn,uss
c         CALL TICK(DUMTIME,DUMDOY,DUMGMT,UYR,UMO,UDA,UHR,UMN,USS,
c     &    -21600)
C	print *,'after uyr,umo,uda,uhr,umn,uss',uyr,umo,uda,uhr,umn,uss 
C 
C
C  If the time is 12 or later the file is time stamped
C  with the next day.  So check for that first
 
c	print *,'at astart, var=',var
	   
         OPEN(90,FILE='temp',FORM='FORMATTED',
     &	   ACCESS='DIRECT',RECL=80)
c	print *,'writing catchdir',catchdir
         WRITE(90,96,REC=1) catchDIR
         READ(90,92,REC=1) (FBASE(I),I=1,40)

         WRITE(90,96,REC=1) var
         READ(90,92,REC=1) (FBASEVAR(I),I=1,40)
        CC=0
        DO I=1,40
          IF(FBASEVAR(I).EQ.(' ').AND.CC.EQ.0) CC=I-1
        ENDDO
c	print *,'fbasevar=',fbasevar,'var=',var
         if(UYR.ge.2000) SUYR=UYR-2000
	   if(UYR.le.1999) SUYR=UYR-1900
 
         WRITE(90,98,rec=1)'/',UYR,'/'
         READ(90,99,REC=1) FDIR
         DO I=1,6
          IF(FDIR(I).EQ.(' ')) FDIR(I)='0'
         ENDDO
C         IBREAK=10     !all files have 4 digit years now
c	print *,'writing'
	     WRITE(90,94,REC=1) UYR,UMO
           READ(90,95,REC=1) FTIME
           DO I=1,6
             IF(FTIME(I).EQ.(' ')) FTIME(I)='0'	     
           ENDDO	   
      	 
  
        C=0
        DO I=1,40
          IF(FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1 
        ENDDO
	   WRITE(90,90,REC=1) (FBASE(I),I=1,C),(FBASEVAR(I),I=1,CC),
     &  (FDIR(I),I=1,6),
     &  (FTIME(I),I=1,6)

       READ(90,93,REC=1) NAME

c	print *,'trying to open ',NAME
	 CLOSE(90)
	 RETURN
	 END     
         	    	   
	  
