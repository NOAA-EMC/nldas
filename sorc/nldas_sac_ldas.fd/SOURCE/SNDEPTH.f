C **  SUBROUTINE CALCULATES SNOW DEPTH
C
      SUBROUTINE SNDEPTH(WE,SLIQ,DFALL,SGSLOS,SRFRZ,
     +                   TA,DTA,IDT,SH,DS,TSNOW)

C **  WE     - WATER EQUIVALENT, MM
C **  SLIQ   - LIQUID WATER CONTENT, MM
C **  DFALL  - NEW SNOWFALL, MM
C **  SGSLOS - GROUND SNOW MELT, MM
C **  SRFRZ  - REFROZEN MELT/RAIN WATER, MM
C **  TA     - AIR TEMPERATURE, C
C **  DTA    - AIR TEMPERATURE CHANGE FOR THE 6HR TIME INTERVAL  
C **  IDT    - TIME STEP, HR
C **  SH     - SNOW DEPTH, CM
C **  DS     - SNOW DENSITY, G/CM3
C **  TSNOW  - AVERAGE SNOW TEMPERATURE, CELSIUS

      DT=IDT
C  ADJUST SNOW DENSITY DUE TO SNOWFALL
      DHC=0.
      SDN=0.0
      TSNEW=TA
      IF(DFALL .GT. 0.) THEN
C  CALCILATE NEW SNOW FALL DEPTH/DENSITY
       CALL SNEW(TA,DFALL,DHC,SDN)
       CALL SNOWT(DHC,SDN,DFALL,0.,TA,DTA,TSNEW,0.)
       CALL SNOWPACK(DFALL,DT,DHC,SDN,0.,0.,0.,TSNEW)
      ENDIF
       
      IF(SH .GT. 0.0001) THEN
C  CALCULATE OLD SNOW COMPACTION/METAMORPHISM
       SHN=DHC+SH
       DSN=(SDN*DHC+DS*SH)/SHN
       CALL SNOWT(SHN,DSN,WE,SLIQ,TA,DTA,TSNOW,DHC)
       CALL SNOWPACK(WE,DT,SH,DS,SLIQ,DFALL,SRFRZ,TSNOW)
C  ACCOUNT FOR GROUND SNOW MELT
       IF(SGSLOS .GT. 0.) SH=SH-0.1*SGSLOS/DS
       IF(SH .LT. 0.) SH=0.0
C  COMBINE NEW SNOW FALL AND OLD SNOWPACK
       TSNOW=(TSNOW*SH+TSNEW*DHC)/(SH+DHC)
       SH=SH+DHC
       DS=0.1*WE/SH
      ELSE
C  THERE WAS NO SNOWPACK BEFORE THIS SNOW FALL
       SH=DHC
       DS=SDN
       IF(SGSLOS .GT. 0.) SH=SH-0.1*SGSLOS/DS
       IF(SH .LT. 0.) SH=0.0
       TSNOW=TSNEW
      ENDIF 
C      print*,'SNDPT,DS,TSNOW,WE,SLIQ,DFALL,SRFRZ,TA,DTA=',SH,DS,TSNOW,
C     &        WE,SLIQ,DFALL,SRFRZ,TA,DTA
      
      RETURN
      END
