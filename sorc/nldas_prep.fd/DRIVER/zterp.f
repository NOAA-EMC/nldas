
        SUBROUTINE ZTERP (IFLAG,LAT,LON,BTIME,ETIME,
     &          MBTIME,JULIANB,weight1,weight2,czbegdata,
     &          czenddata,czmodel,LDAS,GRID,avgangle)


**********************************************************************************
*
*       PROGRAMMER: Brian Cosgrove    ORGANIZATION: NASA/GSFC   DATE: 10/2/98
*
*
*       ABSTRACT:
*
*       This subroutine is based, in part, on modified subroutines from
*       Jean C. Morrill of the GSWP project.  The program temporally interpolates
*       time averge or instantaneous data to that needed by the model
*       at the current timestep.  It does this by combining a linear interpolation
*       approach with a solar zenith angle approach in a fashion suitable for use
*       with data such as short wave radiation values.  
*       It cannot be used with input data points which are more
*       than 24 hours apart.  The program outputs two weights which can then be
*       applied to the original data to arrive at the interpolated data.  If
*       IFLAG=0, then WEIGHT1 is the weight which should be applied to the
*       time averaged data (from the time period which the model is currently in)
*       to arrive at the interpolated value and weight 2 is not used at all.  If
*       IFLAG=1, then WEIGHT1 should be applied to the original instantaneous
*       data located just prior to the model time step, and WEIGHT2 should be
*       applied to the original instantaneous data located just after the model
*       time step.
*       i.e.  (IF IFLAG=0)   interpolated data = (WEIGHT1 * time averaged data from
*                                                time period that model is currently in)
*       i.e.  (IF IFLAG=1) interp. data = (WEIGHT1*past data)+(WEIGHT2*future data)
**********************************************************************************
*
*
*       PROGRAM HISTORY LOG:
*        10/2/98 Brian Cosgrove
*        6/28/00 Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
*                not a hard-wired undefined value of -999.999.  Also changed
*                call to ZTERP subroutine so that LDAS and GRID mod files
*                are brought in and can be accessed in this subroutine
*        2/27/01 Brian Cosgrove; Added czmodel into call for ZTERP subroutine
*        10/7/01 Urszula Jambor; Added conditional to prevent stop in code when
*                AGRMET data are available, but both endtime cos(zen.) are small.
*
**********************************************************************************
*
*       INPUT ARGUMENT LIST:
*
*       IFLAG   -Int., flag specifing whether input data is
*               time averaged (IFLAG=0) or
*               instantaneous (IFLAG=1)
*       LAT     -Real, latitude (deg) of current data point
*       LON     -Real, longitude (deg) of current data point
*       BTIME   -Real, beginning time of orig. avg. data (IFLAG=0) or
*               time of original instantaneous data point which
*               is located at or just prior to the current model
*                time step (IFLAG=1).  Expects GMT time in hour
*               fractions.  i.e., 6:30 Z would be 6.5.
*       ETIME   -Real, endding time of orig. avg. data (IFLAG=0) or
*               time of original instantaneous data point which
*               is located at or just after the current model
*               time step (IFLAG=1).    Expects GMT time in hour
*               fractions.  i.e., 6:30 Z would be 6.5.
*       MBTIME  -Real, time of current model time step.
*               Expects GMT time in hour fractions
*               i.e., 6:30 Z would be 6.5.
*       JULIANB -Int., Julian day upon which BTIME falls
*
*
**********************************************************************************
*
*       OUTPUT ARGUMENT LIST:
*
*       WEIGHT1 -Real, weight applied to original time averaged
*               data (IFLAG=0) or
*               weight applied to orig instantaneous data point
*               located just prior to the current model time step
*       WEIGHT2 -Real, weight applied to orig instantaneous
*               data point located just after the current model
*               time step (IFLAG=1)
*               If IFLAG=0, then this weight is meaningless and
*               should not be used
*
**********************************************************************************
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS
      type (griddec) GRID(LDAS%NC,LDAS%NR)

        REAL LAT,LON,MBTIME,INTERVAL
        REAL BTIME,ETIME,WEIGHT1,WEIGHT2
        REAL CZMODEL,CZAVGDATA,TOTANGLE,AVGANGLE
        REAL CZENDDATA,CZBEGDATA,LHOUR,GMT
        REAL LMBTIME,LETIME,LBTIME,WEIGHTE,WEIGHTB
        REAL DIFFEB,DIFFBM,DIFFEM
        INTEGER ZONEETIME,ZONEBTIME,ZONEMBTIME,IFLAG
        INTEGER ZONE,JULIANB,JULIANE,JULIANMB,JULIANTEMP,I




C       This section contains hardwired data that will be supplied by main program.
C       These values were chosen arbitrarily and exist simply to check the
C       functioning of the program.
C
C       Initialize variables
C
        I=1
        TOTANGLE=0
        WEIGHTE=LDAS%UDEF
        WEIGHTB=LDAS%UDEF
        WEIGHT1=LDAS%UDEF
        WEIGHT2=LDAS%UDEF
        CZBEGDATA=LDAS%UDEF
        CZENDDATA=LDAS%UDEF
        CZMODEL=LDAS%UDEF
        GMT=BTIME
        JULIANE=JULIANB
        JULIANMB=JULIANB
        JULIANTEMP=JULIANB

        IF (MBTIME.LT.BTIME) JULIANMB=JULIANMB+1
        IF (ETIME.LE.BTIME) JULIANE=JULIANE+1

C       First case, IFLAG=0 (Time average input, instantaneous output)
C
C       Compute time interval, here arbitrarily divided into 36 parts
        IF (IFLAG.EQ.0) THEN
                CALL LOCALTIME (MBTIME,LON,LHOUR,ZONE)
                CALL COSZENITH(LON,LAT,LHOUR,ZONE,JULIANMB,CZMODEL)
		IF (CZMODEL.EQ.0) THEN
			WEIGHT1=0
			RETURN
		ENDIF

                IF (ETIME.GT.BTIME) THEN
                        INTERVAL = ((ETIME-BTIME)/36.0)
		ELSEIF (ETIME.LT.BTIME) THEN
                        INTERVAL = (((24-BTIME)+(ETIME))/36.0)
                ELSE 
			INTERVAL=24.0/36.0
		ENDIF

C       Compute cosine of zenith angle for each time interval
C
                DO WHILE (I.LE.37)

                        IF (I.GT.1) THEN
                          IF ((GMT+INTERVAL).LT.24) THEN
                                  GMT=GMT+INTERVAL
			  ELSE
                                  GMT=(INTERVAL-(24-GMT))
                                  JULIANTEMP=JULIANTEMP+1
                          ENDIF
                        ENDIF
                        CALL LOCALTIME (GMT,LON,LHOUR,ZONE)
                        CALL COSZENITH(LON,LAT,LHOUR,ZONE,
     &                     JULIANTEMP,CZAVGDATA)
                        TOTANGLE=TOTANGLE+CZAVGDATA
                        I=I+1
                ENDDO
	
C       Compute average cosine of zenith angle and also
C       weight which will be applied to original data (WEIGHT1)
C
                AVGANGLE=(TOTANGLE/37.0)
                IF (AVGANGLE.EQ.0) THEN
                        WEIGHT1=0
                ELSE
                WEIGHT1=(CZMODEL/AVGANGLE)
                ENDIF

        ENDIF

C       Second case:  IFLAG=1 (instantaneous input and output)
        IF (IFLAG.eq.1) THEN
C
C       Compute local times and cosine (zenith angle)
C
                CALL LOCALTIME (BTIME,LON,LBTIME,ZONEBTIME)
                IF (LBTIME.GT.BTIME) JULIANB=JULIANB-1
                CALL LOCALTIME (ETIME,LON,LETIME,ZONEETIME)
                IF (LETIME.GT.ETIME) JULIANE=JULIANE-1
                CALL LOCALTIME (MBTIME,LON,LMBTIME,ZONEMBTIME)
                IF (LMBTIME.GT.MBTIME) JULIANMB=JULIANMB-1
                CALL COSZENITH (LON,LAT,LBTIME,ZONEBTIME,
     &              JULIANB,CZBEGDATA)
                CALL COSZENITH (LON,LAT,LETIME,ZONEETIME,
     &              JULIANE,CZENDDATA)
                CALL COSZENITH (LON,LAT,LMBTIME,ZONEMBTIME,
     &              JULIANMB,CZMODEL)

C       Decision tree to deal with contingencies
C       If COS(zenith angle at current model time =0, weight =0
C       If COS(zenith angle =0 at beg. and end times, PROBLEM, STOP
C       Otherwise use beginning and ending data to calculate weight
C
                IF (CZMODEL.LE.0.01) THEN
                        WEIGHT1=0
                        WEIGHT2=0
                ELSE
                IF ((CZBEGDATA.GT.0.01).or.(CZENDDATA.GT.0.01)) THEN
                        IF (CZBEGDATA.LE.0.01) THEN
                                WEIGHT1=0
                                WEIGHT2=(CZMODEL/CZENDDATA)
                        ENDIF
                        IF (CZENDDATA.LE.0.01) THEN
                                WEIGHT1=(CZMODEL/CZBEGDATA)
                                WEIGHT2=0
                        ENDIF
                        IF((CZENDDATA.GT.0.01).and.
     1                   (CZBEGDATA.GT.0.01))THEN

                        IF (BTIME.LE.MBTIME) THEN
                                DIFFBM=MBTIME-BTIME
                        ELSE
                                DIFFBM=24-BTIME+MBTIME
                        ENDIF

                        IF (ETIME.GE.MBTIME) THEN
                                DIFFEM=ETIME-MBTIME
                        ELSE
                                DIFFEM=24-MBTIME+ETIME
                        ENDIF

                        IF (ETIME.GT.BTIME) THEN
                                DIFFEB=ETIME-BTIME
                        ELSEIF (ETIME.EQ.BTIME) THEN

                                DIFFEB=24
                        ELSE
                                DIFFEB=24-BTIME+ETIME
                        ENDIF
                        WEIGHTE=(DIFFBM/DIFFEB)
                        WEIGHTB=(DIFFEM/DIFFEB)

                        WEIGHT1=((CZMODEL/CZBEGDATA)*WEIGHTB)
                        WEIGHT2=((CZMODEL/CZENDDATA)*WEIGHTE)

                        ENDIF
                ELSE IF (LDAS%AGRMETSW .EQ. 1) THEN
!                        WRITE(79,*) 'CZB & CZE: ',CZBEGDATA,CZENDDATA
!                        WRITE(79,*) 'Need to ignore AGRMET data'
                        print*,'Ignore AGRMET SW, fall back to model SW' 
                ELSE
                        WRITE(79,*)'NO DATA TO INTERPOLATE TO/FROM'
                        WRITE(79,*)'BEGINNING AND ENDING DATA BOTH = 0'
                        WRITE(79,*)'STOPPING!!!'

                        print*, 'NO DATA TO INTERPOLATE TO/FROM'
                        print*, 'BEGINNING AND ENDING DATA BOTH = 0'
                        print*, 'STOPPING!!!'                        

                        STOP
                ENDIF
                ENDIF




        ENDIF
 
	RETURN
        END





C---------------------------------------------------------------------C
      SUBROUTINE COSZENITH (LON,LATD,LHOUR,ZONE,JULIAN,CZENITH)
C
C    
C     The purpose is to calculate the following:
C
C        1)  Day angle (GAMMA)
C
C        2)  Solar DEClination
C
C        3)  Equation of time
C
C        4)  Local apparent time
C
C        5)  Hour angle
C
C        6)  Cosine of zenith angle
C
C     All equations come from "An Introduction to
C     Solar Radition" By Muhammad Iqbal, 1983.
C
C
C-----------------------------------------------------------------------
      implicit none
C---------------------------Local variable------------------------------
C
      INTEGER
     $     ZONE,                ! time zone (1-24) GMT=12
     $     JULIAN               ! julian day

C
      REAL
     $     CZENITH,             ! cosine of zenith angle (radians)
     $     DECd,                ! solar declination (degrees)
     $     DEC,                 ! solar declination (radians)
     $     et,                  ! equation of time (minutes)
     $     GAMMA,               ! day angle (radians)
     $     LATime,              ! local apparent time
     $     LCORR,               ! longitudical correction
     $     LHOUR,               ! local standard time
     $     LON,                 ! local longitude (deg)
     $     LLAT,                ! local latitude in radians
     $     LATD ,               ! local latitude in degrees
     $     LS,                  ! standard longitude (deg)
     $     OMEGAD,              ! omega in degrees
     $     OMEGA ,              ! omega in radians
     $     PI,                  ! universal constant PI [-]
     $     ZENITH,              ! zenith angle(radians)
     $     ZEND                 ! zenith angle(degress)
C
C     Neither ZENITH nor ZEND are necessary for this program.
C     I originally used them as checks, and left them here in
C     case anyone else had a use for them.
************************************************************************
C
C     1)  Day angle GAMMA (radians) page 3

      PI= 3.141592              ! universal constant PI
      GAMMA=2*PI*(JULIAN-1)/365.
************************************************************************
C     2) Solar declination (assumed constant for a 24 hour period)  page 7
C     in radians
C
      DEC=(0.006918-0.399912*COS(GAMMA)+0.070257*SIN(GAMMA)
     $     -0.006758*COS(2*GAMMA)+0.000907*SIN(2*GAMMA)
     $     -0.002697*COS(3*GAMMA)+0.00148*SIN(3*GAMMA))
      DECd=DEC*(180./PI)
C
C     maximum error 0.0006 rad (<3'), leads to error of less than 1/2 degree
C     in ZENITH angle
************************************************************************^M
C     3)  Equation of time  page 11

      et=(0.000075+0.001868*COS(GAMMA)-0.032077*SIN(GAMMA)
     $     -0.014615*COS(2*GAMMA)-0.04089*SIN(2*GAMMA))*229.18
C    
************************************************************************^M
C     4) Local apparent time  page 13
C
C     LS     standard longitude (nearest 15 degree meridian)
C     LON     local longitude
C     LHOUR  local standard time
C     LATIME local apparent time
C     LCORR  longitudunal correction (minutes)
C
      LS=((ZONE-1)*15)-180.
      LCORR=4.*(LS-LON)*(-1)
      LATIME=LHOUR+LCORR/60.+et/60.
      IF (LATIME.LT.0.) LATIME=LATIME+24
      IF (LATIME.GT.24.) LATIME=LATIME-24
************************************************************************
C     5) Hour angle OMEGA  page 15
C
C     hour angle is zero at noon, postive in the morning
C     It ranges from 180 to -180
C
C     OMEGAD is OMEGA in degrees, OMEGA is in radians
C
      OMEGAD=(LATime-12.)*(-15.)
      OMEGA=OMEGAD*PI/180.
C
************************************************************************
C     6)  Zenith angle page 15
C
C     CZENITH cosine of zenith angle (radians)
C     LATD=local latitude in degrees
C     LLAT=local latitude in radians
C
      LLAT=LATD*PI/180.
      CZENITH=SIN(DEC)*SIN(LLAT)+COS(DEC)*COS(LLAT)*COS(OMEGA)
      CZENITH=AMAX1(0.,CZENITH)
      ZENITH=ASIN(CZENITH)
      ZEND=180.*ZENITH/PI
************************************************************************
C
      RETURN
      END




        SUBROUTINE LOCALTIME (GMT,LON,LHOUR,ZONE)


        REAL
     $  GMT,                 ! GMT time (0-23)
     $  LON,                 ! longitude in degrees
     $  CHANGE,              ! the change in number of hours between
     $  LHOUR                ! local hour (0-23) 0= midnight, 23= 11:00 p.m.

        INTEGER
     &  I,                   ! working integer
     $  JULIAN,              ! day of the year (1-365)
     $  ZONE                 ! time zone (1-24)


C
C     Determine into which time ZONE (15 degree interval) the
C     longitude falls.
C
      DO I=1,25
         IF (LON.LT.(-187.5+(15*i))) THEN
            ZONE=I
            IF (ZONE.eq.25) ZONE=1
            EXIT
         END IF
      END DO
C
C     Calculate change (in number of hours) from GMT time to
C     local hour.  Change will be negative for zones < 13 and
C     positive for zones > 13.
C
C     There is also a correction for LHOUR < 0 and LHOUR > 23
C     to LHOUR between 0 and 23.
C
      IF (ZONE.LT.13) THEN
         CHANGE=ZONE-13
         LHOUR=GMT+CHANGE
      ELSEIF (ZONE.eq.13) THEN
         LHOUR=GMT
      ELSE
         CHANGE=ZONE-13
         LHOUR=GMT+CHANGE
      END IF
      IF (LHOUR.LT.0) LHOUR=LHOUR+24
      IF (LHOUR.GT.23) LHOUR=LHOUR-24
      RETURN
      END
C







