	SUBROUTINE GRIBOUT_FORCE_GRIB2(ldas,grid,fbase,fyrmodir)


***********************************************************************
C  PROGRAM:  GRIBLDAS           WRITE LDAS OUTPUT IN GRIB FORMAT
C  PRGMMR:  MARSHALL & LOHMANN  ORG:  W/NP20
C
C  ABSTRACT:  TAKES HR1LY OUTPUT FILES FROM LDAS MODELS (IN BINARY
C             FORMAT - BIG ENDIAN) AND MAKES GRIB RECORDS FOR ALL FIELDS
C             IN THE FILE.  NOTE ORDER OF FIELDS IN EACH BINARY FILE
C             MUST CORRESPOND TO ORDER OF FIELDS IN KPDS.tbl.  TIME
C             STAMP INFORMATION MUST ALSO APPEAR (YYYYMMDDHH) IN THE
C             TITLE OF EACH BINARY HR1LY OUTPUT FILE.  THIS INFO IS
C             USED IN CONSTRUCTING THE TIME-RELATED PDS OCTETS WHICH
C             DO NOT APPEAR IN KPDS.tbl   
C
C  PROGRAM HISTORY LOG:
C  00-10-15 MARSHALL ORIGINAL TEMPLATE FOR NOAH MODEL
C  00-11-01 LOHMANN  FINALIZED AND CHANGED TO SUBROUTINE
c  00-11-16 Cosgrove Altered for use with NASA LDAS model driver
c  00-11-17 Cosgrove Added ability to output on time intervals other
c           than hourly.  ALSO, ADDED PRINT STATEMENT IN BEGINNING
C           TO PRINT BLANK LINE...THIS KEEPS LDAS DRIVER FROM 
C           CRASHING IN THIS SUBROUTINE....THERE APPARENTLY IS A
C           MEMORY-RELATED BUG IN THE DRIVER, BUT EXTENSIVE DEBUGGING
C           FAILED TO FIND IT AND THE PRINT STATMENT WAS ADDED AS A 
C           STOP-GAP MEASURE TO ALLOW THE LDAS DRIVER TO RUN
c  00-11-17 Cosgrove Added checks and fixes for specific humidity and precip
c           output
c  05 Feb 2000:  Brian Cosgrove; Changed directory structure of output, added
checks for CAPE 
c  02-05-15 Jambor Chianged LOGICAL to LOGICAL*1 to match new 
c           GRIB libraries
c 13-10-28  Youlong Xia; revised code to output grib2 format
C***********************************************************************

! Declare modules and data structures
      use ldas_module      ! ldas non-model-specific 1-d variables
      use grid_module      ! ldas non-model-specific grid variables
      implicit none
      type (ldasdec) ldas
      type (griddec) grid(ldas%nc,ldas%nr)

	logical found
	integer bound,dist,xx,yy,x,y,nx,ny
      INTEGER c,r,t,m,cc,jj,counter
      INTEGER SS1,TS,MN1,HR1,DA1,MO1,YR1,TS1,DOYT,DOYP
      INTEGER NSOLD,DOY1
      PARAMETER (NSOLD=8)
      CHARACTER*8 TODAY, YESTERDAY
      CHARACTER*1 TOD(8),YES(8)
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13),FTIMED(4)
      CHARACTER*1  FYRMODIR(80),FSUBFT(80),EXPCODE(3),FTIMEC(8)
      CHARACTER*1  FSUBFG(80),EXPARRAY(80),FTIMEB(10),FSUBGB(18)
      CHARACTER*1 GRIBF(256),GRIBF2(256)

      CHARACTER*13 VNAME 
      !FOR DEBUGGING
      CHARACTER*80  CZNAME
      CHARACTER*4   SYEAR
      CHARACTER*2   SMONTH,SDAY,SHOUR

      INTEGER luga, lugb,I, J, K, N, IRET, KPDS(13), 
     &        NLDAS, LENGDS, IG, JRET, NFIELDS
       real :: lat(ldas%nc,ldas%nr),lon(ldas%nc,ldas%nr)
       real :: czmgrid(ldas%nc,ldas%nr),radgrid(ldas%nc,ldas%nr)
       real :: gmtm, tmp, tmpcz, fill, doyd, doy
       real :: zw1, zw2, czb, cze, czm, undef
       real :: angle, deg2rad, maxrad, czmean
       integer rad,ifound,irad,i1,j1,ix,jx
       integer p_angle, c_angle, numarcs


      LOGICAL*1  LDASMASK(LDAS%NC,LDAS%NR)

      CHARACTER CENT*2, CHOUR*2, GDS(400)
      CHARACTER*256 GRIBFILEA,GRIBFILEB

      REAL DUMMYGMT,DUMMY(LDAS%NC,LDAS%NR),rho,vap,newtemp
      REAL DUMMY1(LDAS%NC,LDAS%NR),DUMMY2(LDAS%NC,LDAS%NR)
	REAL DUMMY0(LDAS%NC,LDAS%NR),DUMMY3(LDAS%NC,LDAS%NR)
      REAL*8 DUMMYTIME 
	NX=464
	NY=224
!	Construct GRIB logical type land/sea mask from LDAS integer land/sea mask

      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
          IF((GRID(C,R)%FIMASK.EQ.1).OR.(GRID(C,R)%FIMASK.EQ.2)
     &  .OR.(GRID(C,R)%FIMASK.EQ.3)) THEN
            LDASMASK(C,R)=.TRUE.
          ELSE
            LDASMASK(C,R)=.FALSE.
          ENDIF
       ENDDO
      ENDDO
        HR1=LDAS%HR
        DA1=LDAS%DA
        MO1=LDAS%MO
        YR1=LDAS%YR
        MN1=LDAS%MN
        SS1=0
	TS1=-3600*24
	DUMMYGMT=1.0
        DUMMYTIME=1.0

!	Construct Grib output file names
109    FORMAT(I4,I2,I2)
110    FORMAT(8A1)
111    FORMAT(I4,I2,I2,I2)
112    FORMAT(10A1)
113    FORMAT(I4)
114    FORMAT(4A1)

       OPEN(95,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
        WRITE(95,109,REC=1)YR1,MO1,DA1
        READ(95,110,REC=1)TOD
        DO I=1,8
         IF(TOD(I).EQ.(' '))TOD(I)='0'
        ENDDO
	TODAY=TOD(1)//TOD(2)//TOD(3)//TOD(4)//TOD(5)
     &  //TOD(6)//TOD(7)//TOD(8)

      CALL TICK(DUMMYTIME,DOY1,DUMMYGMT,YR1,MO1,DA1,HR1,MN1,SS1,TS1)

        WRITE(95,109,REC=1)YR1,MO1,DA1
        READ(95,110,REC=1)YES
        DO I=1,8
         IF(YES(I).EQ.(' '))YES(I)='0'
        ENDDO
        YESTERDAY=YES(1)//YES(2)//YES(3)//YES(4)//YES(5)
     &  //YES(6)//YES(7)//YES(8)

        DO I=1,13
         KPDS(I)=0
        ENDDO
C     read past header info in KPDS.tbl
      OPEN (UNIT = 30, FILE ='KPDS_completemosaic.force.narr.grib2.tbl')
      DO K = 1, 42
         READ(30,*)
      END DO

C     MAKE OUTPUT FILENAME A

 92    FORMAT(80A1)
 93    FORMAT(A256)
 91    FORMAT(256A1)
 94    FORMAT(A80)
 95    format(80A,A4,3A1,A5,4A1,A1,8A1,A1,10A1,15A1)
101    FORMAT(A18)
102    FORMAT (40I)
103    FORMAT (3I)

        WRITE(95,111,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
        READ(95,112,REC=1)FTIMEB
        DO I=1,10
         IF(FTIMEB(I).EQ.(' '))FTIMEB(I)='0'
        ENDDO

        WRITE(95,109,REC=1)LDAS%YR,LDAS%MO,LDAS%DA
        READ(95,110,REC=1)FTIMEC
        DO I=1,8
         IF(FTIMEC(I).EQ.(' '))FTIMEC(I)='0'
        ENDDO

        WRITE(95,113,REC=1)LDAS%YR
        READ(95,114,REC=1)FTIMED
        DO I=1,4
         IF(FTIMED(I).EQ.(' '))FTIMED(I)='0'
        ENDDO

        WRITE(95,101,REC=1)'.nldasforce-a.grb2'

        READ(95,92,REC=1) (FSUBGB(I),I=1,18)
        C=0
       DO I=1,40
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO

C      IF ( (HR1.LE. 23).AND.(HR1.GE.13) ) THEN
         WRITE (CHOUR, '(i2.2)') HR1
        WRITE(95,103,REC=1)LDAS%EXPCODE
c       WRITE(95,102,REC=1)LDAS%EXPCODE
        READ(95,92,REC=1)EXPARRAY
        COUNTER=1
        DO I=1,40
         IF(EXPARRAY(I).NE.(' ')) THEN
          EXPCODE(COUNTER)=EXPARRAY(I)
          COUNTER=COUNTER+1
         ENDIF
        ENDDO

        WRITE(95,95,REC=1)(FBASE(I),I=1,C),'/EXP',EXPCODE,'/FOR/',
     &  FTIMED,'/',FTIMEC,'/',
     &   FTIMEB,(FSUBGB(I),I=1,18)
        READ(95,93,REC=1) GRIBFILEA

C     MAKE OUTPUT FILENAME B

        WRITE(95,111,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
        READ(95,112,REC=1)FTIMEB
        DO I=1,10
         IF(FTIMEB(I).EQ.(' '))FTIMEB(I)='0'
        ENDDO

        WRITE(95,109,REC=1)LDAS%YR,LDAS%MO,LDAS%DA
        READ(95,110,REC=1)FTIMEC
        DO I=1,8
         IF(FTIMEC(I).EQ.(' '))FTIMEC(I)='0'
        ENDDO

        WRITE(95,113,REC=1)LDAS%YR
        READ(95,114,REC=1)FTIMED
        DO I=1,4
         IF(FTIMED(I).EQ.(' '))FTIMED(I)='0'
        ENDDO

        WRITE(95,101,REC=1)'.nldasforce-b.grb2'
        READ(95,92,REC=1) (FSUBGB(I),I=1,18)
        C=0
       DO I=1,40
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO

C      IF ( (HR1.LE. 23).AND.(HR1.GE.13) ) THEN
         WRITE (CHOUR, '(i2.2)') HR1
        WRITE(95,103,REC=1)LDAS%EXPCODE
c       WRITE(95,102,REC=1)LDAS%EXPCODE
        READ(95,92,REC=1)EXPARRAY
        COUNTER=1
        DO I=1,40
         IF(EXPARRAY(I).NE.(' ')) THEN
          EXPCODE(COUNTER)=EXPARRAY(I)
          COUNTER=COUNTER+1
         ENDIF
        ENDDO

        WRITE(95,95,REC=1)(FBASE(I),I=1,C),'/EXP',EXPCODE,'/FOR/',
     &  FTIMED,'/',FTIMEC,'/',
     &   FTIMEB,(FSUBGB(I),I=1,18)
        READ(95,93,REC=1) GRIBFILEB
        CLOSE (95)
C     OPEN GRIB FILE

      luga = 40
      lugb = 41
C --------open NLDAS a & b file --------------------------------------------
      CALL BAOPEN (luga, GRIBFILEA, IRET)
      CALL BAOPEN (lugb, GRIBFILEB, IRET)
C     WRITE EACH FIELD (IN GRIB) TO HR1LY OUTPUT FIELD
 15   FORMAT (29x, 7I6) 
      NLDAS = LDAS%NC * LDAS%NR
 
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(14)
         ENDDO
        ENDDO

      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     &DUMMY,KPDS,luga,IRET) 
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = TMP'
         STOP
      ELSE
C         PRINT*, 'wrote TMP'
      END IF

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(15)
         ENDDO
        ENDDO
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call gRIb2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SPFH'
         STOP
      ELSE
C         PRINT*, 'wrote SPFH'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(3)
         ENDDO
        ENDDO
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = PRESN'
         STOP
      ELSE
C         PRINT*, 'wrote PRESN'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(1)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = UGRD'
         STOP
      ELSE
C         PRINT*, 'wrote UGRD'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(2)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = VGRD'
         STOP
      ELSE
C         PRINT*, 'wrote VGRD'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(20)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DSWRF'
         STOP
      ELSE
C         PRINT*, 'wrote DSWRF'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(21)
         ENDDO
        ENDDO

      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DLWRF'
         STOP
      ELSE
C         PRINT*, 'wrote DLWRF'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%ETADATA1(4)/3.0
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = APCPN'
         STOP
      ELSE
C         PRINT*, 'wrote APCPN'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY1(I,J)=
     &  GRID(I,J)%ETADATA1(5)
          DUMMY2(I,J)=
     &  GRID(I,J)%ETADATA1(4)

	if (dummy1(i,j).gt.dummy2(i,j)) dummy1(i,j)=dummy2(i,j)

	if (dummy2(i,j).gt.0) then
          DUMMY(I,J)=
     &  dummy1(i,j)/dummy2(i,j)
	else
	  DUMMY(I,J)=0.0
	endif

          DUMMY1(I,J)=
     &  GRID(I,J)%ETADATA1(5)/3.0
          DUMMY2(I,J)=
     &  GRID(I,J)%ETADATA1(4)/3.0

        if (dummy1(i,j).gt.dummy2(i,j)) dummy1(i,j)=dummy2(i,j)

	DUMMY3(I,J)=GRID(I,J)%FORCING(4)*3600.0
        if (DUMMY3(I,J).le.0) dummy(i,j)=0.0
        if (DUMMY3(I,J).le.-999) dummy3(i,j)=-999
        dummy0(i,j)=dummy(i,j)

         ENDDO
        ENDDO

         !  REPLACE MISSING convective/total ratio POINTS WITH 
         !  NEAREST NEIGHBOR IF NEAREST NEIGHBOR
         !  WITHIN ABOUT 240 KM ("BOUND").  ELSE, REPLACE WITH Zero VALUE.

          rad=464
          deg2rad = 0.01745329

!        open (unit=45,file='out.bin',form='unformatted',
!     &          status='unknown')
!       write(45) dummy2
!       close(45)
!        stop

          DO y=1, ldas%nr
          DO x=1, ldas%nc
        if (dummy(x,y).gt.1.0) print *,'first dummy, dummy1,dummy2',
     &  dummy(x,y),dummy1(x,y),dummy2(x,y)

! out2 new and out is orig array
!        out2(i,j) = out(i,j)
!        IF (out(i, j) .EQ. undef) THEN
         IF ((GRID(X,Y)%FORCING(4).gt.0).and.
     &             (dummy2(x,y).le.0.0)) THEN

           ifound = 0
! search radius is rad
           DO irad = 1, rad
             numarcs = irad * 8
             p_angle = 360000 / numarcs
             DO c_angle=0,360000-p_angle,p_angle
               angle = c_angle / 1000.00
               j1 = y + INT(COS(angle*deg2rad)*irad)
               i1 = x + INT(SIN(angle*deg2rad)*irad)
               jx = j1
               ix = i1
               IF (jx .LE. 0 ) jx = 1
               IF (jx .GE. ldas%nr ) jx = ldas%nr
               IF (ix .LE. 0 ) ix = 1
               IF (ix .GT. ldas%nc ) ix = ldas%nc
               !write(*,'(f8.2,5i4)') angle,i1,j1,ix,jx,irad
!               IF(out(ix, jx) .NE. undef ) then
               IF (dummy2(ix,jx).gt.0.0) THEN
               tmp = dummy1(ix,jx)/dummy2(ix,jx)
	        if (tmp.gt.1.0) print *,'tmp, dummy1,dummy2',
     &  tmp,dummy1(ix,jx),dummy2(ix,jx)

                 ifound = ifound + 1
               ENDIF
! number of points to find
               IF(ifound .GE. 1) EXIT
             ENDDO
             IF (ifound .GE. 1) EXIT
           ENDDO
           IF (ifound .GE. 1) THEN
             dummy(X,Y) = tmp
           ELSE
              dummy(x,y)=1000.0
           ENDIF
         ELSE
!         dummy(x,y)=dummy2(x,y)

         ENDIF !IF MISSING DATA AT POINT

        if (dummy(x,y).gt.1.0) print *,'dummy, dummy1,dummy2',
     &  dummy(x,y),dummy1(x,y),dummy2(x,y)

          ENDDO !NC
         ENDDO !NR
 
       READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',',
     &  'FIELD = ACPCP Fraction'
        STOP
      ELSE
C         PRINT*, 'wrote ACPCP FRACTION'
      END IF

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%ETADATA1(5)/3.0
         ENDDO
        ENDDO

      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ACPCP'
         STOP
      ELSE
C         PRINT*, 'wrote ACPCP'
      END IF

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(6)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = CAPE'
         STOP
      ELSE
C         PRINT*, 'wrote CAPE'
      END IF

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%ETADATA1(19)/3.0
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = PEVAP'
         STOP
      ELSE
C         PRINT*, 'wrote PEVAP'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
         ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(17)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET) 
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SFEXC'
         STOP
      ELSE
C         PRINT*, 'wrote SFEXC'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(7)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = Hybrid TMP'
         STOP
      ELSE
C         PRINT*, 'wrote Hybrid TMP'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(10)
         ENDDO
        ENDDO
     
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET) 
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = 
     &  ',HR1,',','FIELD = Hybrid SPFH'
         STOP
      ELSE
C         PRINT*, 'wrote Hybrid SPFH'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(8)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = 
     &  ',HR1,',','FIELD = Hybrid PRES'
         STOP
      ELSE
C         PRINT*, 'wrote Hybrid PRES'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(11)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = 
     &  ',HR1,',','FIELD = Hybrid UGRD'
         STOP
      ELSE
C         PRINT*, 'wrote Hybrid UGRD'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(12)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = 
     &  ',HR1,',','FIELD = Hybrid VGRD'
         STOP
      ELSE
C         PRINT*, 'wrote Hybrid VGRD'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(9)
         ENDDO
        ENDDO
      
      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,lugb,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = Hybrid HGT'
         STOP
      ELSE
C         PRINT*, 'wrote Hybrid HGT'
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=GRID(I,J)%FORCING(4)*3600.0
         ENDDO
        ENDDO

      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',',
     &  'FIELD = Merged Precip APCP'
         STOP
      ELSE
C         PRINT*, 'wrote Merged Precip APCP' 
      END IF
        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DUMMY(I,J)=LDAS%UDEF
         ENDDO
        ENDDO

        !INITIALIZE TIME DATA FOR ZTERP
        HR1=LDAS%HR
        DA1=LDAS%DA
        MO1=LDAS%MO
        YR1=LDAS%YR
        MN1=LDAS%MN
        SS1=0
	TS1=0
	DUMMYGMT=1.0
        DUMMYTIME=1.0
        DOYT = 0
        DOYP = 0

        !CALL TICK TO GET DOY AND GMT TIME FOR ZTERP
        CALL TICK(DUMMYTIME,DOYP,DUMMYGMT,YR1,MO1,DA1,HR1,MN1,SS1,TS1)

!        OPEN(UNIT=80,FILE=CZNAME,FORM="UNFORMATTED")

        !PRINT*,"TDATE1: ",DOYP,DUMMYGMT,DUMMYTIME,CZNAME
        !PRINT*,"TDATE2: ",YR1,MO1,DA1,HR1,MN1,SS1,TS1

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
          DOYT = DOYP
          !PRINT*,"TINFO1: ",DOYT,DOYP,DUMMYGMT,CZB,CZM,CZE
          CALL ZTERP(1,GRID(I,J)%LAT,GRID(I,J)%LON,DUMMYGMT,
     +      DUMMYGMT,DUMMYGMT,DOYT,ZW1,ZW2,CZB,CZE,CZM,LDAS,
     +      GRID,CZMEAN)
          !PRINT*,"TINFO2: ",DOY,DOY1,DUMMYGMT,CZB,CZM,CZE
          !CZMGRID(I,J) = CZM
          DUMMY(I,J)=GRID(I,J)%FORCING(22)*GRID(I,J)%FORCING(20)
          MAXRAD = 1367.0 * CZM
          ! Allow 1% over MAXRAD for diffuse cases
          MAXRAD = MAXRAD + MAXRAD*(0.01) 
          !RADGRID(I,J) = MAXRAD
          IF( DUMMY(I,J) .GT. MAXRAD ) THEN
            PRINT*,"BCRAD TOO HIGH (X,Y,VALUE): ",I,J,DUMMY(I,J)
            PRINT*,"CORRECTED TO: ",MAXRAD
            DUMMY(I,J) = MAXRAD
          ENDIF
         ENDDO
        ENDDO

      !WRITE(80) CZMGRID
      !CLOSE(80)

      READ (30, *) VNAME, (KPDS(I), I=1,13)
      call grib2_wrt_g2func(LDASMASK,HR1,LDAS%WRITEINTM,YESTERDAY,TODAY,
     & DUMMY,KPDS,luga,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DSWRF-BC'
         STOP
      ELSE
C         PRINT*, 'wrote DSWRF-BC'
      END IF

      CALL BACLOSE (luga, JRET)
      CALL BACLOSE (lugb, JRET)
      CLOSE (30)
      END
