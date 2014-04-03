	SUBROUTINE GRIBOUT_GLBFORCE (ldas,grid,fbase,fyrmodir)

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
c  05 Feb 2000:  Brian Cosgrove; Changed directory structure of output, added checks
c                for CAPE 
c  02-02-06 Jambor Renamed gribout_glbforce for use with global domains.
c           Changed output variables to correspond with 9 forcing 
c           variables.  Allow for different output time intervals.  
c           Modified allowed range for output parameters.
c           Need to find appropriate kpds(2) process ID codes for
c           GEOS and observed precipitation sources.
c  02-05-15 Jambor Changed LOGICAL to LOGICAL*1 to match new 
c           GRIB libraries
c  30-07-02 Gottschalck altered if/endif structures with multiple LDAS%GPCPSRC instead of GPCP
C***********************************************************************

! Declare modules and data structures
	use ldas_module		! ldas non-model-specific 1-d variables
	use grid_module		! ldas non-model-specific grid variables
	implicit none
	type (ldasdec) ldas
	type (griddec) grid(ldas%nc,ldas%nr)

	INTEGER c,r, COUNTER
	INTEGER SS1,TS,MN1,HR1,DA1,MO1,YR1,TS1,DOY1
	CHARACTER*8 TODAY, YESTERDAY
	CHARACTER*1 TOD(8),YES(8)
	CHARACTER*1 FBASE(40),FYRMODIR(80)
	CHARACTER*1 EXPCODE(3),EXPARRAY(80)
	CHARACTER*1 FTIMEB(10),FTIMEC(8),FTIMED(4),FSUBGB(14)

	INTEGER LUGB, I, J, K, IRET, origKPDS2, KPDS(25), KGDS(22)
	INTEGER NLDAS, LENGDS, IG, JRET

	LOGICAL*1  LDASMASK(LDAS%NC,LDAS%NR)
	
	CHARACTER CHOUR*2, GDS(400)
	CHARACTER*256 GRIBFILE

	REAL DUMMYGMT,DUMMY(LDAS%NC,LDAS%NR),rho,vap,newtemp
	REAL*8 DUMMYTIME 

!	Construct GRIB logical type land/sea mask from LDAS integer land/sea mask

	DO C=1,LDAS%NC
	   DO R=1,LDAS%NR
	      IF((GRID(C,R)%FIMASK.EQ.1).OR.(GRID(C,R)%FIMASK.EQ.2)
	1	   .OR.(GRID(C,R)%FIMASK.EQ.3)) THEN
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
 109	FORMAT(I4,I2,I2)
 110	FORMAT(8A1)
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
	1    //TOD(6)//TOD(7)//TOD(8)

	CALL TICK(DUMMYTIME,DOY1,DUMMYGMT,YR1,MO1,DA1,HR1,MN1,SS1,TS1)

        WRITE(95,109,REC=1)YR1,MO1,DA1
        READ(95,110,REC=1)YES
        DO I=1,8
	   IF(YES(I).EQ.(' '))YES(I)='0'
        ENDDO
        YESTERDAY=YES(1)//YES(2)//YES(3)//YES(4)//YES(5)
	1    //YES(6)//YES(7)//YES(8)
        DO I=1,25
	   KPDS(I)=0
        ENDDO
        DO I=1,22
	   KGDS(I)=0
        ENDDO

C     NON-CHANGING PDS ELEMENTS
	KPDS(1)=222		!ID FOR PRODUCTS
	if (LDAS%FORCE==1) then
	   KPDS(2)=082		!ID FOR GDAS MODEL
	else if (LDAS%FORCE==2) then
	   KPDS(2)=000		!ID FOR GEOS MODEL?? check with DAO
	else if (LDAS%FORCE==7) then
	   KPDS(2)=200		!ID FOR ECMWF MODEL
	else
	   KPDS(2)=000          !no ID
	   write(*,*) 'Invalid global model forcing'
	end if
	origKPDS2=KPDS(2)
	KPDS(4)=192		!BMS FLAG... DOn't worry about this.
	KPDS(12)=0		!assume output time minute always = 0
	KPDS(13)=1		!forecast time unit (hours)
	KPDS(17)=INT((LDAS%WRITEINTF*3600)/LDAS%TS) !number of time steps 
				!in averaged/accum variables
	KPDS(18)=0		!GRIB version -- left as 0 in NCEP products
	KPDS(19)=1		!version number of KPDS.tbl for LDAS.  
	KPDS(20)=0		!none missing from averages/accumulations
	KPDS(23)=222		!GSFC ID# 222
	KPDS(24)=0		!Does not apply to LDAS output
	KPDS(25)=0		!not used

C     DEFINE THE GDS FOR THE LDAS GRID (THIS ARRAY NEVER CHANGES FROM
C     FIELD TO FIELD. ONLY HAS TO BE DEFINED ONCE, UNLIKE THE PDS).

	IF (LDAS%DOMAIN==1) THEN
	   IG = 110		!NLDAS GRID, SUBROUTINE w3fi71.f in w3lib   
	   CALL MAKGDS(IG, KGDS, GDS, LENGDS, IRET)
	ELSE
	   DO J=1,22
	      KGDS(J)=LDAS%LDAS_KGDS(J)
	   ENDDO
	ENDIF

C     read past header info in KPDS.tbl
	OPEN (UNIT = 30, FILE = './SRC/KPDS_glbforce.tbl')
	DO K = 1, 43
	   READ(30,*)
	END DO

C     MAKE OUTPUT FILENAME
 92	FORMAT(80A1)
 93	FORMAT(A256)
 91	FORMAT(256A1)
 94	FORMAT(A80)
 95    format(80A,A4,3A1,A5,4A1,A1,8A1,A1,10A1,18A1)
 101	FORMAT(A12)
 102	FORMAT (40I)
 103	FORMAT (3I)

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

c        WRITE(95,101,REC=1)'.lsmforce_nasa'
        WRITE(95,101,REC=1)'.FORCING.GRB'

        READ(95,92,REC=1) (FSUBGB(I),I=1,12)
        C=0
	DO I=1,40
	   IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
	ENDDO

	WRITE (CHOUR, '(i2.2)') HR1
        WRITE(95,103,REC=1)LDAS%EXPCODE
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
     &   FTIMEB,(FSUBGB(I),I=1,12)
        READ(95,93,REC=1) GRIBFILE
        CLOSE (95)

C     OPEN GRIB FILE

	LUGB = 40
	CALL BAOPEN (LUGB, GRIBFILE, IRET)

C     WRITE EACH FIELD (IN GRIB) TO OUTPUT FILE
 15	FORMAT (29x, 7I6) 

	NLDAS = LDAS%NC * LDAS%NR

        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(1)
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF ((DUMMY(I,J).GT.353.0).OR.(DUMMY(I,J).LT.213.0)) THEN
		    PRINT*,'CORRECTING, VALUE OUT OF RANGE,
	1		 TEMP(',I,',',J,')=',DUMMY(I,J)
		    newtemp=SQRT(SQRT(GRID(I,J)%FORCING(4)/5.67E-08))
		    dummy(i,j)=newtemp
		    print *,'corrected temp value=',dummy(i,j)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,
	1    DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = TMP'
	   STOP
	ELSE
C       PRINT*, 'wrote TMP'
	END IF

        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(2)
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF (DUMMY(I,J).LE.0.0000001) THEN
		    PRINT*,'CORRECTING, VALUE OUT OF RANGE, 
	1		 SPFH(',I,',',J,')=',DUMMY(I,J)
		    vap=(1.0/100.0) * 6.11 * 
	1		 10**((7.5*(GRID(I,J)%FORCING(1)-273.15))/
	2		 (237.7+(GRID(I,J)%FORCING(1)-273.15)))
		    
		    rho=0.62197 * ( vap / ((GRID(I,J)%FORCING(7)/100.0) 
	1		 + vap*(.062197-1) ) )
		    dummy(i,j)=rho
		    print *,'corrected spfh value=',dummy(i,j)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM
	
C       MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
C	print*, KPDS
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,
	1    DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SPFH'
	   STOP
	ELSE
C       PRINT*, 'wrote SPFH'
	END IF

	
        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(5)
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF ((DUMMY(I,J).GT.75.0).OR.(DUMMY(I,J).LT.-75.0)) THEN
		    PRINT*,'VALUE OUT OF RANGE, UGRD(',I,',',J,')=',DUMMY(I,J)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM
	
C       MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,
	1    DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = UGRD'
	   STOP
	ELSE
C       PRINT*, 'wrote UGRD'
	END IF
	
        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(6)
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF ((DUMMY(I,J).GT.75.0).OR.(DUMMY(I,J).LT.-75.0)) THEN
		    PRINT*,'VALUE OUT OF RANGE, VGRD(',I,',',J,')=',DUMMY(I,J)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM
	
C     MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,
	1    DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = VGRD'
	   STOP
	ELSE
C       PRINT*, 'wrote VGRD'
	END IF
	

        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(7)
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF ((DUMMY(I,J).GT.110000.0).OR.
	1	      (DUMMY(I,J).LT.5000.0)) THEN
		    PRINT*,'VALUE OUT OF RANGE, PRES(',I,',',J,')=',DUMMY(I,J)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM
	
C       MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,
	1    DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = PRES'
	   STOP
	ELSE
C       PRINT*, 'wrote PRES'
	END IF

        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(3)
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF ((DUMMY(I,J).GT.1360.0).OR.(DUMMY(I,J).LT.0.0)) THEN
		    PRINT*,'VALUE OUT OF RANGE,DSWRF(',I,',',J,')=',DUMMY(I,J)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
	IF (LDAS%AGRMETSW>0) THEN
	   KPDS(2)=087		!ID FOR AGRMET PRODUCT
	ELSE
	   KPDS(2)=origKPDS2    !Maintain model ID
	END IF
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C       MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DSWRF'
	   STOP
	ELSE
C       PRINT*, 'wrote DSWRF'
	END IF
	
        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(4)
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF ((DUMMY(I,J).GT.800.0).OR.(DUMMY(I,J).LT.50.0)) THEN
		    PRINT*,'VALUE OUT OF RANGE,DLWRF(',I,',',J,')=',DUMMY(I,J)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
	IF (LDAS%AGRMETLW>0) THEN
	   KPDS(2)=026		!ID FOR RTNEPH PRODUCT
	ELSE
	   KPDS(2)=origKPDS2    !Maintain model ID
	END IF
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM
	
C       MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,
	1    DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DLWRF'
	   STOP
	ELSE
C       PRINT*, 'wrote DLWRF'
	END IF
	
	
        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(8)*(LDAS%WRITEINTF)*3600.0
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF (DUMMY(I,J).LT.0.0) THEN
		    PRINT*,'VALUE OUT OF RANGE,APCP(',I,',',J,')=',DUMMY(I,J)
		    DUMMY(I,J)=0.0
		    PRINT *,'set APCP to',dummy(i,j)
		 ENDIF
		 IF (DUMMY(I,J).GT.200.0) THEN
		    PRINT*,'VALUE OUT OF RANGE,APCP(',I,',',J,')=',DUMMY(I,J)
		    DUMMY(I,J)=200.0
		    PRINT *,'set APCP to',dummy(i,j)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
        IF (LDAS%GPCPSRC(1) > 0 .OR. LDAS%GPCPSRC(2) > 0
     &      .OR. LDAS%GPCPSRC(3) > 0) THEN
	   KPDS(2)=009		!ID FOR GPCP PRODUCT?? no process # known
	ELSE
	   KPDS(2)=origKPDS2    !Maintain model ID
	END IF
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM
	
C       MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,
	1    DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = APCP'
	   STOP
	ELSE
C       PRINT*, 'wrote APCP'
	END IF
	
	
        DO J=1,LDAS%NR
	   DO I=1,LDAS%NC
	      DUMMY(I,J)=GRID(I,J)%FORCING(9)*(LDAS%WRITEINTF)*3600.0
	      IF ((GRID(I,J)%FMASK.GE.1.0).AND.
	1	   (DUMMY(I,J).NE.LDAS%UDEF)) THEN
		 IF (DUMMY(I,J).LT.0.0) THEN
		    PRINT*,'VALUE OUT OF RANGE,ACPCP(',I,',',J,')=',DUMMY(I,J)
		    DUMMY(I,J)=0.0
		    PRINT *,'set ACPCP to',dummy(i,j)
		 ENDIF
		 IF (DUMMY(I,J).GT.200.0) THEN
		    PRINT*,'VALUE OUT OF RANGE,ACPCP(',I,',',J,')=',DUMMY(I,J)
		    DUMMY(I,J)=200.0 
		    PRINT *,'set ACPCP to',dummy(i,j)
		 ENDIF
	      ENDIF
	   ENDDO
        ENDDO
        IF (LDAS%GPCPSRC(1) > 0 .OR. LDAS%GPCPSRC(2) > 0
     &      .OR. LDAS%GPCPSRC(3) > 0) THEN
	   KPDS(2)=009		!ID FOR GPCP PRODUCT?? no process # known
	ELSE
	   KPDS(2)=origKPDS2    !Maintain model ID
	END IF
	READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
	1    KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM
	
	
C       MAKE TIME DEPENDENT PDS PARAMETERS
	CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
	CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,
	1    DUMMY,IRET)
	IF (IRET .NE. 0) THEN
	   PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ACPCP'
	ELSE
C       PRINT*, 'wrote ACPCP'
C       STOP
	END IF
	
	
	CALL BACLOSE (LUGB, JRET)
	CLOSE (30)
	END
