	subroutine gribout (ldas,grid,tile,mos,fbase,fyrmodir)

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
C  01-01-19 Cosgrove Added CLOSE statement
C  01-04-05 Cosgrove Changed definition of KPDS 17 (number of timesteps
C           in accumulation or average)...made it variable instead of 
C           hardwired 
C  01-09-04 Cosgrove Updated--Changed output directory structure, enabled
C           statistical output, changed real to real*8 for tick call.  
C           Volumetric soil moisture output is commented out right now,
C           but could be enabled to output 5cm,25cm,60cm,75cm,tot-col,
C           layer1,layer2 and layer3 vol. soil mst.
C  02-02-05 Jambor Removed misleading code regarding KPDS(13) values
C           that had been commented out.
C  02-05-15 Jambor; Changed LOGICAL to LOGICAL*1,matches new GRIB libraries
C  02-06-04 Jambor Reverted to using old version of KPDS table for
C           global domain runs for consistency with existing output.
C  02-07-12 Cosgrove; Made old version of KPDS table the default used
C           for all domains >0.  This should be changed to >1 when 
C           performing -official group- NLDAS runs so that NLDAS 
C           official KPDS table is used.
C***********************************************************************

! Declare modules and data structures
      use ldas_module      ! ldas non-model-specific 1-d variables
      use tile_module      ! ldas non-model-specific tile variables
      use grid_module      ! ldas non-model-specific grid variables
      use mos_module
      implicit none
      type (ldasdec) ldas
      type (tiledec) tile(ldas%nch)
      type (mosdec)  mos(ldas%nch)
      type (griddec) grid(ldas%nc,ldas%nr)

      INTEGER c,r,t,m,cc,jj,counter
      INTEGER SS1,TS,MN1,HR1,DA1,MO1,YR1,TS1,DOY1
      INTEGER NSOLD
      PARAMETER (NSOLD=8)
      CHARACTER*8 TODAY, YESTERDAY
      CHARACTER*1 TOD(8),YES(8)
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13),FTIMED(4)
      CHARACTER*1  FYRMODIR(18),FSUBFT(80),EXPCODE(3),FTIMEC(8)
      CHARACTER*1  FSUBFG(80),EXPARRAY(80),FTIMEB(10),FSUBGB(11)
      CHARACTER*1 GRIBF(256),GRIBF2(256)

      INTEGER LUGB, I, J, K, N, IRET, KPDS(25), KGDS(22),
     &        NLDAS, LENGDS, IG, JRET, NFIELDS

      LOGICAL*1  LDASMASK(LDAS%NC,LDAS%NR)

      CHARACTER CENT*2, CHOUR*2, GDS(400)
      CHARACTER*256 GRIBFILE

      REAL    VMEAN,VSTDEV,VMIN,VMAX
      REAL    NSWRS(LDAS%NC,LDAS%NR) 
      REAL    NLWRS(LDAS%NC,LDAS%NR)
      REAL    LHTFL(LDAS%NC,LDAS%NR) 
      REAL    SHTFL(LDAS%NC,LDAS%NR) 
      REAL    GFLUX(LDAS%NC,LDAS%NR) 
      REAL    DSWRF(LDAS%NC,LDAS%NR)
      REAL    DLWRF(LDAS%NC,LDAS%NR)
      REAL    ARAIN(LDAS%NC,LDAS%NR)
      REAL    EVP(LDAS%NC,LDAS%NR)
      REAL    SSRUN(LDAS%NC,LDAS%NR)
      REAL    BGRUN(LDAS%NC,LDAS%NR) 
      REAL    AVSFT(LDAS%NC,LDAS%NR)
      REAL    ALBDO(LDAS%NC,LDAS%NR) 
      REAL    WEASD(LDAS%NC,LDAS%NR)
      REAL    CWAT(LDAS%NC,LDAS%NR)
      REAL    SOILT(nsold,LDAS%NC,LDAS%NR)
      REAL    SOILM(nsold,LDAS%NC,LDAS%NR)
      REAL    SOILV(nsold,LDAS%NC,LDAS%NR)
      REAL    MSTAV(nsold,LDAS%NC,LDAS%NR)
      REAL    PEVPR(LDAS%NC,LDAS%NR)
      REAL    VEG(LDAS%NC,LDAS%NR)
      REAL    LAI(LDAS%NC,LDAS%NR)
      REAL    SNOD(LDAS%NC,LDAS%NR)
      REAL    SNOHF(LDAS%NC,LDAS%NR)
      REAL    ASNOW(LDAS%NC,LDAS%NR)
      REAL    SNOM(LDAS%NC,LDAS%NR)
      REAL    TRANS(LDAS%NC,LDAS%NR)
      REAL    EVBS(LDAS%NC,LDAS%NR)
      REAL    EVCW(LDAS%NC,LDAS%NR)
      REAL    SBSNO(LDAS%NC,LDAS%NR)
      REAL    ACOND(LDAS%NC,LDAS%NR)
      REAL    CCOND(LDAS%NC,LDAS%NR)
      REAL    SNOC(LDAS%NC,LDAS%NR)
      REAL    TMP(LDAS%NC,LDAS%NR)
      REAL    SPFH(LDAS%NC,LDAS%NR)
      REAL    UGRD(LDAS%NC,LDAS%NR)
      REAL    VGRD(LDAS%NC,LDAS%NR)
      REAL    PRES(LDAS%NC,LDAS%NR)
      REAL    ACPCP(LDAS%NC,LDAS%NR)
      REAL DUMMYGMT,DUMMY(LDAS%NC,LDAS%NR)
      REAL*8 DUMMYTIME

!	Initialize variables to zero
	DO j=1,LDAS%NR
	DO i=1,LDAS%NC
         NSWRS(I,J)=0.0
         NLWRS(I,J)=0.0
         LHTFL(I,J)=0.0
         SHTFL(I,J)=0.0
         GFLUX(I,J)=0.0
         DSWRF(I,J)=0.0
         DLWRF(I,J)=0.0
         ARAIN(I,J)=0.0
         EVP(I,J)=0.0
         SSRUN(I,J)=0.0
         BGRUN(I,J)=0.0
         AVSFT(I,J)=0.0
         ALBDO(I,J)=0.0
         WEASD(I,J)=0.0
         CWAT(I,J)=0.0
	 DO k=1,nsold
          SOILV(k,i,j)=0.0
          SOILT(k,i,j)=0.0
          SOILM(k,i,j)=0.0
          MSTAV(k,i,j)=0.0
	 ENDDO
         PEVPR(I,J)=0.0
         VEG(I,J)=0.0
         LAI(I,J)=0.0
         SNOD(I,J)=0.0
         SNOHF(I,J)=0.0
         ASNOW(I,J)=0.0
         SNOM(I,J)=0.0
         DUMMY(I,J)=0.0
         EVCW(I,J)=0.0
         TRANS(I,J)=0.0
         EVBS(I,J)=0.0
         SBSNO(I,J)=0.0
         ACOND(I,J)=0.0
         CCOND(I,J)=0.0
         SNOC(I,J)=0.0
         TMP(I,J)=0.0
         SPFH(I,J)=0.0
         UGRD(I,J)=0.0
         VGRD(I,J)=0.0
         PRES(I,J)=0.0
         ACPCP(I,J)=0.0
	ENDDO
	ENDDO
!	Construct GRIB logical type land/sea mask from LDAS integer land/sea mask
      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
          IF(GRID(C,R)%IMASK.EQ.1) THEN
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
        DO I=1,25
         KPDS(I)=0
        ENDDO
        DO I=1,22
         KGDS(I)=0
        ENDDO
C     NON-CHANGING PDS ELEMENTS
      KPDS(1)=222               !ID FOR GSFC PRODUCTS
      KPDS(2)=222               !ID FOR Mosaic MODEL (change value for other models)
      KPDS(4)=192               !BMS FLAG... DOn't worry about this.
      KPDS(12)=0                !assume output time minute always = 0
      KPDS(13)=1                !Forecast Time Unit (Hours)
      KPDS(17)=INT((LDAS%WRITEINTM*3600.0)/LDAS%TS) !number of time steps in
                                                    !averaged/accum variables
      KPDS(18)=0                !GRIB version -- left as 0 in NCEP products
      KPDS(19)=1                !version number of KPDS.tbl for LDAS.  
      KPDS(20)=0                !none missing from averages/accumulations (always4)
      KPDS(23)=222              !GSFC ID#
      KPDS(24)=0                !DOes not apply to LDAS output
      KPDS(25)=0                !not used

C          IF(GRID(C,R)%IMASK.EQ.0)Gtmp(C,R)=LDAS%UDEF  !Set Water to undefined
       DO T=1,LDAS%NCH
!	Time Averaged
        NSWRS(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   NSWRS(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(18)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        NLWRS(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   NLWRS(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(19)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        LHTFL(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   LHTFL(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(21)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        SHTFL(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SHTFL(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(20)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        GFLUX(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   GFLUX(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(22)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        SNOHF(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SNOHF(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(45)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        DSWRF(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   DSWRF(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(37)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        DLWRF(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   DLWRF(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(38)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
!	Accumulated
        ASNOW(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   ASNOW(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(28)*TILE(T)%FGRD
        ARAIN(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   ARAIN(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(27)*TILE(T)%FGRD
        EVP(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   EVP(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(23)*TILE(T)%FGRD
        SSRUN(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SSRUN(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(33)*TILE(T)%FGRD
        BGRUN(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   BGRUN(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(46)*TILE(T)%FGRD
        SNOM(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SNOM(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(29)*TILE(T)%FGRD
!	Instantaneous
        AVSFT(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   AVSFT(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(1)*TILE(T)%FGRD
        ALBDO(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   ALBDO(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(39)*TILE(T)%FGRD
        WEASD(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   WEASD(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(5)*TILE(T)%FGRD
        CWAT(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   CWAT(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(4)*TILE(T)%FGRD
        SOILT(1,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILT(1,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(2)*TILE(T)%FGRD
        SOILV(1,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILV(1,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(61)*TILE(T)%FGRD
        SOILV(2,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILV(2,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(62)*TILE(T)%FGRD
        SOILV(3,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILV(3,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(63)*TILE(T)%FGRD
        SOILV(4,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILV(4,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(64)*TILE(T)%FGRD
        SOILV(5,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILV(5,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(60)*TILE(T)%FGRD
        SOILV(6,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILV(6,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(14)*TILE(T)%FGRD
        SOILV(7,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILV(7,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(15)*TILE(T)%FGRD
        SOILV(8,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILV(8,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(16)*TILE(T)%FGRD
        SOILM(1,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILM(1,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(13)*TILE(T)%FGRD
        SOILM(2,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILM(2,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(52)*TILE(T)%FGRD
        SOILM(3,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILM(3,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(51)*TILE(T)%FGRD
        SOILM(4,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILM(4,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(10)*TILE(T)%FGRD
        SOILM(5,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILM(5,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(11)*TILE(T)%FGRD
        SOILM(6,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SOILM(6,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(12)*TILE(T)%FGRD
        MSTAV(1,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   MSTAV(1,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(17)*TILE(T)%FGRD
        MSTAV(2,TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   MSTAV(2,TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(53)*TILE(T)%FGRD
!	Time Averaged
        EVCW(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   EVCW(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(24)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        TRANS(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   TRANS(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(26)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)

        EVBS(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   EVBS(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(25)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
        SBSNO(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SBSNO(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(40)*
     2    TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
!	Instantaneous

        ACOND(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   ACOND(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(43)*TILE(T)%FGRD
        CCOND(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   CCOND(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(44)*TILE(T)%FGRD
        VEG(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   VEG(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(35)*TILE(T)%FGRD
        LAI(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   LAI(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(34)*TILE(T)%FGRD
        SNOD(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SNOD(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(49)*TILE(T)%FGRD
        SNOC(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SNOC(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(50)*TILE(T)%FGRD
        TMP(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   TMP(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(54)*TILE(T)%FGRD
        SPFH(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   SPFH(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(55)*TILE(T)%FGRD
        UGRD(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   UGRD(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(56)*TILE(T)%FGRD
        VGRD(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   VGRD(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(57)*TILE(T)%FGRD
        PRES(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   PRES(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(58)*TILE(T)%FGRD
! 	Accumulation
        ACPCP(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   ACPCP(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%TOTRET(59)*TILE(T)%FGRD

       ENDDO
C     DEFINE THE GDS FOR THE LDAS GRID (THIS ARRAY NEVER CHANGES FROM
C     FIELD TO FIELD. ONLY HAS TO BE DEFINED ONCE, UNLIKE THE PDS).

      IF (LDAS%DOMAIN==1) THEN
        IG = 236    !FOR LDAS GRID IN SUBROUTINE w3fi71.f in w3lib

        CALL MAKGDS(IG, KGDS, GDS, LENGDS, IRET)
      ELSE
        DO J=1,22
          KGDS(J)=LDAS%LDAS_KGDS(J)
        ENDDO
      ENDIF

C     read past header info in KPDS.tbl

c############################################################
C     IMPORTANT  IMPORTANT  

C     For most runs, use this line to ensure maximum GrADS compatibility
C     with GrADS Z layers 1 to 6 corresponding to total col,rt zone,
C     top1m,layer1,layer2,layer3
C     IF (LDAS%DOMAIN > 0) THEN

C     For OFFICIAL GROUP NLDAS runs, use this line to comply with the 
C     NLDAS GRIB soil layer code information, i.e., 20099 would
C     be used to indicate a 2meter total column variable, with
C     GrADS Z layers 1 to 6 corresponding to totcol, level3, top1m,
C     root zone, level 2 , level 1
C     IF (LDAS%DOMAIN > 1) THEN

      IF (LDAS%DOMAIN > 0) THEN
        OPEN (UNIT = 30, FILE = './SRC/KPDS_completemosaic.old.tbl')
      ELSE
        OPEN (UNIT = 30, FILE = './SRC/KPDS_completemosaic.tbl')
      ENDIF
      DO K = 1, 42
         READ(30,*)
      END DO

C     MAKE OUTPUT FILENAME
 92    FORMAT(80A1)
 93    FORMAT(A256)
 91    FORMAT(256A1)
 94    FORMAT(A80)
 95    format(80A,A4,3A1,A5,4A1,A1,8A1,A1,10A1,11A1)
101    FORMAT(A11)
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


        WRITE(95,101,REC=1)'.mosaic.grb'
        READ(95,92,REC=1) (FSUBGB(I),I=1,11)
        C=0
       DO I=1,40
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO



C      IF ( (HR1.LE. 23).AND.(HR1.GE.13) ) THEN
         WRITE (CHOUR, '(i2.2)') HR1
        WRITE(95,103,REC=1)LDAS%EXPCODE
c	WRITE(95,102,REC=1)LDAS%EXPCODE
        READ(95,92,REC=1)EXPARRAY
	COUNTER=1
        DO I=1,40
         IF(EXPARRAY(I).NE.(' ')) THEN
          EXPCODE(COUNTER)=EXPARRAY(I)
          COUNTER=COUNTER+1
	 ENDIF
        ENDDO

c        WRITE(95,95,REC=1)(FBASE(I),I=1,C),FYRMODIR,'/',
c     &  'LDAS.E',EXPCODE,'.',
c     &   TODAY,CHOUR,'.MOSgrib'
        WRITE(95,95,REC=1)(FBASE(I),I=1,C),'/EXP',EXPCODE,'/MOS/',
     &  FTIMED,'/',FTIMEC,'/',
     &   FTIMEB,(FSUBGB(I),I=1,11)
	READ(95,93,REC=1) GRIBFILE
c        print *,'gribfile=',gribfile,'end',LDAS%EXPCODE,'end'
        CLOSE (95)
C      ELSE 
C         WRITE (CHOUR, '(i2.2)') HR1
C        WRITE(95,95,REC=1)(FBASE(I),I=1,C),FYRMODIR,'/',
C     &   YESTERDAY,CHOUR,'_MOSAIC.LDASGRIB'
C        READ(95,93,REC=1) GRIBFILE
C      END IF

C     OPEN GRIB FILE

      LUGB = 40
      CALL BAOPEN (LUGB, GRIBFILE, IRET)

C     WRITE EACH FIELD (IN GRIB) TO HR1LY OUTPUT FIELD
 15   FORMAT (29x, 7I6) 

      NLDAS = LDAS%NC * LDAS%NR

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,NSWRS,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = NSWRS'
         STOP
      ELSE
C         PRINT*, 'wrote NSWRS'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,NLWRS,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = NLWRS'
         STOP
      ELSE
C         PRINT*, 'wrote NLWRS'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,LHTFL,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = LHTFL'
         STOP
      ELSE
C         PRINT*, 'wrote LHTFL'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SHTFL,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SHTFL'
         STOP
      ELSE
C         PRINT*, 'wrote SHTFL'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,GFLUX,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = GFLUX'
         STOP
      ELSE
C         PRINT*, 'wrote GFLUX'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOHF,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOHF'
         STOP
      ELSE
C         PRINT*, 'wrote SNOHF'
      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DSWRF,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DSWRF'
         STOP
      ELSE
C         PRINT*, 'wrote DSWRF'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DLWRF,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DLWRF'
         STOP
      ELSE
C         PRINT*, 'wrote DLWRF'
      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ASNOW,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ASNOW'
         STOP
      ELSE
C         PRINT*, 'wrote ASNOW'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ARAIN,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ARAIN'
         STOP
      ELSE
C         PRINT*, 'wrote ARAIN'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,EVP,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = EVP'
         STOP
      ELSE
C         PRINT*, 'wrote EVP'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SSRUN,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SSRUN'
         STOP
      ELSE
C         PRINT*, 'wrote SSRUN'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,BGRUN,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = BGRUN'
         STOP
      ELSE
C         PRINT*, 'wrote BGRUN'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOM,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOM'
         STOP
      ELSE
C         PRINT*, 'wrote SNOM'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,AVSFT,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = AVSFT'
         STOP
      ELSE
C         PRINT*, 'wrote AVSFT'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ALBDO,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ALBDO'
         STOP
      ELSE
C         PRINT*, 'wrote ALBDO'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,WEASD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = WEASD'
         STOP
      ELSE
C         PRINT*, 'wrote WEASD'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,CWAT,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = CWAT'
         STOP
      ELSE
C         PRINT*, 'wrote CWAT'
      END IF




      DO N = 1, 1
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILT(N,I,J)
            END DO
         END DO
         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +        KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C        MAKE TIME DEPENDENT PDS PARAMETERS
         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = SOILT',N
            STOP
         ELSE
C         PRINT*, 'wrote SOILT ', N
         END IF
      END DO

c      DO N = 1, 8
c         DO J = 1,LDAS%NR
c            DO I = 1,LDAS%NC
c               DUMMY(I,J) = SOILV(N,I,J)
c            END DO
c         END DO
c         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
c     +        KPDS(15), KPDS(16), KPDS(22)
c        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM
c
C        MAKE TIME DEPENDENT PDS PARAMETERS
c         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
c         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
c         IF (IRET .NE. 0) THEN
c            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = SOILV',N
c            STOP
c         ELSE
C         PRINT*, 'wrote SOILV ',N
c         END IF
c      END DO




      DO N = 1, 6
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILM(N,I,J)
            END DO
         END DO
         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +        KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C        MAKE TIME DEPENDENT PDS PARAMETERS
         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = SOILM',N
            STOP
         ELSE
C         PRINT*, 'wrote SOILM ',N
         END IF
      END DO

      DO N = 1, 2
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = MSTAV(N,I,J)*100.0
            END DO
         END DO
         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +        KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C        MAKE TIME DEPENDENT PDS PARAMETERS
         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = MSTAV',N
            STOP
         ELSE
C         PRINT*, 'wrote MSTAV ',N
         END IF
      END DO


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,EVCW,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = EVCW'
         STOP
      ELSE
C         PRINT*, 'wrote EVCW'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,TRANS,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = TRANS'
         STOP
      ELSE
C         PRINT*, 'wrote TRANS'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,EVBS,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = EVBS'
         STOP
      ELSE
C         PRINT*, 'wrote EVBS'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SBSNO,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SBSNO'
         STOP
      ELSE
C         PRINT*, 'wrote SBSNO'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ACOND,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ACOND'
         STOP
      ELSE
C         PRINT*, 'wrote ACOND'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,CCOND,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = CCOND'
         STOP
      ELSE
C         PRINT*, 'wrote CCOND'
      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,VEG,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = VEG'
         STOP
      ELSE
C         PRINT*, 'wrote VEG'
      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,LAI,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = LAI'
         STOP
      ELSE
C         PRINT*, 'wrote LAI'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOD'
         STOP
      ELSE
C         PRINT*, 'wrote SNOD'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOC,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOC'
         STOP
      ELSE
C         PRINT*, 'wrote SNOC'
      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,TMP,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = TMP'
         STOP
      ELSE
C         PRINT*, 'wrote TMP'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SPFH,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SPFH'
         STOP
      ELSE
C         PRINT*, 'wrote SPFH'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,UGRD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = UGRD'
         STOP
      ELSE
C         PRINT*, 'wrote UGRD'
      END IF
      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,VGRD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = VGRD'
         STOP
      ELSE
C         PRINT*, 'wrote VGRD'
      END IF
      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,PRES,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = PRES'
         STOP
      ELSE
C         PRINT*, 'wrote PRES'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%WRITEINTM


C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ACPCP,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ACPCP'
      ELSE
C         PRINT*, 'wrote ACPCP'
C         STOP
      END IF

      CALL BACLOSE (LUGB, JRET)
      CLOSE (30)
!=== Write statistical output
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,NSWRS,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'NSWRS (W/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,NLWRS,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'NLWRS (W/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,LHTFL,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'LHTFL (W/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,SHTFL,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SHTFL (W/m2)',VMEAN,VSTDEV,VMIN,VMAX


       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,GFLUX,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'GFLUX (W/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,SNOHF,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SNOHF (W/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DSWRF,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'DSWRF (W/m2)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DLWRF,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'DLWRF (W/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,ASNOW,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'ASNOW (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,ARAIN,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'ARAIN (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,EVP,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'EVP (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,SSRUN,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SSRUN (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,BGRUN,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'BGRUN (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,SNOM,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SNOM (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,AVSFT,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'AVSFT (K)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,ALBDO,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'ALBDO (fraction)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,WEASD,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'WEASD (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,CWAT,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'CWAT (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,CWAT,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'CWAT (Kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILT(1,I,J)
            END DO
         END DO
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SOILT 1 (K)',VMEAN,VSTDEV,VMIN,VMAX


          DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILV(1,I,J)
            END DO
         END DO

C       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
C     &  VMEAN,VSTDEV,VMIN,VMAX)
C       WRITE(65,999)'SOILV 1',VMEAN,VSTDEV,VMIN,VMAX

         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILM(1,I,J)
            END DO
         END DO

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SOILM 1 (kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILM(2,I,J)
            END DO
         END DO
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SOILM 2 (kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILM(3,I,J)
            END DO
         END DO

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SOILM 3 (kg/m2)',VMEAN,VSTDEV,VMIN,VMAX
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILM(4,I,J)
            END DO
         END DO

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SOILM 4 (kg/m2)',VMEAN,VSTDEV,VMIN,VMAX
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILM(5,I,J)
            END DO
         END DO
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SOILM 5 (kg/m2)',VMEAN,VSTDEV,VMIN,VMAX
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILM(6,I,J)
            END DO
         END DO

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SOILM 6 (kg/m2)',VMEAN,VSTDEV,VMIN,VMAX

         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = MSTAV(1,I,J)*100.0
            END DO
         END DO

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'MSTAV 1 (fraction)',VMEAN,VSTDEV,VMIN,VMAX

         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = MSTAV(2,I,J)*100.0
            END DO
         END DO

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,DUMMY,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'MSTAV 2 (fraction)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,EVCW,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'EVCW (W/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,TRANS,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'TRANS (W/m2)',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,EVBS,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'EVBS (W/m2)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,SBSNO,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SBSNO (W/m2)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,ACOND,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'ACOND (m/s)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,CCOND,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'CCOND (m/s)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,VEG,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'VEG (fraction)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,LAI,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'LAI',VMEAN,VSTDEV,VMIN,VMAX

       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,SNOD,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SNOD (m)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,SNOC,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SNOC (fraction)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,TMP,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'TMP (K)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,SPFH,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'SPFH (kg/kg)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,UGRD,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'UGRD (m/s)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,VGRD,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'VGRD (m/s)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,PRES,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'PRES (mb)',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSGMOS(LDAS,GRID,LDAS%UDEF,ACPCP,LDAS%NC,LDAS%NR,
     &  VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(65,999)'ACPCP (kg/m2)',VMEAN,VSTDEV,VMIN,VMAX
999    FORMAT(1X,A18,4F14.3)

      END
