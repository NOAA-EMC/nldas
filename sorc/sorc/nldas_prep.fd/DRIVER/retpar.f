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
! retrad.f:
!
! DESCRIPTION:
!  Opens, reads, interpolates and overlays PAR radiation forcing.
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data
!
!
! REVISION HISTORY:
!  20  Dec 2000: Brian Cosgrove; Initial code based on retrad.f
!  06  Feb 2002: Brian Cosgrove; Modified to use with new expanded radiation domain
!  15  May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new 
!                GRIB libraries
!=========================================================================
        

	Subroutine RETPAR(ORDER,LDAS,GRID,NAME1,FERROR,FLAG,yr,mo)
	
	USE ldas_module      	! LDAS non-model-specIFic 1-D variables
      	USE grid_module      	! LDAS non-model-specIFic grid variables
      	IMPLICIT NONE
      	TYPE (ldasdec) LDAS
      	TYPE (griddec) GRID(LDAS%NC,LDAS%NR)

!=== Local Variables =====================================================
        CHARACTER*80 NAME1
        CHARACTER*1 dirnom(42),infile(24),dirnom1(42),infile1(24)
        CHARACTER*24 ndummy
	INTEGER FLAG            !INPUT DATA, 1=PAR 
        INTEGER ORDER  		!RETRIEVE PARDATA FOR TIME 1 or TIME 2
        INTEGER DATAFLAG   	!1=edas/eta, 0=ncep/hourly
        INTEGER FERROR		!0=No radiation data found
                             	!1=Found PAR
                             	!2=Found undefined data
        INTEGER FTYPE      	!file type, EDAS,ETA3hr, ETA6hr
        INTEGER NR,NC   	!#of Rows & Columns in LDAS Grid
        INTEGER XSIZE,YSIZE,SUCCESSFLAG,C,R
        INTEGER SIXOCTET,NINEOCTET,TENOCTET,TWEN1OCT
        INTEGER VARADJ
        INTEGER ERRORFLAG
        PARAMETER (XSIZE=111)
        PARAMETER (YSIZE=51)
        REAL VARFIELD(LDAS%NC,LDAS%NR),EMPTYVF(LDAS%NC,LDAS%NR)
	REAL OUTDATA(LDAS%NC,LDAS%NR)
        INTEGER NLDAS, NX, NY, IP, IPOPT(20),
     +          IRET,I,X,Y,CR,
     +          KGDSPAR(200),NPAR,IPPAR,IPPAROPT(20),KMPAR,
     +          IBIPAR,nparnew,npareight,
     +          IBO, NO,yr,mo,
     +          KGDSLDAS(200),NHX,NHY,COUNT
C
        PARAMETER (NLDAS = 103936, NX = 464, NY = 224,
     +             NHX=111,NHY=51,NPAR=5661,nparnew=7381,
     &  npareight=88641)
        REAL LAT(NPAR),LON(NPAR),latnew(nparnew),lonnew(nparnew)
	real lateight(npareight),loneight(npareight)
C
        REAL PARLDAS1D(NLDAS),RLAT(NLDAS),RLON(NLDAS)
        REAL PARRAD1D(NPAR)
        REAL PARRAD1DNEW(NPARNEW)
	real PARRAD1DEIGHT(NPAREIGHT)

        LOGICAL*1 LBPAR(NPAR), LO(NLDAS),lbparnew(nparnew)
	LOGICAL*1 LBpareight(npareight)
!=== End Variable Definition =============================================
!=== Set GRIB octet values to those suitable for radiation
!=== Initialize variables
	SIXOCTET=154
	NINEOCTET=204
	TENOCTET=1
	TWEN1OCT=0
	VARADJ=0
	FERROR=0

!=== If using PAR data, open and read in PAR files
!=== If error reading file, goto 100, FERROR=0
	IF (FLAG.EQ.1) THEN
         OPEN (UNIT=70,FILE=NAME1,STATUS='OLD',ERR=100)
        IF ( (YR.GE.2002).OR.((YR.EQ.2001).AND.(MO.GE.7))) THEN
         DO CR=1,NPARnew
           READ (70,'(3F8.1)') LATnew(CR),LONnew(CR),
     &     PARRAD1Dnew(CR)
         ENDDO
	ENDIF
        IF ( (YR.eq.2001).AND.((MO.GE.1).AND.(MO.LE.6))) THEN
         DO CR=1,NPAR
           READ (70,'(3F8.1)') LAT(CR),LON(CR),
     &     PARRAD1D(CR)
         ENDDO
        ENDIF
        IF (YR.le.2000) THEN
         DO CR=1,NPAReight
           READ (70,'(3F11.3)') LATeight(CR),LONeight(CR),
     &     PARRAD1Deight(CR)
         ENDDO
        ENDIF



	 CLOSE (70)
        
        DO I=1,200
         KGDSLDAS(I)=0
         KGDSPAR(I)=0
        ENDDO

!=== Interpolate PAR data from .5 to 1/8th degree

         KGDSLDAS(1)=0
         KGDSLDAS(2)=464
         KGDSLDAS(3)=224
         KGDSLDAS(4)=25063
         KGDSLDAS(5)=-124938
         KGDSLDAS(6)=128
         KGDSLDAS(7)=52938
         KGDSLDAS(8)=-67063
         KGDSLDAS(9)=125
         KGDSLDAS(10)=125
         KGDSLDAS(11)=64
         KGDSLDAS(12)=0
         KGDSLDAS(13)=0
         KGDSLDAS(14)=0
         KGDSLDAS(15)=0
         KGDSLDAS(16)=0
         KGDSLDAS(17)=0
         KGDSLDAS(18)=0
         KGDSLDAS(19)=0
         KGDSLDAS(20)=255
         KGDSLDAS(21)=0
         KGDSLDAS(22)=0

        IF  (YR.LE.2000) THEN
         KGDSPAR(1)=0
         KGDSPAR(2)=441
         KGDSPAR(3)=201
         KGDSPAR(4)=25000
         KGDSPAR(5)=-125000
         KGDSPAR(6)=128
         KGDSPAR(7)=50000
         KGDSPAR(8)=-70000
         KGDSPAR(9)=125
         KGDSPAR(10)=125
         KGDSPAR(11)=64
         KGDSPAR(12)=0
         KGDSPAR(13)=0
         KGDSPAR(14)=0
         KGDSPAR(15)=0
         KGDSPAR(16)=0
         KGDSPAR(17)=0
         KGDSPAR(18)=0
         KGDSPAR(19)=0
         KGDSPAR(20)=255
         KGDSPAR(21)=0
         KGDSPAR(22)=0
         IPPAR = 0
         IPPAROPT(1) = 2
         IPPAROPT(2) = -1
         KMPAR = 1
         IBIPAR = 1
         DO I=1,NPAReight
           IF (PARRAD1Deight(I).GT.1.0) THEN
            LBPAReight(I)=.TRUE.
           ELSE
            LBPAReight(I)=.FALSE.
           ENDIF
         ENDDO
         DO I=1,NLDAS
           LO(I)=.TRUE.
         ENDDO
         CALL IPOLATES (IPPAR,IPPAROPT,KGDSPAR,KGDSLDAS,NPAReight,
     &                   NLDAS,KMPAR,IBIPAR,
     +                   LBPAReight,PARRAD1Deight,NO,RLAT,RLON,IBO,LO,
     &                   PARLDAS1D,IRET)

        ENDIF  !endif the data from 2000 or before

        IF ( (YR.eq.2001).and.((MO.ge.1).AND.(MO.LE.6))) THEN
         KGDSPAR(1)=0
         KGDSPAR(2)=111
         KGDSPAR(3)=51
         KGDSPAR(4)=25000
         KGDSPAR(5)=-125000
         KGDSPAR(6)=128
         KGDSPAR(7)=50000
         KGDSPAR(8)=-70000
         KGDSPAR(9)=500
         KGDSPAR(10)=500
         KGDSPAR(11)=64
         KGDSPAR(12)=0
         KGDSPAR(13)=0
         KGDSPAR(14)=0
         KGDSPAR(15)=0
         KGDSPAR(16)=0
         KGDSPAR(17)=0
         KGDSPAR(18)=0
         KGDSPAR(19)=0
         KGDSPAR(20)=255
         KGDSPAR(21)=0
         KGDSPAR(22)=0
         IPPAR = 0 
         IPPAROPT(1) = 4
         IPPAROPT(2) = -1
         KMPAR = 1
         IBIPAR = 1
         DO I=1,NPAR
           IF (PARRAD1D(I).GT.1.0) THEN
            LBPAR(I)=.TRUE.
           ELSE
            LBPAR(I)=.FALSE.
           ENDIF
         ENDDO
         DO I=1,NLDAS
           LO(I)=.TRUE.
         ENDDO
         CALL IPOLATES (IPPAR,IPPAROPT,KGDSPAR,KGDSLDAS,NPAR,
     &                   NLDAS,KMPAR,IBIPAR,
     +                   LBPAR,PARRAD1D,NO,RLAT,RLON,IBO,LO,
     &                   PARLDAS1D,IRET)

	ENDIF  !endif the data from January to June 2001

        IF ( (YR.GE.2002).OR.((YR.EQ.2001).AND.(MO.GE.7))) THEN
         KGDSPAR(1)=0
         KGDSPAR(2)=121
         KGDSPAR(3)=61
         KGDSPAR(4)=24000
         KGDSPAR(5)=-126000
         KGDSPAR(6)=128
         KGDSPAR(7)=54000
         KGDSPAR(8)=-66000
         KGDSPAR(9)=500
         KGDSPAR(10)=500
         KGDSPAR(11)=64
         KGDSPAR(12)=0
         KGDSPAR(13)=0
         KGDSPAR(14)=0
         KGDSPAR(15)=0
         KGDSPAR(16)=0
         KGDSPAR(17)=0
         KGDSPAR(18)=0
         KGDSPAR(19)=0
         KGDSPAR(20)=255
         KGDSPAR(21)=0
         KGDSPAR(22)=0
         IPPAR = 0
         IPPAROPT(1) = 4
         IPPAROPT(2) = -1
         KMPAR = 1
         IBIPAR = 1
         DO I=1,NPARNEW
           IF (PARRAD1DNEW(I).GT.1.0) THEN
            LBPARNEW(I)=.TRUE.
           ELSE
            LBPARNEW(I)=.FALSE.
           ENDIF
         ENDDO
         DO I=1,NLDAS
           LO(I)=.TRUE.
         ENDDO
         CALL IPOLATES (IPPAR,IPPAROPT,KGDSPAR,KGDSLDAS,NPARNEW,
     &                   NLDAS,KMPAR,IBIPAR,
     +                   LBPARNEW,PARRAD1DNEW,NO,RLAT,RLON,IBO,LO,
     &                   PARLDAS1D,IRET)

        ENDIF  !endif the data from July 2001 or after

C       Transfer 1D data to 2D array
         COUNT = 0
         DO Y = 1, NY
           DO X = 1, NX
             OUTDATA(X,Y) = PARLDAS1D(X+COUNT)
           END DO
           COUNT = COUNT + NX
         END DO

	 FERROR=2

!=== If PAR data is not undefined, set FERROR to 1 
	 DO R=1,LDAS%NR
          DO C=1,LDAS%NC
	   IF(OUTDATA(C,R).GE.0.0) FERROR=1
	  ENDDO
	 ENDDO 

  100  	 CONTINUE

	  
!=== If PAR data is defined, transfer it into the proper
!=== output array.  If it's undefined, fill the output array with
!=== undefined values

	 IF (FERROR.EQ.1) THEN
	  DO R=1,LDAS%NR
	   DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
	     GRID(C,R)%PARDATA1=OUTDATA(C,R)
	    ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%PARDATA2=OUTDATA(C,R)
	    ENDIF
	   ENDDO
	  ENDDO
	 ELSE
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
             GRID(C,R)%PARDATA1=LDAS%UDEF
            ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%PARDATA2=LDAS%UDEF
            ENDIF
           ENDDO
          ENDDO
	 ENDIF
	ENDIF  !END the flag.eq.1 loop


	RETURN
	END

