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
! retradbias.f:
!
! DESCRIPTION:
!  Opens, reads, interpolates and overlays radiation forcing.
!    WORKS ONLY ON NLDAS 1/8th Degree GRID
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data
!
!
! REVISION HISTORY:
!  28  Oct 1999: Brian Cosgrove; Initial code
!  11 Apr 2000: Brian Cosgrove; changed code to use Forcing Mask (With inland
!               water filled in).  Deleted unused variables.
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!                not a hard-wired undefined value of -999.9 and -999.0
!  05 Sep 2001: Brian Cosgrove; changed code so to use GRIB interpolation
!               package and to process BRTTMP data too
!  06  Feb 2002: Brian Cosgrove; Modified to use with new expanded radiation domain
!  15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new 
!                GRIB libraries
!  15 May 2007: Chuck Alonge; Modified to use with GOES/NARR Bias Correction Ratios
!=========================================================================
        
	Subroutine RETRADBC(ORDER,LDAS,GRID,NAME1,FERROR)
	
	USE ldas_module      	! LDAS non-model-specIFic 1-D variables
      	USE grid_module      	! LDAS non-model-specIFic grid variables
      	IMPLICIT NONE
      	TYPE (ldasdec) LDAS
      	TYPE (griddec) GRID(LDAS%NC,LDAS%NR)

!=== Local Variables =====================================================
        CHARACTER*80 NAME1
        CHARACTER*1 dirnom(42),infile(24),dirnom1(42),infile1(24)
        CHARACTER*24 ndummy
	INTEGER FLAG            !INPUT DATA, 1=PINKER 2=NESDIS
        INTEGER DATAFLAG   	!1=edas/eta, 0=ncep/hourly
        INTEGER FERROR		!0=No radiation data found
                             	!1=Found Pinker or NESDIS or BRTTMP
                             	!2=Found undefined data
        INTEGER FTYPE      	!file type, EDAS,ETA3hr, ETA6hr
        INTEGER NR,NC   	!#of Rows & Columns in LDAS Grid
        INTEGER SUCCESSFLAG,C,R
        INTEGER VARADJ,yr,mo
        INTEGER ERRORFLAG, ORDER
        INTEGER NLDAS, NX, NY, IP, IPOPT(20),
     +          IRET,I,X,Y,CR,
     +          KGDSPINK(200),NPINK,IPPINK,IPPINKOPT(20),KMPINK,
     +          IBIPINK,npinknew,npinkeight,
     +          IBO, NO,
     +          KGDSLDAS(200),COUNT
C
        PARAMETER (NLDAS = 103936, NX = 464, NY = 224,
     +             NPINK=5661,npinknew=7381,npinkeight=88641)
        REAL LAT(NPINK),LON(NPINK),latnew(npinknew),lonnew(npinknew)
	real lateight(npinkeight),loneight(npinkeight)

C
        REAL VARFIELD(NX,NY),EMPTYVF(LDAS%NC,LDAS%NR)
	REAL OUTDATA(LDAS%NC,LDAS%NR)

        REAL PINKERLDAS1D(NLDAS),RLAT(NLDAS),RLON(NLDAS)
        REAL PINKERRAD1D(NPINK),pinkerrad1dnew(npinknew)
	REAL PINKERRAD1DEIGHT(NPINKEIGHT)
        LOGICAL*1 LBPINK(NPINK), LBPINKNEW(NPINKNEW),LO(NLDAS)
	LOGICAL*1 LBPINKEIGHT(NPINKEIGHT)

!=== End Variable Definition =============================================

!=== Initialize variables
	VARADJ=0
	FERROR=0


!=== CHECK FOR FORCING GENERATION ON NLDAS GRID (ONLY VALID IN THIS CASE)
        IF(LDAS%NC .EQ. NX .AND. LDAS%NR .EQ. NY) THEN

!=== If using Bias Ratio data, open and read in Bias correction files
!=== If error reading file, goto 100, FERROR=0
!        WRITE(*,*) "IN RETRADBC"
!        WRITE(*,*) "FILENAME: ", NAME1
        OPEN (UNIT=70,FILE=NAME1,STATUS='OLD',FORM="UNFORMATTED")
        READ (70) VARFIELD
	CLOSE (70)
      
!=== If Bias correction data is not undefined, set FERROR to 1 
	 DO R=1,LDAS%NR
          DO C=1,LDAS%NC
	   IF(VARFIELD(C,R).GT.0.0) THEN
             OUTDATA(C,R) = VARFIELD(C,R)
             FERROR=1
           ELSE
             OUTDATA(C,R) = 1.0
             FERROR=1
           ENDIF
	  ENDDO
	 ENDDO 
	  
!=== If bias correction data is defined, transfer it into the proper
!=== output array.  If it's undefined, fill the output array with
!=== undefined values

	 IF (FERROR.EQ.1) THEN
	  DO R=1,LDAS%NR
	   DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
	     GRID(C,R)%RADBCDATA1=OUTDATA(C,R)
	    ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%RADBCDATA2=OUTDATA(C,R)
	    ENDIF
	   ENDDO
	  ENDDO
         ELSE ! NO VALID DATA MAKE RATIO 1
	  DO R=1,LDAS%NR
	   DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
	     GRID(C,R)%RADBCDATA1=1.0
	    ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%RADBCDATA2=1.0
	    ENDIF
	   ENDDO
	  ENDDO
         ENDIF

        ELSE
          WRITE(*,*) "STOP: CAN ONLY USE THIS BIAS" 
          WRITE(*,*) "      CORRECTION ON NLDAS GRID"
          STOP
        ENDIF

	RETURN
	END

