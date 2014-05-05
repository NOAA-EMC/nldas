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
!  Opens, reads, interpolates and overlays radiation forcing.
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
!=========================================================================
        

	Subroutine RETRAD(ORDER,LDAS,GRID,NAME1,FERROR,FLAG,yr,mo)
	
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
        INTEGER ORDER  		!RETRIEVE PINKDATA FOR TIME 1 or TIME 2
        INTEGER DATAFLAG   	!1=edas/eta, 0=ncep/hourly
        INTEGER FERROR		!0=No radiation data found
                             	!1=Found Pinker or NESDIS or BRTTMP
                             	!2=Found undefined data
        INTEGER FTYPE      	!file type, EDAS,ETA3hr, ETA6hr
        INTEGER NR,NC   	!#of Rows & Columns in LDAS Grid
        INTEGER SUCCESSFLAG,C,R
        INTEGER SIXOCTET,NINEOCTET,TENOCTET,TWEN1OCT
        INTEGER VARADJ,yr,mo
        INTEGER ERRORFLAG
        REAL VARFIELD(LDAS%NC,LDAS%NR),EMPTYVF(LDAS%NC,LDAS%NR)
	REAL OUTDATA(LDAS%NC,LDAS%NR)
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
        REAL PINKERLDAS1D(NLDAS),RLAT(NLDAS),RLON(NLDAS)
        REAL PINKERRAD1D(NPINK),pinkerrad1dnew(npinknew)
	REAL PINKERRAD1DEIGHT(NPINKEIGHT)
        LOGICAL*1 LBPINK(NPINK), LBPINKNEW(NPINKNEW),LO(NLDAS)
	LOGICAL*1 LBPINKEIGHT(NPINKEIGHT)

!=== End Variable Definition =============================================

!=== Set GRIB octet values to those suitable for radiation
!=== Initialize variables
	SIXOCTET=154
	NINEOCTET=204
	TENOCTET=1
	TWEN1OCT=0
	VARADJ=0
	FERROR=0

!=== If using Pinker data, open and read in Pinker files
!=== If error reading file, goto 100, FERROR=0
	IF ((FLAG.EQ.1).OR.(FLAG.EQ.4)) THEN
         OPEN (UNIT=70,FILE=NAME1,STATUS='OLD',ERR=100)
	COUNT=1

	IF ( (YR.GE.2002).OR.((YR.EQ.2001).AND.(MO.GE.7))) THEN
        DO CR=1,NPINKNEW
           READ (70,'(3F8.1)') LATNEW(CR),LONNEW(CR),
     &     PINKERRAD1DNEW(CR)
        ENDDO
	ENDIF

        IF ( (YR.EQ.2001).AND.((MO.GE.1).AND.(MO.LE.6))) THEN
	DO CR=1,NPINK
           READ (70,'(3F8.1)') LAT(CR),LON(CR),
     &     PINKERRAD1D(CR)
        ENDDO
	ENDIF

        IF (YR.LE.2000) THEN
        DO CR=1,NPINKeight
           READ (70,'(3F11.3)') LATeight(CR),LONeight(CR),
     &     PINKERRAD1Deight(CR)
        ENDDO
        ENDIF

	 CLOSE (70)
      
 

        DO I=1,200
         KGDSLDAS(I)=0
         KGDSPINK(I)=0
        ENDDO

!=== Interpolate Pinker data from .5 to 1/8th degree
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

        IF (YR.LE.2000) THEN
         KGDSPINK(1)=0
         KGDSPINK(2)=441
         KGDSPINK(3)=201
         KGDSPINK(4)=25000
         KGDSPINK(5)=-125000
         KGDSPINK(6)=128
         KGDSPINK(7)=50000
         KGDSPINK(8)=-70000
         KGDSPINK(9)=125
         KGDSPINK(10)=125
         KGDSPINK(11)=64
         KGDSPINK(12)=0
         KGDSPINK(13)=0
         KGDSPINK(14)=0
         KGDSPINK(15)=0
         KGDSPINK(16)=0
         KGDSPINK(17)=0
         KGDSPINK(18)=0
         KGDSPINK(19)=0
         KGDSPINK(20)=255
         KGDSPINK(21)=0
         KGDSPINK(22)=0
         IPPINK = 0
         IPPINKOPT(1) = 2
         IPPINKOPT(2) = -1
         KMPINK = 1
         IBIPINK = 1
         DO I=1,NPINKeight
           IF (PINKERRAD1Deight(i).GT.1.0) THEN
            LBPINKeight(i)=.TRUE.
           ELSE
            LBPINKeight(i)=.FALSE.
           ENDIF
        if (flag.eq.4) then
           IF (PINKERRAD1Deight(i).GT.200.0) THEN
            LBPINKeight(i)=.TRUE.
           ELSE
            LBPINKeight(i)=.FALSE.
           ENDIF
        endif
         ENDDO

        if (flag.eq.4) then
        do y=1,npinkeight
        if ((pinkerrad1deight(y).le.200).and.
     &  (lbpinkeight(y).eqv..true.)) then
        print*,'brtm=',y,pinkerrad1deight(y)
        stop
        endif
        enddo
        endif
         CALL IPOLATES (IPPINK,IPPINKOPT,KGDSPINK,KGDSLDAS,NPINKeight,
     &                   NLDAS,KMPINK,IBIPINK,
     +                   LBPINKeight,PINKERRAD1Deight,NO,RLAT,RLON,
     &                   IBO,LO,
     &                   PINKERLDAS1D,IRET)



        ENDIF !endif the data is from 2000 or before

        IF ( (YR.eq.2001).and.((MO.ge.1).AND.(MO.LE.6))) THEN
         KGDSPINK(1)=0
         KGDSPINK(2)=111
         KGDSPINK(3)=51
         KGDSPINK(4)=25000
         KGDSPINK(5)=-125000
         KGDSPINK(6)=128
         KGDSPINK(7)=50000
         KGDSPINK(8)=-70000
         KGDSPINK(9)=500
         KGDSPINK(10)=500
         KGDSPINK(11)=64
         KGDSPINK(12)=0
         KGDSPINK(13)=0
         KGDSPINK(14)=0
         KGDSPINK(15)=0
         KGDSPINK(16)=0
         KGDSPINK(17)=0
         KGDSPINK(18)=0
         KGDSPINK(19)=0
         KGDSPINK(20)=255
         KGDSPINK(21)=0
         KGDSPINK(22)=0
         IPPINK = 0
         IPPINKOPT(1) = 4
         IPPINKOPT(2) = -1
         KMPINK = 1
         IBIPINK = 1
         DO I=1,NPINK
           IF (PINKERRAD1D(I).GT.1.0) THEN
            LBPINK(I)=.TRUE.
           ELSE
            LBPINK(I)=.FALSE.
           ENDIF
        if (flag.eq.4) then
           IF (PINKERRAD1D(I).GT.200.0) THEN
            LBPINK(I)=.TRUE.
           ELSE
            LBPINK(I)=.FALSE.
           ENDIF
	endif
         ENDDO

	if (flag.eq.4) then
	do y=1,npink
	if ((pinkerrad1d(y).le.200).and.(lbpink(y).eqv..true.)) then
	print*,'brtm=',y,pinkerrad1d(y)
	stop
	endif
	enddo
	endif
         CALL IPOLATES (IPPINK,IPPINKOPT,KGDSPINK,KGDSLDAS,NPINK,
     &                   NLDAS,KMPINK,IBIPINK,
     +                   LBPINK,PINKERRAD1D,NO,RLAT,RLON,IBO,LO,
     &                   PINKERLDAS1D,IRET)

	

        ENDIF !endif the data is from Jan 2001 to June 2001 

        IF ( (YR.GE.2002).OR.((YR.EQ.2001).AND.(MO.GE.7))) THEN
         KGDSPINK(1)=0
         KGDSPINK(2)=121
         KGDSPINK(3)=61
         KGDSPINK(4)=24000
         KGDSPINK(5)=-126000
         KGDSPINK(6)=128
         KGDSPINK(7)=54000
         KGDSPINK(8)=-66000
         KGDSPINK(9)=500
         KGDSPINK(10)=500
         KGDSPINK(11)=64
         KGDSPINK(12)=0
         KGDSPINK(13)=0
         KGDSPINK(14)=0
         KGDSPINK(15)=0
         KGDSPINK(16)=0
         KGDSPINK(17)=0
         KGDSPINK(18)=0
         KGDSPINK(19)=0
         KGDSPINK(20)=255
         KGDSPINK(21)=0
         KGDSPINK(22)=0
         IPPINK = 0
         IPPINKOPT(1) = 4
         IPPINKOPT(2) = -1
         KMPINK = 1
         IBIPINK = 1
         DO I=1,NPINKNEW
           IF (PINKERRAD1DNEW(I).GT.1.0) THEN
            LBPINKNEW(I)=.TRUE.
           ELSE
            LBPINKNEW(I)=.FALSE.
           ENDIF
        if (flag.eq.4) then
           IF (PINKERRAD1DNEW(I).GT.200.0) THEN
            LBPINKNEW(I)=.TRUE.
           ELSE
            LBPINKNEW(I)=.FALSE.
           ENDIF
        endif
         ENDDO

        if (flag.eq.4) then
        do y=1,npinknew
        if ((pinkerrad1dnew(y).le.200).and.
     &  (lbpinknew(y).eqv..true.)) then
        print*,'brtm=',y,pinkerrad1dnew(y)
        stop
        endif
        enddo
        endif
         CALL IPOLATES (IPPINK,IPPINKOPT,KGDSPINK,KGDSLDAS,NPINKNEW,
     &                   NLDAS,KMPINK,IBIPINK,
     +                   LBPINKNEW,PINKERRAD1DNEW,NO,RLAT,RLON,IBO,LO,
     &                   PINKERLDAS1D,IRET)
        ENDIF !endif the data is from July 2001 or after


        if (flag.eq.4) then
        do y=1,nldas
        if ((pinkerldas1d(y).le.200).and.(pinkerldas1d(y).gt.
     &  -10000)) then
c        print*,'brtminterp=',y,pinkerldas1d(y)
        endif
        enddo
        endif

C       Transfer 1D data to 2D array
         COUNT = 0
         DO Y = 1, NY
           DO X = 1, NX
             OUTDATA(X,Y) = PINKERLDAS1D(X+COUNT)
           END DO
           COUNT = COUNT + NX
         END DO
         FERROR=2





!=== If Pinker data is not undefined, set FERROR to 1 
	 DO R=1,LDAS%NR
          DO C=1,LDAS%NC
	   IF(OUTDATA(C,R).GT.0.0) FERROR=1
	  ENDDO
	 ENDDO 

  100  	 CONTINUE

	  
!=== If Pinker data is defined, transfer it into the proper
!=== output array.  If it's undefined, fill the output array with
!=== undefined values

	IF (FLAG.EQ.1) THEN	
	 IF (FERROR.EQ.1) THEN
	  DO R=1,LDAS%NR
	   DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
	     GRID(C,R)%PINKDATA1=OUTDATA(C,R)
	    ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%PINKDATA2=OUTDATA(C,R)
	    ENDIF
	   ENDDO
	  ENDDO
	 ELSE
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
             GRID(C,R)%PINKDATA1=LDAS%UDEF
            ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%PINKDATA2=LDAS%UDEF
            ENDIF
           ENDDO
          ENDDO
	 ENDIF
	ENDIF

        IF (FLAG.EQ.4) THEN
         IF (FERROR.EQ.1) THEN
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
             GRID(C,R)%BRTTMPDATA1=OUTDATA(C,R)
            ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%BRTTMPDATA2=OUTDATA(C,R)
            ENDIF
           ENDDO
          ENDDO
         ELSE
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
             GRID(C,R)%BRTTMPDATA1=LDAS%UDEF
            ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%BRTTMPDATA2=LDAS%UDEF
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDIF

	ENDIF  !END the flag.eq.1 or flag.eq.4 loop

!=== If using NESDIS data, open and read in NESDIS GRIB file
        IF (FLAG.EQ.2) THEN

!         CALL UNGRIBRAD(LDAS,GRID,NAME1,LDAS%NC,LDAS%NR,
!     &   LDAS%NCOLD,LDAS%NROLD,
!     &   NINEOCTET,TWEN1OCT,TENOCTET,
!     &   FTYPE,VARADJ,
!     &   VARFIELD,EMPTYVF,ERRORFLAG,GRID%FMASK,DATAFLAG,
!     &   SIXOCTET)

	 FERROR=0

!=== If ungribbing was successful, set FERROR to 2
	 IF (ERRORFLAG.eq.0) then
         FERROR=2

!=== Set FERROR=1 if field is not undefined

         DO R=1,LDAS%NR
          DO C=1,LDAS%NC
           IF(VARFIELD(C,R).GE.0.001) FERROR=1
          ENDDO
         ENDDO


!=== If NESDIS data is defined, transfer it into the proper
!=== output array.  If it's undefined, fill the output array with
!=== undefined values
         IF (FERROR.EQ.1) THEN
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
            IF (GRID(C,R)%NESDATA1.LT.0) GRID(C,R)%NESDATA1=0.0
            IF (GRID(C,R)%NESDATA2.LT.0) GRID(C,R)%NESDATA2=0.0
            IF (ORDER.EQ.1) THEN
             GRID(C,R)%NESDATA1=VARFIELD(C,R)
            ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%NESDATA2=VARFIELD(C,R)
            ENDIF
           ENDDO
          ENDDO
         ELSE
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
             GRID(C,R)%NESDATA1=LDAS%UDEF
            ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%NESDATA2=LDAS%UDEF
            ENDIF
           ENDDO
          ENDDO
         ENDIF
	ENDIF  ! End 'IF ERRORFLAG.eq.0' loop
       ENDIF ! End 'IF FLAG.EQ.2' loop





!=== If using BRTTMP data, open and read in NESDIS GRIB file
        IF (FLAG.EQ.3) THEN
!=== Set GRIB octet values to those suitable for brightness temperature
!=== Initialize variables
        SIXOCTET=154
        NINEOCTET=118
        TENOCTET=1
        TWEN1OCT=0
        VARADJ=0
        FERROR=0
!         CALL UNGRIBRAD(LDAS,GRID,NAME1,LDAS%NC,LDAS%NR,
!     &   LDAS%NCOLD,LDAS%NROLD,
!     &   NINEOCTET,TWEN1OCT,TENOCTET,
!     &   FTYPE,VARADJ,
!     &   VARFIELD,EMPTYVF,ERRORFLAG,GRID%FMASK,DATAFLAG,
!     &   SIXOCTET)
         FERROR=0

!=== If ungribbing was successful, set FERROR to 2
         IF (ERRORFLAG.eq.0) then
         FERROR=2
!=== Set FERROR=1 if field is not undefined


         DO R=1,LDAS%NR
          DO C=1,LDAS%NC
           IF(VARFIELD(C,R).GE.0.001) FERROR=1
          ENDDO
         ENDDO

!=== If BRTTMP data is defined, transfer it into the proper
!=== output array.  If it's undefined, fill the output array with
!=== undefined values
         IF (FERROR.EQ.1) THEN
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
            IF (GRID(C,R)%BRTTMPDATA1.LT.0) GRID(C,R)%BRTTMPDATA1=
     &                   LDAS%UDEF
            IF (GRID(C,R)%BRTTMPDATA2.LT.0) GRID(C,R)%BRTTMPDATA2=
     &                   LDAS%UDEF
            IF (ORDER.EQ.1) THEN
             GRID(C,R)%BRTTMPDATA1=VARFIELD(C,R)
            ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%BRTTMPDATA2=VARFIELD(C,R)
            ENDIF
           ENDDO
          ENDDO
         ELSE
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
            IF (ORDER.EQ.1) THEN
             GRID(C,R)%BRTTMPDATA1=LDAS%UDEF
            ELSEIF (ORDER.EQ.2) THEN
             GRID(C,R)%BRTTMPDATA2=LDAS%UDEF
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDIF  ! End 'IF ERRORFLAG.eq.0' loop
       ENDIF ! End 'IF FLAG.EQ.2' loop




	RETURN
	END

