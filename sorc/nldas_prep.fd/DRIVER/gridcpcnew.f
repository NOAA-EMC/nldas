	SUBROUTINE GRIDCPCNEW (HIGGINSNAME,CPC,NLDAS,FLAG)
C
C***********************************************************************
C  PROGRAM:  GRIDCPC                REGRID CPC PREC TO LDAS GRID
C  PRGMMR:  MARSHALL                ORG:  W/NP20
C
C  ABSTRACT:  THIS PROGRAM REGRIDS THE CPC (WAYNE HIGGIN'S) GAUGE-ONLY
C  DAILY (12Z - 12Z) PRECIP PRODUCT, TO THE 1/8 DEGREE LDAS GRID. 
C  REGRIDDING IS PERFORMED BY DOING A DIRECT ASSIGNMENT OF A 1/4 DEGREE
C  PIXEL TO THE 4 1/8 DEGREE PIXELS CONTAINED THEREIN. THE GRID IS 
C  RETURNED TO THE CALLING PROGRAM, PRECFORC.f, WHERE IT IS USED IN
C  CONJUNCTION WITH WEIGHTS FROM STAGE IV TO COMPUTE HOURLY PRECIPS. 
C
C  PROGRAM HISTORY LOG:
C  98-12-30  MARSHALL  ORIGINAL CODING
C  00-02-28  Cosgrove  Added 'FLAG' variable to program to indicate
C                      whether successful precip retrieval occurred
C  00-11-14  Cosgrove  Changed LOGICAL *1 to LOGICAL
C  01-03-08  Cosgrove  Added ERR=77 line to READ statement
C  02-05-15  Jambor    Changed LOGICAL to LOGICAL*1 to match new 
C                      GRIB libraries
C  12-04-18 Xia        Used 0.5 degree CPC unified global precipitation
C***********************************************************************
C
C  DECLARATIONS
C
	IMPLICIT NONE
C
	INTEGER I, NCPC, NLDAS, LDASGRD, LDASGDS(22), MISG, 
     +          GDS(22), LENGDS, IRET, CPCGRD, CPCGDS(22), IP,
     +          IPOPT(20), KM, IBI, NO, IBO,  LUGB, NN
C
	INTEGER FLAG           !0=Failed to get CPC Precip
			       !1=Successful opening of file

	PARAMETER (NCPC = 720*360)
C
	REAL DUMM(NCPC)
	REAL RAWPREC(NCPC), CPC(NLDAS), RLAT(NLDAS), RLON(NLDAS)
C
	LOGICAL*1 RAWBIT(NCPC), LO(NLDAS)
	CHARACTER*80 HIGGINSNAME
C
C  DATA STATEMENTS AND CONSTANTS
C
        FLAG=0
	KM = 1
	IBI = 1
	LDASGRD = 110
c	CPCGRD = 238
	IP = 3
!         DO I=1,1
!            IPOPT(I)=17
!         ENDDO
!         DO I=2,20
!            IPOPT(I)=-1
!         ENDDO

         DO I=1,20
            IPOPT(I)=-1
         ENDDO

         IPOPT(1)=4

C
C  SETUP GDS OCTETS FOR INPUT(CPC) AND OUTPUT(LDAS) GRIDS
C
	CALL MAKGDS (LDASGRD, LDASGDS, GDS, LENGDS, IRET)
	IF (IRET .NE. 0) THEN
	  PRINT*, 'MAKLDASGDS FAILED, IRET= ',IRET
	  STOP
	END IF
C
C	write(*,*) LDASGDS
        DO I=1,22
         CPCGDS(I)=0
        ENDDO
         CPCGDS(1)=0
         CPCGDS(2)=720
         CPCGDS(3)=360
         CPCGDS(4)=-89750
         CPCGDS(5)=000
         CPCGDS(6)=128  
         CPCGDS(7)=89750
         CPCGDS(8)=360000
         CPCGDS(9)=500
         CPCGDS(10)=500
         CPCGDS(11)=64
         CPCGDS(12)=0
         CPCGDS(13)=0
         CPCGDS(14)=0
         CPCGDS(15)=0
         CPCGDS(16)=0
         CPCGDS(17)=0
         CPCGDS(18)=0
         CPCGDS(19)=0
         CPCGDS(20)=255
         CPCGDS(21)=0
         CPCGDS(22)=0

C initialize CPC data
        DO I=1,NLDAS
        CPC(I)=-999.0
        ENDDO
C
C  READ IN THE ORIGINAL GRID AND SETUP A BITMAP FOR IPOLATES.  IF VALUE
C  LESS THAN 0.254 mm, THEN RESET TO ZERO.  
C
C direct access read
C        open (unit=59,file=HIGGINSNAME,
C     &  FORM='UNFORMATTED',
C     &  STATUS='OLD',access='direct',recl=321*201*4,err=77)
C        read (59,rec=1) rawprec
C        CLOSE(59) 
C sequential read
        open (unit=59,file=HIGGINSNAME,FORM='UNFORMATTED',err=77) 
        read (59) rawprec
        CLOSE(59)
C ---------- see if all data are undefined ------------------------
        NN=0
        DO I = 1, NCPC
        IF(rawprec(I).GT.0.0) NN=NN+1
        ENDDO

        IF(NN.NE.0) THEN
C --- convert the data into mm -------------------------------------
 	DO I = 1, NCPC
       	  RAWBIT(I) = .TRUE.
          RAWPREC(I) = RAWPREC(I)/10.0

          IF (RAWPREC(I) .LT. 0.1) THEN
            RAWPREC(I) = 0
          END IF
   
	END DO
C
C  REGRID THE DATA WITH IPOLATES   
C
	CALL IPOLATES (IP,IPOPT,CPCGDS,LDASGDS,NCPC,NLDAS,
     +                 KM,IBI,RAWBIT,RAWPREC,NO,RLAT,RLON,
     +                 IBO,LO,CPC,IRET)
        
        IF (IRET .EQ. 0) FLAG=1
        ENDIF       ! not all CPC P values are undefined
 77     continue
C	PRINT *,'flag=',flag
	RETURN
	END

