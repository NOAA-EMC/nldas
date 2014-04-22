	SUBROUTINE GRIDCPC (HIGGINSNAME,CPC,NLDAS,FLAG)
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
C***********************************************************************
C
C  DECLARATIONS
C
	IMPLICIT NONE
C
	INTEGER I, NCPC, NLDAS, LDASGRD, LDASGDS(22), MISG, 
     +          GDS(22), LENGDS, IRET, CPCGRD, CPCGDS(22), IP,
     +          IPOPT(20), KM, IBI, NO, IBO,  LUGB
C
	INTEGER FLAG           !0=Failed to get CPC Precip
			       !1=Successful opening of file

	PARAMETER (NCPC = 4131)
C

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


        DO I=1,22
         CPCGDS(I)=0
        ENDDO
         CPCGDS(1)=0
         CPCGDS(2)=81
         CPCGDS(3)=51
         CPCGDS(4)=10000
         CPCGDS(5)=-140000
         CPCGDS(6)=128  
         CPCGDS(7)=60000
         CPCGDS(8)=-60000
         CPCGDS(9)=1000
         CPCGDS(10)=1000
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
	
!	CALL MAKGDS (CPCGRD, CPCGDS, GDS, LENGDS, IRET)
!	IF (IRET .NE. 0) THEN
!	  PRINT*, 'MAKCPCGDS FAILED, IRET= ',IRET
!	  STOP
!	END IF
C
C  READ IN THE ORIGINAL GRID AND SETUP A BITMAP FOR IPOLATES.  IF VALUE
C  LESS THAN 0.254 mm, THEN RESET TO ZERO.  
C

! direct access read
!        open (unit=59,file=HIGGINSNAME,
!     &  FORM='UNFORMATTED',
!     &  STATUS='OLD',access='direct',recl=81*51*4,err=77)
!        read (59,rec=1) rawprec
!	CLOSE(59)

! sequential read
        open (unit=59,file=HIGGINSNAME,
     &  FORM='UNFORMATTED',
     &  STATUS='OLD',err=77)
        read (59) rawprec
        CLOSE(59)

 	DO I = 1, NCPC
       	  RAWBIT(I) = .TRUE.
C          RAWPREC(I) = RAWPREC(I)*25.4    
          IF (RAWPREC(I) .LT. 0.254) THEN
            RAWPREC(I) = 0
          END IF
	END DO
C
C  REGRID THE DATA WITH IPOLATES   
C
	CALL IPOLATES (IP,IPOPT,CPCGDS,LDASGDS,NCPC,NLDAS,
     +                 KM,IBI,RAWBIT,RAWPREC,NO,RLAT,RLON,
     +                 IBO,LO,CPC,IRET)
	IF (IRET .NE. 0) THEN
	  PRINT*, 'IPOLATES CPC FAILED, IRET= ',IRET
	  STOP
	END IF
C
	FLAG=1
 77     continue
c	PRINT *,'flag=',flag
	RETURN
	END

