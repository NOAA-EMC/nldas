      SUBROUTINE GRIDCPC (CPC,NLDAS)

C***********************************************************************
C  PROGRAM:  GRIDCPC                REGRID CPC PREC TO LDAS GRID
C  PRGMMR:   YOULONG XIA            ORG:  EMC Land Group
C
C  ABSTRACT:  THIS PROGRAM REGRIDS THE CPC OPERATIONAL GLOBAL GAUGE-ONLY
C  DAILY (12Z - 12Z) PRECIP PRODUCT (1/8TH DEGREE), TO THE NLDAS DOAMIN.
C  THE GRID IS RETURNED TO THE CALLING PROGRAM, PRECFORC.f, WHERE IT IS USED
C  IN CONJUNCTION WITH WEIGHTS FROM STAGE IV TO COMPUTE HOURLY PRECIPS.
C
C  PROGRAM HISTORY LOG:
C  2015-12-15  YOULONG XIA  ORIGINAL CODING
C***********************************************************************
            
	IMPLICIT NONE
C read cpc global daily 0.125 degree precipitation and cut NLDAS domain  
	INTEGER I, NCPC, NLDAS, LDASGRD, LDASGDS(22), MISG, 
     +          GDS(22), LENGDS, IRET, CPCGRD, CPCGDS(22), IP,
     +          IPOPT(20), KM, IBI, NO, IBO,  LUGB
C
	PARAMETER (NCPC = 2881*1441)
C        PARAMETER (NLDAS= 464*224)
C ---------------- for NLDAS US domain ---------------------------------- 
        REAL DUMM(NCPC)
	REAL RAWPREC(NCPC), CPC(NLDAS), RLAT(NLDAS), RLON(NLDAS)
        REAL CPCNLDAS(NLDAS) 
C
	LOGICAL*1 RAWBIT(NCPC), LO(NLDAS)
        LOGICAL THERE
        character*48 filename
C
C  DATA STATEMENTS AND CONSTANTS
C
C  INITIALIZE INTERPOLATED CPC PRECIPITTAION  
        DO I=1,NLDAS
        CPC(I)=-999.0
        ENDDO

	KM = 1
	IBI = 1
	LDASGRD = 110
	IP = 3

         DO I=1,20
            IPOPT(I)=-1
         ENDDO

         IPOPT(1)=4

C
C  SETUP GDS OCTETS FOR INPUT(CPC) AND OUTPUT(LDAS) GRIDS
C
	CALL MAKGDS (LDASGRD, LDASGDS, GDS, LENGDS, IRET)
	IF (IRET .NE. 0) THEN
	  PRINT*, 'MAKLDASGDS FAILED for gridcpc.f,IRET= ',IRET
	  STOP
	ENDIF
C
        DO I=1,22
         CPCGDS(I)=0
        ENDDO
         CPCGDS(1)=0
         CPCGDS(2)=2881
         CPCGDS(3)=1441
         CPCGDS(4)=-90000
         CPCGDS(5)=0
         CPCGDS(6)=128  
         CPCGDS(7)=90000
         CPCGDS(8)=360000
         CPCGDS(9)=125
         CPCGDS(10)=125
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
C
C  READ IN THE ORIGINAL GRID AND SETUP A BITMAP FOR IPOLATES.  IF VALUE
C  LESS THAN 0.01 mm, THEN RESET TO ZERO.  

C sequential read
c inquiry if daily gauge precipitation exists
       
        INQUIRE( FILE='cpc_125global.bin', EXIST=THERE )
        IF ( THERE ) THEN        
        open(59,file='cpc_125global.bin', access='direct',recl=2881*1441,convert='little_endian',err=77)
        read (59,rec=1) rawprec
        read (59,rec=2) dumm
        CLOSE(59)

 	DO I = 1, NCPC
       	  RAWBIT(I) = .TRUE.
          RAWPREC(I) = RAWPREC(I)*1.0   
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
	IF (IRET .EQ. 0) THEN
C --------------------- daily cpc precip done! -----------------------
        DO I=1,NLDAS
        IF(CPC(I).GE.0.0) THEN
        CPC(I)=CPC(I)*0.1/25.4
        ENDIF
        ENDDO
        ENDIF       ! data are sucessully interpolated 
       
        ENDIF       ! file exists

  77    CONTINUE      

	RETURN
	END

