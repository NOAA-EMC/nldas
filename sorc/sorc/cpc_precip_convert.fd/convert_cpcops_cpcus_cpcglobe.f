      PROGRAM  MAIN

C     Convert CPC operational precipitation to CPC NLDAS domain
C     Youlong Xia, EMC, July 2011            
C    C  DECLARATIONS
C
	IMPLICIT NONE
C -------- for daily cpc us precipitation (0.125 deg resolution)----------   
	INTEGER I, NCPC, NLDAS, LDASGRD, LDASGDS(22), MISG, 
     +          GDS(22), LENGDS, IRET, CPCGRD, CPCGDS(22), IP,
     +          IPOPT(20), KM, IBI, NO, IBO,  LUGB
C
	INTEGER FLAG           !0=Failed to get CPC Precip
			       !1=Successful opening of file

	PARAMETER (NCPC = 601*241)
        PARAMETER (NLDAS= 464*224)
C ------ for daily cpc global precipitation (0.25 deg resolution) ---------
        INTEGER MCPC
        
        PARAMETER (MCPC=720*360)
        
        REAL RAWGP(MCPC), DUMMG(MCPC)
C ---------------- for NLDAS US domain ---------------------------------- 
        REAL DUMM(NCPC)
	REAL RAWPREC(NCPC), CPC(NLDAS), RLAT(NLDAS), RLON(NLDAS)
C
	LOGICAL*1 RAWBIT(NCPC), LO(NLDAS)
        character*48 filename
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
          CALL ERREXIT(IRET)
	  STOP
	END IF
C
 
        DO I=1,22
         CPCGDS(I)=0
        ENDDO
         CPCGDS(1)=0
         CPCGDS(2)=601
         CPCGDS(3)=241
         CPCGDS(4)=20000
         CPCGDS(5)=-130000
         CPCGDS(6)=128  
         CPCGDS(7)=50125
         CPCGDS(8)=-54875
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
	
C	CALL MAKGDS (CPCGRD, CPCGDS, GDS, LENGDS, IRET)
C	IF (IRET .NE. 0) THEN
C	  PRINT*, 'MAKCPCGDS FAILED, IRET= ',IRET
C	  STOP
C	END IF
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
        
      open(59,file='cpcops.bin',access='direct',recl=601*241,
     &convert='little_endian')
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
	IF (IRET .NE. 0) THEN
	  PRINT*, 'IPOLATES CPC FAILED, IRET= ',IRET
          CALL ERREXIT(IRET)
	  STOP
	END IF
C
	FLAG=1
 77     continue
C	PRINT *,'flag=',flag
C ----------------- CONVERT 0.1MM INTO INCH -------------------------------
        DO I=1,NLDAS
        IF(CPC(I).GE.0.0) THEN
        CPC(I)=CPC(I)*0.1/25.4
        ENDIF
        ENDDO
       
        open(20,file='cpcldas.bin',form='unformatted') 
        write(20) CPC
        close(20)
C --------------------- daily cpc us done! -----------------------------------
        open (unit=59,file='cpc_global.bin',access='direct',
     &recl=720*360,convert='little_endian')
        read (59,rec=1) RAWGP
        read (59,rec=2) DUMMG
        CLOSE(59)
        
        open(59,file='cpc_globe.bin',form='unformatted')
        write(59) RAWGP
        close(59)

	stop
	END

