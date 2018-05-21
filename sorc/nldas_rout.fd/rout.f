      PROGRAM rout

      IMPLICIT NONE

C     Routing Code Version 2.0, Dag Lohmann, 03/16/02
C     READ NOAH GRIB2 OUTPUT AND OUTPUT GRIB2 STREAMFLOW FILE
C     Modified byYoulong Xia, 11/07/13
C
C     NX    -- grid points in west-east direction
C     NY    -- grid points in south-north direction
C     DT    -- time step in seconds 
C     DUH   -- length of the internal unit-hydrograph in days
C     LUH   -- length of the internal unit-hydrograph in DT
C     MAX_T -- maximum travel time of water through one grid box in seconds
C     LTR   -- length of transport unit-hydrograph in DT

      INTEGER NX
      INTEGER NY
      INTEGER DUH
      INTEGER LUH
      INTEGER DT
      INTEGER MAX_T
      INTEGER LTR
      INTEGER NLDAS
      INTEGER HOUR, HOUR1
      INTEGER jret
      
      PARAMETER (NX    = 464)
      PARAMETER (NY    = 224)
      PARAMETER (NLDAS = NX*NY)
      PARAMETER (DT    = 3600)
      PARAMETER (DUH   = 2)
      PARAMETER (LUH   = DUH*86400/DT)
      PARAMETER (MAX_T = 8*3600)
      PARAMETER (LTR   = MAX_T / DT)

      INTEGER ORDER(4,NX*NY)

      LOGICAL*1 LB(NX,NY)
      LOGICAL*1 LDASMASK(NX,NY)
      REAL SURFACE_RUNOFF(NX,NY)
      REAL BASEFLOW(NX,NY)
      REAL PRECIP(NX,NY)
      REAL UH_INTERN(LUH,NX,NY)
      REAL UH_TRANS(LTR,NX,NY)
      REAL RUNOFF_INTERN(LUH,NX,NY)
      REAL RUNOFF_TRANS(LTR,NX,NY)
      REAL RUNOFF_INTERN_sum(NX,NY)
      REAL RUNOFF_TRANS_sum(NX,NY)
      REAL STREAMFLOW(NX,NY)
      REAL AREA(NX,NY)
      
      INTEGER LAND_SEA(NX,NY) 
      INTEGER JPDT1(200), JPDT2(200)      

      CHARACTER*8 TODAY, YESTERDAY
      CHARACTER*8 DAYS, DAYS0
      CHARACTER*4 YEARS, YEARS0
      CHARACTER*11  CNTRFL
      CHARACTER*150 initial_1, initial_2, order_file, uh_file_1, uh_file_2
      CHARACTER*150 filename
      CHARACTER*11 ENVVAR
      
      CHARACTER*100 FILEOUT, routing_out
      CHARACTER*2   hourchar
      CHARACTER*2   hourchar1 
      INTEGER       I, J, K, length, TIME, ORDER_N
C     READ LDAS MASK
      OPEN(30, FILE='UMDunifiedmask19990327.bin',FORM='UNFORMATTED')
      READ(30) LAND_SEA
      CLOSE(30)

C     DEFINE THE OUTPUT MASK FOR GRIB FILES.
      LDASMASK = .TRUE.
      WHERE (LAND_SEA .EQ. 0)
         LDASMASK = .FALSE.
      END WHERE

      OPEN(10,FILE='directories.txt')
      READ(10,*) initial_1
      READ(10,*) initial_2
      READ(10,*) order_file
      READ(10,*) uh_file_1
      READ(10,*) uh_file_2
      WRITE(*,*) initial_1, initial_2, order_file, uh_file_1, uh_file_2

      CALL read_real_spatial(uh_file_1,LUH,NX,NY,UH_INTERN)
      CALL read_real_spatial(uh_file_2,LTR,NX,NY,UH_TRANS)
      CALL read_real_spatial(initial_1,LUH,NX,NY,RUNOFF_INTERN)
      CALL read_real_spatial(initial_2,LTR,NX,NY,RUNOFF_TRANS)
      CALL read_int_spatial(order_file,4,NX*NY,1,ORDER)

      DO I = 1, NX*NY
         IF (ORDER(1,I) .EQ. 0) THEN
            ORDER_N = I - 1
            EXIT
         END IF
      END DO

      CALL CALC_AREA(NX,NY,AREA,DT)
      OPEN(17, FILE = 'days0.txt')
      READ(17,*) years0, days0

      OPEN(17, FILE = 'days.txt')
       READ(17,*) years, days

      do i=1,200
      jpdt1(i)=-9999 ! initialize for this parameter for wildcard
      jpdt2(i)=-9999
      enddo
C    for surface runoff
      jpdt1(1)=0
      jpdt1(2)=6
      jpdt1(10)=1     
      jpdt1(11)=0
      jpdt1(12)=0
C     for baseflow
      jpdt2(1)=0
      jpdt2(2)=5
      jpdt2(10)=1     
      jpdt2(11)=0
      jpdt2(12)=0

         DO j = 0, 23
            HOUR=J+1
            HOUR1=J
            WRITE(hourchar,'(I2.2)') HOUR
CJZ            filename = 'fort.'//hourchar//char(0)
            ENVVAR = 'FORT'//hourchar//char(0)
            CALL GET_ENVIRONMENT_VARIABLE(ENVVAR,FILENAME)
            TODAY = days
            YESTERDAY = days0
C            WRITE(*,*) filename
CYX            CALL READ_GRIB(NLDAS,SURFACE_RUNOFF,FILENAME,LB,59,235,0)
CYX            CALL READ_GRIB(NLDAS,BASEFLOW,FILENAME,LB,59,234,0)
            CALL READ_NOAH_GRIB2(NLDAS,SURFACE_RUNOFF,FILENAME,59,jpdt1)       
            CALL READ_NOAH_GRIB2(NLDAS,BASEFLOW,FILENAME,59,jpdt2)
c            write(*,'(f7.4)') SURFACE_RUNOFF(375,108)
c            write(*,'(f7.4)') BASEFLOW(375,108)
            CALL ROUT_IT(NX,NY,LUH,LTR,SURFACE_RUNOFF,BASEFLOW,
     &           STREAMFLOW,RUNOFF_INTERN,RUNOFF_TRANS,UH_INTERN,UH_TRANS,
     &           ORDER,ORDER_N,AREA)

C -------------- OUTPUT HOURLY STREAMFLOW GRIB FILE ----------------------    
         WHERE (STREAMFLOW .LT. 0)
         LDASMASK = .FALSE.
         END WHERE
    
         CALL GRIBOUT(NX,NY,LDASMASK,HOUR1,YESTERDAY,TODAY,STREAMFLOW)
         END DO

      WRITE(*,*) 'END OF HOUR LOOP'
C -------- OUTPUT RIVER ROUTING MODEL INITIALS ---------------------------- 
      WRITE(81) RUNOFF_INTERN
      WRITE(82) RUNOFF_TRANS
      END


