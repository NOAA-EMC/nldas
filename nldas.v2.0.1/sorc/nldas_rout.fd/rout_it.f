      SUBROUTINE ROUT_IT(NX,NY,LUH,LTR,SURFACE_RUNOFF,BASEFLOW,
     &     STREAMFLOW,RUNOFF_INTERN,RUNOFF_TRANS,UH_INTERN,UH_TRANS,
     &     ORDER,ORDER_N,AREA)

      IMPLICIT NONE

C     NX    -- grid points in west-east direction
C     NY    -- grid points in south-north direction
C     DT    -- time step in seconds 
C     LUH   -- length of the internal unit-hydrograph in DT
C     LTR   -- length of transport unit-hydrograph in DT

      INTEGER NX
      INTEGER NY
      INTEGER LUH
      INTEGER LTR
      INTEGER NOB
      INTEGER ORDER_N

      INTEGER ORDER(4,NX*NY)

      REAL SURFACE_RUNOFF(NX,NY)
      REAL BASEFLOW(NX,NY)
      REAL STREAMFLOW(NX,NY)
      REAL UH_INTERN(LUH,NX,NY)
      REAL UH_TRANS(LTR,NX,NY)
      REAL RUNOFF_INTERN(LUH,NX,NY)
      REAL RUNOFF_TRANS(LTR,NX,NY)
      REAL RUNOFF_IN(NX,NY)
      REAL AREA(NX,NY)
      REAL DT

      INTEGER N, I, J, IX, IY, IXX, IYY

      DO J = 1, NY
         DO I = 1, NX
            BASEFLOW(I,J) = BASEFLOW(I,J) * AREA(I,J)/3.6
            SURFACE_RUNOFF(I,J) = SURFACE_RUNOFF(I,J) * AREA(I,J)/3.6
         END DO
      END DO

C      WRITE(*,*) 'ORDER_N = ', ORDER_N
      DO N = 1,ORDER_N       
         IX  = ORDER(1,N)
         IY  = ORDER(2,N)
         RUNOFF_IN(IX,IY) = 0.0
      END DO
      STREAMFLOW = -9999.0

      DO N = 1,ORDER_N         
         IX  = ORDER(1,N)
         IY  = ORDER(2,N)
         IXX = ORDER(3,N)
         IYY = ORDER(4,N)

         IF (BASEFLOW(IX,IY) .LT. 0.0) THEN
C            WRITE(*,*) IX, IY,BASEFLOW(IX,IY) 
            BASEFLOW(IX,IY) = 0.0
         END IF

         IF (SURFACE_RUNOFF(IX,IY) .LT. 0.0) THEN
            WRITE(*,*) IX, IY, SURFACE_RUNOFF(IX,IY)
            SURFACE_RUNOFF(IX,IY) = 0.0
         END IF

         DO I = 1,LUH
            RUNOFF_INTERN(I,IX,IY) = RUNOFF_INTERN(I,IX,IY)
     &           + (SURFACE_RUNOFF(IX,IY) + BASEFLOW(IX,IY))
     &           * UH_INTERN(I,IX,IY)
         END DO
         RUNOFF_IN(IXX,IYY) = RUNOFF_IN(IXX,IYY) + 
     &                        RUNOFF_INTERN(1,IX,IY)

         DO I = 1,LTR
            RUNOFF_TRANS(I,IX,IY) = RUNOFF_TRANS(I,IX,IY) + 
     &                              UH_TRANS(I,IX,IY) *
     &                              RUNOFF_IN(IX,IY)
         END DO
         RUNOFF_IN(IXX,IYY) = RUNOFF_IN(IXX,IYY) +  
     &                        RUNOFF_TRANS(1,IX,IY)
      END DO

      DO N = 1,ORDER_N      
         IX  = ORDER(1,N)
         IY  = ORDER(2,N)
         STREAMFLOW(IX,IY) = RUNOFF_TRANS(1,IX,IY) + RUNOFF_INTERN(1,IX,IY)

         DO I = 2,LUH
            RUNOFF_INTERN(I-1,IX,IY) = RUNOFF_INTERN(I,IX,IY)
         END DO
         RUNOFF_INTERN(LUH,IX,IY) = 0.0

         DO I = 2,LTR
            RUNOFF_TRANS(I-1,IX,IY) = RUNOFF_TRANS(I,IX,IY)
         END DO
         RUNOFF_TRANS(LTR,IX,IY) = 0.0

      END DO

      END



