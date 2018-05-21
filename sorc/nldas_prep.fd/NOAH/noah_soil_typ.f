      SUBROUTINE SOILTYPE(NC,NR,SAND,CLAY,SOILTYP)

!  This subroutine uses the percentages of sand and clay   
!  derived from the global soils dataset of Reynolds, 
!  Jackson, and Rawls [1999], to convert to Zobler
!  soil class values to be used in NOAH LSM v2.5 in LDAS.
!   (Original code by Matt Rodell, 3/7/01)
!
!  28 Apr 2002, K Arsenault:  Added NOAH LSM to LDAS

      implicit none

      integer :: NC,NR
      integer :: SOILTYP(NC,NR)
      integer :: I,J

      real :: SA,CL
      real :: SAND(NC,NR),CLAY(NC,NR)

      DO J=1,NR
        DO I=1,NC

!     Test for ocean points.
         IF (CLAY(I,J) .lt. 0.00) THEN
           SOILTYP(I,J) = -99
         ELSE
           CL = CLAY(I,J)
           SA = SAND(I,J)
         ENDIF

!     Identify texture class.

         IF (CL .lt. 0.23) THEN
            IF (SA .lt. 0.50) THEN
              SOILTYP(I,J) = 8		! Loam
            ELSE
              IF (SA .lt. 0.75) THEN
                SOILTYP(I,J) = 4        ! Sandy Loam
              ELSE
                SOILTYP(I,J) = 1        ! Loamy Sand
              END IF
            END IF
!--------
         ELSE 

            IF (CL .lt. 0.28) THEN
              IF (SA .lt. 0.45) THEN
                SOILTYP(I,J) = 8        ! Loam
              ELSE
                SOILTYP(I,J) = 7        ! Sandy Clay Loam
              ENDIF
!-----------
            ELSE

              IF (CL .lt. 0.37) THEN
                IF (SA .lt. 0.2) THEN
                  SOILTYP(I,J) = 2      ! Silty Clay Loam
                ELSE
                  IF (SA .lt. 0.43) THEN
                    SOILTYP(I,J) = 6    ! Clay Loam
                  ELSE
                    SOILTYP(I,J) = 7    ! Sandy Clay Loam
                  END IF
                END IF
!-------------
              ELSE

                IF (CL .lt. 0.41) THEN
                  IF (SA .lt. 0.2) THEN
                    SOILTYP(I,J) = 2   ! Silty Clay Loam
                  ELSE
                     IF (SA .lt. 0.43) THEN
                       SOILTYP(I,J) = 6    ! Clay Loam
                     ELSE
                       SOILTYP(I,J) = 5    ! Sandy Clay
                     END IF
                  END IF

                ELSE

                  IF (SA .lt. 0.43) THEN
                     SOILTYP(I,J) = 3      ! Light Clay
                  ELSE
                     SOILTYP(I,J) = 5      ! Sandy Clay
                  END IF

                END IF
              END IF
            END IF
         END IF

!- End do loops 

         END DO
       END DO

       END

