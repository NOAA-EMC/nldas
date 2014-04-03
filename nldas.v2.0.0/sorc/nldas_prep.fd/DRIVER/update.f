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
! update.f:
!
! DESCRIPTION:
!  Updates PDS information (dates?)
!
!
! REVISION HISTORY:
!  Curtis Marshall, Original Code From NCEP
!  18 Feb 2000: Brian Cosgrove; Used in conjunction with merged precip 
!               product generation
!  15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match
!               new GRIB libraries
!=========================================================================
        



	SUBROUTINE UPDATE (KPDS, LEAP)
C
        INTEGER KPDS(25)
        LOGICAL*1 LEAP
C
        IF (KPDS(11).GE.24) THEN
            KPDS(10)=KPDS(10)+1
            KPDS(11)=KPDS(11)-24
            IF ((((KPDS(9).EQ.9).OR.(KPDS(9).EQ.4).OR.
     +          (KPDS(9).EQ.6).OR.(KPDS(9).EQ.11)).AND.
     +          (KPDS(10).EQ.31)) .OR. (KPDS(10).EQ.32).OR.
     +          ((KPDS(9).EQ.2).AND.(KPDS(10).EQ.29).AND.
     +           (.NOT. LEAP)).OR.((KPDS(9).EQ.2).AND.
     +           (KPDS(10).EQ.30).AND.(LEAP))) THEN
                KPDS(9)=KPDS(9)+1
                KPDS(10)=1
                IF (KPDS(9).EQ.13) THEN
                    KPDS(8)=KPDS(8)+1
                    KPDS(9)=1
                    IF (KPDS(8).EQ.101) THEN
                        KPDS(21)=KPDS(21)+1
                        KPDS(8)=1
                    END IF
                END IF
            END IF
        END IF
C
        RETURN
        END

