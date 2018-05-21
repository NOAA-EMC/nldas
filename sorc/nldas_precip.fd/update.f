        SUBROUTINE UPDATE (KPDS, LEAP)
C
        INTEGER KPDS(25)
        LOGICAL LEAP
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

