!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          MAKEPDS.F
!
!       3-19-02  Cosgrove; Altered code for use with CLM writing
!                intervals
!       5-17-02  Cosgrove; Changed code to use ldas%lsm instead of
!                ldas%rmos and ldas%rclm so that forcing generation
!                would still work when rmos set to zero in card 
!       28-5-02  Arsenault; Altered code for use with NOAH writing
!                intervals
!       03-2-03  Gottschalck; Altered code for CLM v2.0 writing      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      SUBROUTINE MAKEPDS(LDAS,YESTERDAY, BEFOREYESTER, KPDS, HOUR)
      use ldas_module      ! LDAS non-model-specific 1-D variables
      IMPLICIT NONE
      TYPE (LDASDEC) LDAS

      CHARACTER*8 YESTERDAY, BEFOREYESTER
      INTEGER     KPDS(25), HOUR
     
C     SET TIME-RELATED KPDS OCTETS
      IF (KPDS(16) .NE. 0) THEN
            if (ldas%lsm.eq.1) then
              KPDS(11) = HOUR - LDAS%writeintm
            elseif (ldas%lsm.eq.2) then
              KPDS(11) = HOUR - LDAS%writeintc1
            elseif (ldas%lsm.eq.4) then
              KPDS(11) = HOUR - LDAS%writeintn
            elseif (ldas%lsm.eq.5) then
              KPDS(11) = HOUR - LDAS%writeintc2
            endif

            IF (KPDS(11).LT.0) THEN
            if (ldas%lsm.eq.1) then
              KPDS(11)=24-LDAS%writeintm
            elseif (ldas%lsm.eq.2) then
              KPDS(11)=24-LDAS%writeintc1
            elseif (ldas%lsm.eq.4) then
              KPDS(11)=24-LDAS%writeintn
            elseif (ldas%lsm.eq.5) then
              KPDS(11)=24-LDAS%writeintc2
	    endif
              READ (BEFOREYESTER,'(4(I2))') KPDS(21), KPDS(8),
     &        KPDS(9), KPDS(10)
            ELSE
              READ (YESTERDAY,'(4(I2))') KPDS(21), KPDS(8),
     &        KPDS(9), KPDS(10)
            ENDIF
      ELSE
            READ (YESTERDAY,'(4(I2))') KPDS(21), KPDS(8),
     &                                 KPDS(9), KPDS(10)
            KPDS(11) = HOUR
      END IF
  
      IF (KPDS(8) .EQ. 0) THEN
         KPDS(8) = 100
      ELSE
         KPDS(21) = KPDS(21) + 1
      END IF

      END
