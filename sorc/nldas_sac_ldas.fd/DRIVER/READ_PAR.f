      SUBROUTINE READ_PAR(NX,NY,NMONTH,SOILTYP,LAND_SEA,VEGTYP,
     &                    SLOPETYP,SOILDEPTH,ALBEDO,SHDFAC,ELEV_2d,
C                SAC Parameters
     &                    PECLIM,PERATIO,PEADJ,
     &                    UZTWM_2d,UZFWM_2d,UZK_2d,PCTIM_2d,
     &                    ADIMP_2d,RIVA_2d,ZPERC_2d,REXP_2d,
     &                    LZTWM_2d,LZFSM_2d,LZFPM_2d,LZSK_2d,
     &                    LZPK_2d,PFREE_2d,SIDE_2d,RSERV_2d,
C                SNOW17 Parameters
     &                    SCF_2d,MFMAX_2d,MFMIN_2d,UADJ_2d,SI_2d,
     &                    NMF_2d,TIPM_2d,MBASE_2d,PXTEMP_2d,
     &                    PLWHC_2d,DAYGM_2d,ADC_2d,NSOLD,PRCP_MASK)

C      IMPLICIT NONE
C     READS INPUT PARAMETERS FOR NOAH - LSM
C     
      INTEGER NX
      INTEGER NY
      INTEGER NSOLD
      INTEGER NMONTH
C     
      INTEGER LAND_SEA(NX,NY)
      INTEGER NSOIL(NX,NY)
      INTEGER SOILTYP(NX,NY)
      INTEGER VEGTYP(NX,NY)
      INTEGER SLOPETYP(NX,NY)
      REAL    ELEV_2d(NX,NY)
C
      REAL    SOILDEPTH(NSOLD,NX,NY)
      REAL    ALBEDO(NMONTH,NX,NY)
      REAL    SHDFAC(NMONTH,NX,NY)
      REAL    TBOT(NX,NY)
      REAL    PRCP_MASK(NX,NY)
      REAL    PECLIM(NMONTH,NX,NY)
      REAL    PERATIO(NX,NY)
      REAL    PEADJ(NMONTH,NX,NY)

C     SAC PARAMETERS
      REAL  UZTWM_2d(NX,NY),UZFWM_2d(NX,NY)
      REAL  UZK_2d(NX,NY),PCTIM_2d(NX,NY)
      REAL  ADIMP_2d(NX,NY),RIVA_2d(NX,NY)
      REAL  ZPERC_2d(NX,NY),REXP_2d(NX,NY)
      REAL  LZTWM_2d(NX,NY),LZFSM_2d(NX,NY)
      REAL  LZFPM_2d(NX,NY),LZSK_2d(NX,NY)
      REAL  LZPK_2d(NX,NY),PFREE_2d(NX,NY)
      REAL  SIDE_2d(NX,NY),RSERV_2d(NX,NY)

C     SNOW17 PARAMETERS
      REAL SCF_2d(NX,NY),MFMAX_2d(NX,NY)
      REAL MFMIN_2d(NX,NY),UADJ_2d(NX,NY)
      REAL SI_2d(NX,NY),NMF_2d(NX,NY)
      REAL TIPM_2d(NX,NY),MBASE_2d(NX,NY)
      REAL PXTEMP_2d(NX,NY),PLWHC_2d(NX,NY)
      REAL DAYGM_2d(NX,NY)
      REAL ADC_2d(11,NX,NY)
C
C     SAC PARAMETERS
      REAL  UZTWM,UZFWM,UZK,PCTIM
      REAL  ADIMP,RIVA,ZPERC,REXP
      REAL  LZTWM,LZFSM,LZFPM,LZSK
      REAL  LZPK,PFREE,SIDE,RSERV

C     SNOW17 PARAMETERS
      REAL SCF,MFMAX,MFMIN,UADJ
      REAL SI,NMF,TIPM,MBASE
      REAL PXTEMP,PLWHC,DAYGM
      REAL ADC(11)
C
      INTEGER NREAD, I, J, K
      INTEGER ISNO

      WRITE(*,*) 'Read land surface characteristics'

      READ(65) ((LAND_SEA(J,K),J=1,NX),K=1,NY)
      READ(51) ((NSOIL(J,K),J=1,NX),K=1,NY)
      READ(53) ((SOILTYP(J,K),J=1,NX),K=1,NY)
      READ(54) ((VEGTYP(J,K),J=1,NX),K=1,NY)
      READ(55) ((SLOPETYP(J,K),J=1,NX),K=1,NY)
      READ(57) ((TBOT(J,K),J=1,NX),K=1,NY)
      READ(58) (((ALBEDO(I,J,K),I=1,NMONTH),J=1,NX),K=1,NY)
      READ(59) (((SHDFAC(I,J,K),I=1,NMONTH),J=1,NX),K=1,NY)
      READ(60) ((PRCP_MASK(J,K),J=1,NX),K=1,NY)
      READ(61) (((PECLIM(I,J,K),I=1,NMONTH),J=1,NX),K=1,NY)      
      READ(62) ((ELEV_2d(J,K),J=1,NX),K=1,NY)
      READ(63) ((PERATIO(J,K),J=1,NX),K=1,NY)
      READ(64) (((PEADJ(I,J,K),I=1,NMONTH),J=1,NX),K=1,NY)

      WRITE(*,*) 'Read sac parameters'

      READ(81) ((ADIMP_2d(J,K),J=1,NX),K=1,NY)
      READ(82) ((LZFPM_2d(J,K),J=1,NX),K=1,NY)
      READ(83) ((LZFSM_2d(J,K),J=1,NX),K=1,NY)
      READ(84) ((LZPK_2d(J,K),J=1,NX),K=1,NY) 
      READ(85) ((LZSK_2d(J,K),J=1,NX),K=1,NY)
      READ(86) ((LZTWM_2d(J,K),J=1,NX),K=1,NY)
      READ(87) ((PCTIM_2d(J,K),J=1,NX),K=1,NY)
      READ(88) ((PFREE_2d(J,K),J=1,NX),K=1,NY)
      READ(89) ((REXP_2d(J,K),J=1,NX),K=1,NY)
      READ(90) ((RIVA_2d(J,K),J=1,NX),K=1,NY) 
      READ(91) ((RSERV_2d(J,K),J=1,NX),K=1,NY)
      READ(92) ((SIDE_2d(J,K),J=1,NX),K=1,NY)
      READ(93) ((UZFWM_2d(J,K),J=1,NX),K=1,NY)
      READ(94) ((UZK_2d(J,K),J=1,NX),K=1,NY)
      READ(95) ((UZTWM_2d(J,K),J=1,NX),K=1,NY)
      READ(96) ((ZPERC_2d(J,K),J=1,NX),K=1,NY)


C     READ SNOW17 PARAMETERS

      ISNO = 71
C      OPEN (ISNO,FILE='SNOW_PAR',STATUS='OLD') 
      READ (ISNO,*)
      READ (ISNO,*) SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,PXTEMP,
     &              PLWHC,DAYGM
      READ (ISNO,*)
      READ (ISNO,*) (ADC(J),J=1,11)
      CLOSE(ISNO)

      DO J = 1, NY
         DO I = 1, NX
C     SNOW17 PARAMETERS
            SCF_2d(I,J)   = SCF
            MFMAX_2d(I,J) = MFMAX 
            MFMIN_2d(I,J) = MFMIN
            UADJ_2d(I,J)  = UADJ
            SI_2d(I,J)    = SI
            NMF_2d(I,J)   = NMF
            TIPM_2d(I,J)  = TIPM
            MBASE_2d(I,J) = MBASE
            PXTEMP_2d(I,J)= PXTEMP
            PLWHC_2d(I,J) = PLWHC
            DAYGM_2d(I,J) = DAYGM
            DO K = 1, 11
               ADC_2d(K,I,J)= ADC(K)
            END DO
         END DO
      END DO

      RETURN
      END
