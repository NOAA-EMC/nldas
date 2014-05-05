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
! cat_main.f:
!
! DESCRIPTION:
!  CATCHMENT MODEL: this is the main code for the CATCHMENT model - 
!  Significant F90 revisions will be required below this in the future.
!
!
! REVISION HISTORY:
!  4  Apr 2000: Jeffrey Walker; Initial code
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!                not a hard-wired undefined value of -999.
!   3 Aug 2000: Brian Cosgrove; fixed bug in wind vector calculation
!   3 Aug 2000: Brian Cosgrove; fixed bug in precip calculation
!  15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new 
!                GRIB libraries
!=========================================================================

      SUBROUTINE CAT_MAIN(LDAS,GRID,CAT)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      USE cat_module       ! Catchment variables
      IMPLICIT NONE
      type (ldasdec) LDAS
      type (griddec) GRID(LDAS%NC,LDAS%NR)
      type (catdec)  CAT(LDAS%NCATM)

!=== Local Variables =====================================================
      INTEGER :: SEC
      REAL, pointer :: SUNANG(:),RA(:)
      REAL, pointer :: FORC(:,:)
      INTEGER :: I,J
   
      LOGICAL*1 :: BUG=.false.,BUG1=.false.,BUG2=.false.
 
!=== End Variable Definition =============================================

!=== Transform forcing data from grid to catchment
! - 2m Air Temperature
      ALLOCATE (FORC(LDAS%NC,LDAS%NR))
      DO I=1,LDAS%NC
       DO J=1,LDAS%NR
        FORC(I,J)=GRID(I,J)%FORCING(1)
       END DO
      END DO
      CALL TRNSFM(1,LDAS,FORC,CAT%TMP2M,CAT%NFF)

! - 2m Dew Point Temperature
      DO I=1,LDAS%NC
       DO J=1,LDAS%NR
        IF (GRID(I,J)%FORCING(1).EQ.LDAS%UDEF .OR.
     &      GRID(I,J)%FORCING(2).EQ.LDAS%UDEF) THEN
         FORC(I,J)=LDAS%UDEF
        ELSE
         IF (GRID(I,J)%FORCING(2) .LT. 0.00001) THEN
          GRID(I,J)%FORCING(2)=0.00001
         END IF
         FORC(I,J)=5420./LOG(2.53E8*0.622/(GRID(I,J)%FORCING(2)*
     &             ((0.622-0.756*GRID(I,J)%FORCING(2))/
     &             (0.622-1.378*GRID(I,J)%FORCING(2)))*
     &             GRID(I,J)%FORCING(7)/1000.))
        END IF
       END DO
      END DO
      CALL TRNSFM(1,LDAS,FORC,CAT%DEW2M,CAT%NFF)

! - Downward Short Wave Radiation
      DO I=1,LDAS%NC
       DO J=1,LDAS%NR
        FORC(I,J)=GRID(I,J)%FORCING(3)
       END DO
      END DO
      CALL TRNSFM(1,LDAS,FORC,CAT%SWDN,CAT%NFF)

! - Downward Long Wave Radiation
      DO I=1,LDAS%NC
       DO J=1,LDAS%NR
        FORC(I,J)=GRID(I,J)%FORCING(4)
       END DO
      END DO
      CALL TRNSFM(1,LDAS,FORC,CAT%LWDN,CAT%NFF)

! - Wind Speed
      DO I=1,LDAS%NC
       DO J=1,LDAS%NR
        IF (GRID(I,J)%FORCING(5).EQ.LDAS%UDEF .OR.
     &      GRID(I,J)%FORCING(6).EQ.LDAS%UDEF) THEN
         FORC(I,J)=LDAS%UDEF
        ELSE
         FORC(I,J)=SQRT((GRID(I,J)%FORCING(5)**2.)+
     &             (GRID(I,J)%FORCING(6)**2.))
         IF (FORC(I,J) .LT. 0.01) FORC(I,J)=0.01
        END IF
       END DO
      END DO
      CALL TRNSFM(1,LDAS,FORC,CAT%WND,CAT%NFF)

! - Surface Pressure
      DO I=1,LDAS%NC
       DO J=1,LDAS%NR
        FORC(I,J)=GRID(I,J)%FORCING(7)
       END DO
      END DO
      CALL TRNSFM(1,LDAS,FORC,CAT%SFCPRS,CAT%NFF)

! - Total Precipitation
      DO I=1,LDAS%NC
       DO J=1,LDAS%NR
        FORC(I,J)=GRID(I,J)%FORCING(8)*LDAS%TS
       END DO
      END DO
      CALL TRNSFM(1,LDAS,FORC,CAT%TPTP,CAT%NFF)

! - Convective Precipitation
      DO I=1,LDAS%NC
       DO J=1,LDAS%NR
        FORC(I,J)=GRID(I,J)%FORCING(9)*LDAS%TS
       END DO
      END DO
      CALL TRNSFM(1,LDAS,FORC,CAT%CPCP,CAT%NFF)
      DEALLOCATE (FORC)

!=== Compute number of seconds into day
      SEC=(LDAS%MN*60)+(LDAS%HR*60*60)  !total sec into the day

!=== Compute Zenith angle of sun, SUNANG
      ALLOCATE (SUNANG(LDAS%NCATM),RA(LDAS%NCATM))
      CALL ASTRO(
     i           LDAS%YR,LDAS%MO,LDAS%DOY,SEC,CAT%LATT,CAT%LONN,
     i           LDAS%NCATM,
     o           SUNANG,RA)
      DEALLOCATE (RA)

!=== Forecast prognostic states of catchments
      CALL PROCESS(
     i             LDAS%NCATM,FLOAT(LDAS%TS),1.,LDAS%DOY,SUNANG,
     i             CAT%LONN,CAT%LATT,CAT%VEGCLS,CAT%TMP2M,CAT%DEW2M,
     i             CAT%SFCPRS,CAT%CPCP,CAT%TPTP,CAT%LWDN,CAT%SWDN,
     i             CAT%WND,CAT%BF(1),CAT%BF(2),CAT%BF(3),CAT%VGWMAX,
     i             CAT%CDCR(1),CAT%CDCR(2),CAT%PSIS,CAT%BEE,CAT%POROS,
     i             CAT%WPWET,CAT%COND,CAT%GNU,CAT%ARS(1),CAT%ARS(2),
     i             CAT%ARS(3),CAT%ARA(1),CAT%ARA(2),CAT%ARA(3),
     i             CAT%ARA(4),CAT%ARW(1),CAT%ARW(2),CAT%ARW(3),
     i             CAT%ARW(4),CAT%TSA(1),CAT%TSA(2),CAT%TSB(1),
     i             CAT%TSB(2),CAT%ALBEDO(LDAS%MO),CAT%GREEN(LDAS%MO),
     i             CAT%LAI(LDAS%MO),CAT%ZOL(LDAS%MO),CAT%VGD(LDAS%MO),
     i             BUG,BUG1,BUG2,
     u             CAT%TC(1),CAT%TC(2),CAT%TC(3),CAT%QA(1),CAT%QA(2),
     u             CAT%QA(3),CAT%CATDEF,CAT%RZEXC,CAT%SRFEXC,CAT%INT,
     u             CAT%GHT(1),CAT%GHT(2),CAT%GHT(3),CAT%GHT(4),
     u             CAT%GHT(5),CAT%GHT(6),CAT%TSURF,CAT%WESN(1),
     u             CAT%WESN(2),CAT%WESN(3),CAT%HTSN(1),CAT%HTSN(2),
     u             CAT%HTSN(3),CAT%SNDZ(1),CAT%SNDZ(2),CAT%SNDZ(3),
     o             CAT%EVAP,CAT%SHFLX,CAT%RUNOFF,CAT%EINT,CAT%ESOI,
     o             CAT%EVEG,CAT%ESNO,CAT%SNOWTERM,CAT%LWUP,CAT%SMELT,
     o             CAT%BFLOW,CAT%RUNSRF,CAT%AR(1),CAT%AR(2),CAT%RZEQ,
     o             CAT%INFIL,CAT%GHFLX,CAT%TSNOW,CAT%ASNOW,CAT%LHFLX,
     o             CAT%MODALB,CAT%TS(1),CAT%TS(2),CAT%TS(3),CAT%TS(4),
     o             CAT%TS(5),CAT%TS(6),CAT%ASSV(1),CAT%ASSV(2),
     o             CAT%ASSV(3),CAT%ASSV(4),CAT%ASSV(5),CAT%ASSV(6),
     o             CAT%ASSV(7),CAT%ASSV(8),CAT%ASSV(9),CAT%ASSV(10),
     o             CAT%ASSV(11))

      CAT%AR(3)=1.-(CAT%AR(1)+CAT%AR(2))

      DEALLOCATE (SUNANG)

!=== Diagnose soil moisture 
      CALL CALC_MOIST(
     i                LDAS%NCATM,CAT%SRFEXC,CAT%RZEXC,CAT%CATDEF,
     i                CAT%CDCR(1),CAT%CDCR(2),CAT%WPWET,CAT%VGWMAX,
     i                CAT%PSIS,CAT%BEE,CAT%POROS,CAT%ZDPTH(1),
     i                CAT%ZDPTH(2),CAT%ZDPTH(3),CAT%ARA(1),CAT%ARA(2),
     i                CAT%ARA(3),CAT%ARA(4),CAT%ARW(1),CAT%ARW(2),
     i                CAT%ARW(3),CAT%ARW(4),
     o                CAT%SRFMC,CAT%RZMC,CAT%COLMC)

      RETURN
      END

