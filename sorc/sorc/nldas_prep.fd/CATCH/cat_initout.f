!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! cat_initout.f: 
!
! DESCRIPTION:
!  Initialize CATCHMENT output arrays for use in LDAS
!
! REVISION HISTORY:
! 4 Apr 2000: Jeffrey Walker; Initial Code
! 8 Aug 2000: Brian Cosgrove; Fixed error in initialization of 
!             infil and bflow so that sumbflow and suminfil are 
!             now initialized instead of bflow and infil
!=========================================================================

       SUBROUTINE cat_initout (ldas,cat)

! Declare modules and data structures
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE cat_module
      IMPLICIT NONE
      type (ldasdec) LDAS
      type (catdec)  CAT(LDAS%NCATM)

!=== Local variables =====================================================

!=== End Variable List ===================================================

      CAT%SUMTP=0.
      CAT%SUMSNOW=0.
      CAT%SUMQA=0.
      CAT%SUMTC=0.
      CAT%SUMTS6=0.
      CAT%SUMSWUP=0.
      CAT%SUMLWUP=0.
      CAT%SUMSHFLX=0.
      CAT%SUMGHFLX=0.
      CAT%SUMLHFLX=0.
      CAT%SUMCATDEF=0.
      CAT%SUMRZEXC=0.
      CAT%SUMSRFEXC=0.
      CAT%SUMSRFMC=0.
      CAT%SUMRZMC=0.
      CAT%SUMCOLMC=0.
      CAT%SUMAR1=0.
      CAT%SUMAR3=0.
      CAT%SUMINT=0.
      CAT%SUMEINT=0.
      CAT%SUMEVAP=0.
      CAT%SUMEVEG=0.
      CAT%SUMESOI=0.
      CAT%SUMSNDZ=0.
      CAT%SUMSMELT=0.
      CAT%SUMINFIL=0.
      CAT%SUMBFLOW=0.
      CAT%SUMRUNSRF=0.
      CAT%SUMRUNOFF=0.

      CAT%COUNT=0

      END SUBROUTINE cat_initout




