 	subroutine griboutclm2 (ldas,grid,tile,fbase,fyrmodir)

***********************************************************************
C  PROGRAM:  GRIBLDAS           WRITE LDAS OUTPUT IN GRIB FORMAT
C  PRGMMR:  MARSHALL & LOHMANN  ORG:  W/NP20
C
C  ABSTRACT:  TAKES HR1LY OUTPUT FILES FROM LDAS MODELS (IN BINARY
C             FORMAT - BIG ENDIAN) AND MAKES GRIB RECORDS FOR ALL FIELDS
C             IN THE FILE.  NOTE ORDER OF FIELDS IN EACH BINARY FILE
C             MUST CORRESPOND TO ORDER OF FIELDS IN KPDS.tbl.  TIME
C             STAMP INFORMATION MUST ALSO APPEAR (YYYYMMDDHH) IN THE
C             TITLE OF EACH BINARY HR1LY OUTPUT FILE.  THIS INFO IS
C             USED IN CONSTRUCTING THE TIME-RELATED PDS OCTETS WHICH
C             DO NOT APPEAR IN KPDS.tbl   
C
C  PROGRAM HISTORY LOG:
C  00-10-15 MARSHALL ORIGINAL TEMPLATE FOR NOAH MODEL
C  00-11-01 LOHMANN  FINALIZED AND CHANGED TO SUBROUTINE
c  00-11-16 Cosgrove Altered for use with NASA LDAS model driver
c  00-11-17 Cosgrove Added ability to output on time intervals other
c           than hourly.  ALSO, ADDED PRINT STATEMENT IN BEGINNING
C           TO PRINT BLANK LINE...THIS KEEPS LDAS DRIVER FROM 
C           CRASHING IN THIS SUBROUTINE....THERE APPARENTLY IS A
C           MEMORY-RELATED BUG IN THE DRIVER, BUT EXTENSIVE DEBUGGING
C           FAILED TO FIND IT AND THE PRINT STATMENT WAS ADDED AS A 
C           STOP-GAP MEASURE TO ALLOW THE LDAS DRIVER TO RUN
C  01-01-19 Cosgrove Added CLOSE statement
C  01-04-05 Cosgrove Changed definition of KPDS 17 (number of timesteps
C           in accumulation or average)...made it variable instead of 
C           hardwired 
C  01-09-04 Cosgrove Updated--Changed output directory structure, enabled
C           statistical output, changed real to real*8 for tick call.  
C           Volumetric soil moisture output is commented out right now,
C           but could be enabled to output 5cm,25cm,60cm,75cm,tot-col,
C           layer1,layer2 and layer3 vol. soil mst.
C  02-02-05 Jambor Removed misleading code regarding KPDS(13) values
C           that had been commented out.
C  02-03-19 Cosgrove--Based heavily on gribout.f, created griboutclm.f
C           This subroutine outputs LDAS list of CLM variables in GRIB.
C	    Presently, root zone soil moisture is calculated to include
C	    only those layers in which 95% of roots lie
C  02-03-28 Cosgrove--Fixed code so that MSTAV moisture availability
C           variable is correct
C  02-06-30 Gottschalck--Modified from griboutclm.f to work with CLM2
C***********************************************************************

! Declare modules and data structures
      use ldas_module      ! ldas non-model-specific 1-d variables
      use tile_module      ! ldas non-model-specific tile variables
      use grid_module      ! ldas non-model-specific grid variables
      use clmtype
      use clm_varder          ! 1-D CLM variables
      use clm_varcon, ONLY : denh2o, denice, hvap, hsub, hfus, istwet
      use clm_varpar, ONLY : nlevsoi, nlevsno
      use clm_varmap, ONLY : patchvec
      implicit none
      type (ldasdec) ldas
      type (tiledec)      ::  tile(ldas%nch)
      type (griddec)      ::  grid(ldas%nc,ldas%nr)

      INTEGER c,r,t,m,cc,jj,counter,p,pp
      INTEGER SS1,TS,MN1,HR1,DA1,MO1,YR1,TS1,DOY1
      INTEGER NSOLD,flag
      PARAMETER (NSOLD=13)
      CHARACTER*8 TODAY, YESTERDAY
      CHARACTER*1 TOD(8),YES(8)
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13),FTIMED(4)
      CHARACTER*1  FYRMODIR(18),FSUBFT(80),EXPCODE(3),FTIMEC(8)
      CHARACTER*1  FSUBFG(80),EXPARRAY(80),FTIMEB(10),FSUBGB(11)
      CHARACTER*1 GRIBF(256),GRIBF2(256)

      INTEGER LUGB, I, J, K, N, IRET, KPDS(25), KGDS(22),
     &        NLDAS, LENGDS, IG, JRET, NFIELDS

      LOGICAL*1  LDASMASK(LDAS%NC,LDAS%NR)

      CHARACTER CENT*2, CHOUR*2, GDS(400)
      CHARACTER*256 GRIBFILE

      REAL    VMEAN,VSTDEV,VMIN,VMAX

      REAL    NSWRS(LDAS%NC,LDAS%NR) 
      REAL    NLWRS(LDAS%NC,LDAS%NR)
      REAL    LHTFL(LDAS%NC,LDAS%NR) 
      REAL    SHTFL(LDAS%NC,LDAS%NR) 
      REAL    GFLUX(LDAS%NC,LDAS%NR) 
      REAL    DSWRF(LDAS%NC,LDAS%NR)
      REAL    DLWRF(LDAS%NC,LDAS%NR)
      REAL    ARAIN(LDAS%NC,LDAS%NR)
      REAL    EVP(LDAS%NC,LDAS%NR)
      REAL    SSRUN(LDAS%NC,LDAS%NR)
      REAL    BGRUN(LDAS%NC,LDAS%NR) 
      REAL    AVSFT(LDAS%NC,LDAS%NR)
      REAL    ALBDO(LDAS%NC,LDAS%NR) 
      REAL    WEASD(LDAS%NC,LDAS%NR)
      REAL    CWAT(LDAS%NC,LDAS%NR)
      REAL    SOILT(nsold,LDAS%NC,LDAS%NR)
      REAL    SOILM(nsold,LDAS%NC,LDAS%NR)
      REAL    SOILL(nsold,LDAS%NC,LDAS%NR)
      REAL    SOILV(nsold,LDAS%NC,LDAS%NR)
      REAL    MSTAV(nsold,LDAS%NC,LDAS%NR)
      REAL    PEVPR(LDAS%NC,LDAS%NR)
      REAL    VEG(LDAS%NC,LDAS%NR)
      REAL    LAI(LDAS%NC,LDAS%NR)
      REAL    SNOD(LDAS%NC,LDAS%NR)
      REAL    SNOHF(LDAS%NC,LDAS%NR)
      REAL    ASNOW(LDAS%NC,LDAS%NR)
      REAL    SNOM(LDAS%NC,LDAS%NR)
      REAL    TRANS(LDAS%NC,LDAS%NR)
      REAL    EVBS(LDAS%NC,LDAS%NR)
      REAL    EVCW(LDAS%NC,LDAS%NR)
      REAL    SBSNO(LDAS%NC,LDAS%NR)
      REAL    ACOND(LDAS%NC,LDAS%NR)
      REAL    CCOND(LDAS%NC,LDAS%NR)
      REAL    SNOC(LDAS%NC,LDAS%NR)
      REAL    SALBD(LDAS%NC,LDAS%NR)
      REAL    TMP(LDAS%NC,LDAS%NR)
      REAL    SPFH(LDAS%NC,LDAS%NR)
      REAL    UGRD(LDAS%NC,LDAS%NR)
      REAL    VGRD(LDAS%NC,LDAS%NR)
      REAL    PRES(LDAS%NC,LDAS%NR)
      REAL    ACPCP(LDAS%NC,LDAS%NR)
      REAL   SNOWT(LDAS%NC,LDAS%NR)
      REAL   VEGT(LDAS%NC,LDAS%NR)
      REAL   BARET(LDAS%NC,LDAS%NR)
      REAL   RADT(LDAS%NC,LDAS%NR)

      real   aaa,bbb,ccc
      real   snowtemp(ldas%nch),totaldepth(ldas%nch),avgsurft(ldas%nch)
      real   soilwtc(ldas%nch),soilwr(ldas%nch),rootfr
      real   swetwilt(ldas%nch),swetint(ldas%nch)
      real   totmst(ldas%nch),watsat(ldas%nch)
      real   avgwatsat(ldas%nch)

      real   snowtt(ldas%nch),soilmr(ldas%nch),soilm1m(ldas%nch)
      REAL   DUMMYGMT,DUMMY(LDAS%NC,LDAS%NR),depth,factor
      REAL*8 DUMMYTIME

!	Initialize variables to zero
	DO j=1,LDAS%NR
	DO i=1,LDAS%NC
         NSWRS(I,J)=0.0
         NLWRS(I,J)=0.0
         LHTFL(I,J)=0.0
         SHTFL(I,J)=0.0
         GFLUX(I,J)=0.0
         DSWRF(I,J)=0.0
         DLWRF(I,J)=0.0
         ARAIN(I,J)=0.0
         EVP(I,J)=0.0
         SSRUN(I,J)=0.0
         BGRUN(I,J)=0.0
         AVSFT(I,J)=0.0
         ALBDO(I,J)=0.0
         WEASD(I,J)=0.0
         CWAT(I,J)=0.0
	 DO k=1,nsold
          SOILV(k,i,j)=0.0
          SOILT(k,i,j)=0.0
          SOILM(k,i,j)=0.0
          MSTAV(k,i,j)=0.0
	  SOILL(k,i,j)=0.0
	 ENDDO
         PEVPR(I,J)=0.0
         VEG(I,J)=0.0
         LAI(I,J)=0.0
         SNOD(I,J)=0.0
         SNOHF(I,J)=0.0
         ASNOW(I,J)=0.0
         SNOM(I,J)=0.0
         DUMMY(I,J)=0.0
         EVCW(I,J)=0.0
         TRANS(I,J)=0.0
         EVBS(I,J)=0.0
         SBSNO(I,J)=0.0
         ACOND(I,J)=0.0
         CCOND(I,J)=0.0
         SNOC(I,J)=0.0
         SALBD(I,J)=0.0
         TMP(I,J)=0.0
         SPFH(I,J)=0.0
         UGRD(I,J)=0.0
         VGRD(I,J)=0.0
         PRES(I,J)=0.0
         ACPCP(I,J)=0.0
         SNOWT(I,J)=0.0
         VEGT(I,J)=0.0
         BARET(I,J)=0.0
         RADT(I,J)=0.0
	ENDDO
	ENDDO
!	Construct GRIB logical type land/sea mask from LDAS integer land/sea mask
      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
          IF(GRID(C,R)%IMASK.EQ.1) THEN
            LDASMASK(C,R)=.TRUE.
          ELSE
            LDASMASK(C,R)=.FALSE.
          ENDIF
       ENDDO
      ENDDO
        HR1=LDAS%HR
        DA1=LDAS%DA
        MO1=LDAS%MO
        YR1=LDAS%YR
        MN1=LDAS%MN
        SS1=0
	TS1=-3600*24
	DUMMYGMT=1.0
        DUMMYTIME=1.0

!	Construct Grib output file names
109    FORMAT(I4,I2,I2)
110    FORMAT(8A1)
111    FORMAT(I4,I2,I2,I2)
112    FORMAT(10A1)
113    FORMAT(I4)
114    FORMAT(4A1)

       OPEN(95,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
        WRITE(95,109,REC=1)YR1,MO1,DA1
        READ(95,110,REC=1)TOD
        DO I=1,8
         IF(TOD(I).EQ.(' '))TOD(I)='0'
        ENDDO
	TODAY=TOD(1)//TOD(2)//TOD(3)//TOD(4)//TOD(5)
     &  //TOD(6)//TOD(7)//TOD(8)

      CALL TICK(DUMMYTIME,DOY1,DUMMYGMT,YR1,MO1,DA1,HR1,MN1,SS1,TS1)

        WRITE(95,109,REC=1)YR1,MO1,DA1
        READ(95,110,REC=1)YES
        DO I=1,8
         IF(YES(I).EQ.(' '))YES(I)='0'
        ENDDO
        YESTERDAY=YES(1)//YES(2)//YES(3)//YES(4)//YES(5)
     &  //YES(6)//YES(7)//YES(8)
        DO I=1,25
         KPDS(I)=0
        ENDDO
        DO I=1,22
         KGDS(I)=0
        ENDDO
C     NON-CHANGING PDS ELEMENTS
      KPDS(1)=222               !ID FOR GSFC PRODUCTS
      KPDS(2)=222               !ID FOR Mosaic MODEL (change value for other models)
      KPDS(4)=192               !BMS FLAG... DOn't worry about this.
      KPDS(12)=0                !assume output time minute always = 0
      KPDS(13)=1                !Forecast Time Unit (Hours)
      KPDS(17)=INT((LDAS%writeintc2*3600.0)/LDAS%TS) !number of time steps in
                                                    !averaged/accum variables
      KPDS(18)=0                !GRIB version -- left as 0 in NCEP products
      KPDS(19)=1                !version number of KPDS.tbl for LDAS.  
      KPDS(20)=0                !none missing from averages/accumulations (always4)
      KPDS(23)=222              !GSFC ID#
      KPDS(24)=0                !DOes not apply to LDAS output
      KPDS(25)=0                !not used

C          IF(GRID(C,R)%IMASK.EQ.0)Gtmp(C,R)=LDAS%UDEF  !Set Water to undefined
       DO T=1,LDAS%NCH
!	Time Averaged
        NSWRS(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   NSWRS(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%totfsa*
     1    patchvec%wtxy(t)
        NLWRS(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 NLWRS(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%toteflx_lwrad_net*
     1 patchvec%wtxy(t)
        LHTFL(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 LHTFL(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%toteflx_lh_tot*
     1 patchvec%wtxy(t)
        SHTFL(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 SHTFL(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%toteflx_sh_tot*
     1 patchvec%wtxy(t)


        GFLUX(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 GFLUX(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%toteflx_soil_grnd*
     1 patchvec%wtxy(t)


        SNOHF(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 SNOHF(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%totqflx_snomelt*
     1     (-382000.)*patchvec%wtxy(t)/FLOAT(clm(T)%COUNT)
        DSWRF(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 DSWRF(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%totsolisbd*
     1    patchvec%wtxy(t)
        DLWRF(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 DLWRF(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%totforc_lwrad*
     1    patchvec%wtxy(t)
!	Accumulated
        ASNOW(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      ASNOW(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%totsnow*patchvec%wtxy(t)
        ARAIN(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      ARAIN(patchvec%ixy(t),patchvec%jxy(t))+
     1  clm(T)%totrain*patchvec%wtxy(t)
        EVP(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      EVP(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%totqflx_evap*patchvec%wtxy(t)
        SSRUN(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SSRUN(patchvec%ixy(t),patchvec%jxy(t))+(clm(T)%totqflx_surf*
     1         *patchvec%wtxy(t))
        BGRUN(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   BGRUN(patchvec%ixy(t),patchvec%jxy(t))+(clm(T)%totqflx_drain*
     1         patchvec%wtxy(t))
        SNOM(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SNOM(patchvec%ixy(t),patchvec%jxy(t))+(clm(T)%totqflx_snomelt
     1        *patchvec%wtxy(t))
!	Instantaneous

!Snow Temperature Calculation
        snowtemp(t)=0.
         if (clm(t)%itypwat/=istwet)then
          if(clm(t)%snl < 0)then
           totaldepth(t)=0.
           do i=clm(t)%snl+1,0    
! Compute total depth of snow layers
            totaldepth(t)=totaldepth(t)+clm(t)%dz(i)
           enddo

           do i=clm(t)%snl+1,0    
! Compute snow temperature
            snowtemp(t)=snowtemp(t)+(clm(t)%t_soisno(i)*clm(t)%dz(i))
           enddo
           snowtemp(t)=snowtemp(t)/totaldepth(t)
         endif
         if(snowtemp(t).eq.0)snowtemp(t)=ldas%udef
        endif

        SNOWT(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      SNOWT(patchvec%ixy(t),patchvec%jxy(t))+
     1      snowtemp(t)*patchvec%wtxy(t)

        VEGT(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      VEGT(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%t_veg*patchvec%wtxy(t)



        BARET(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      BARET(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%t_grnd*patchvec%wtxy(t)


! SnowTT is the snow surface temperature, i.e. top layer t_soisno


         snowtt(t)=0.
         if (clm(t)%itypwat/=istwet)then
         if(clm(t)%snl < 0)then
           snowtt(t)=clm(t)%t_soisno(clm(t)%snl+1)
         endif
         endif
         if(snowtt(t)==0.)snowtt(t)=-9999.0  !SnowT is undefined when there is no snow

        if(snowtt(t).ne.-9999.0)then
        avgsurft(t)=clm(t)%frac_sno*snowtt(t)+
     &  clm(t)%frac_veg_nosno*clm(t)%t_veg+  
     & (1-(clm(t)%frac_sno+clm(t)%frac_veg_nosno))*clm(t)%t_grnd
        else
        avgsurft(t)=clm(t)%frac_veg_nosno*clm(t)%t_veg+ 
     &         (1-clm(t)%frac_veg_nosno)*clm(t)%t_grnd
        endif
        AVSFT(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      AVSFT(patchvec%ixy(t),patchvec%jxy(t))+
     1    avgsurft(t)*patchvec%wtxy(t)

        RADT(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      RADT(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%t_rad*patchvec%wtxy(t)

        ALBDO(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      ALBDO(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%surfalb_old*patchvec%wtxy(t)*
     1                  100.0
        WEASD(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      WEASD(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%h2osno*patchvec%wtxy(t)

        CWAT(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      CWAT(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%h2ocan*patchvec%wtxy(t)

        SOILT(1,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(1,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(1)*patchvec%wtxy(t)
        SOILT(2,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(2,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(2)*patchvec%wtxy(t)
        SOILT(3,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(3,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(3)*patchvec%wtxy(t)
        SOILT(4,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(4,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(4)*patchvec%wtxy(t)
        SOILT(5,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(5,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(5)*patchvec%wtxy(t)
        SOILT(6,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(6,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(6)*patchvec%wtxy(t)
        SOILT(7,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(7,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(7)*patchvec%wtxy(t)
        SOILT(8,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(8,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(8)*patchvec%wtxy(t)
        SOILT(9,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(9,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(9)*patchvec%wtxy(t)
        SOILT(10,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILT(10,patchvec%ixy(t),patchvec%jxy(t))+
     1    clm(T)%t_soisno(10)*patchvec%wtxy(t)



        SOILM(4,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(4,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(1)+
     1       clm(T)%h2osoi_ice(1))*patchvec%wtxy(t)
        SOILM(5,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(5,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(2)+
     1       clm(T)%h2osoi_ice(2))*patchvec%wtxy(t)
        SOILM(6,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(6,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(3)+
     1       clm(T)%h2osoi_ice(3))*patchvec%wtxy(t)
        SOILM(7,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(7,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(4)+
     1       clm(T)%h2osoi_ice(4))*patchvec%wtxy(t)
        SOILM(8,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(8,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(5)+
     1       clm(T)%h2osoi_ice(5))*patchvec%wtxy(t)
        SOILM(9,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(9,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(6)+
     1       clm(T)%h2osoi_ice(6))*patchvec%wtxy(t)
        SOILM(10,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(10,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(7)+
     1       clm(T)%h2osoi_ice(7))*patchvec%wtxy(t)
        SOILM(11,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(11,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(8)+
     1       clm(T)%h2osoi_ice(8))*patchvec%wtxy(t)

        SOILM(12,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(12,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(9)+
     1       clm(T)%h2osoi_ice(9))*patchvec%wtxy(t)
        SOILM(13,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(13,patchvec%ixy(t),patchvec%jxy(t))+
     1       (clm(T)%h2osoi_liq(10)+
     1       clm(T)%h2osoi_ice(10))*patchvec%wtxy(t)

        SOILM(1,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(1,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(4,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(5,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(6,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(7,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(8,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(9,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(10,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(11,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(12,patchvec%ixy(t),patchvec%jxy(t))+
     1   SOILM(13,patchvec%ixy(t),patchvec%jxy(t))


	rootfr=0.0
        soilmr(t)=0.0
!	Add up soil moisture (mm) in root zone
        do m=1,nlevsoi
          if (rootfr.le.0.95) then
c          soilmr(t)=soilmr(t)+clm(t)%rootfr(m)*(clm(t)%h2osoi_liq(m)+
c     &     clm(t)%h2osoi_ice(m))
          soilmr(t)=soilmr(t)+(clm(t)%h2osoi_liq(m)+
     &     clm(t)%h2osoi_ice(m))

          rootfr=rootfr+clm(t)%rootfr(m)

          endif
        enddo
c          soilmr(t)=soilmr(t)*(1.0/rootfr)

        SOILM(2,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(2,patchvec%ixy(t),patchvec%jxy(t))+
     1        soilmr(t)*patchvec%wtxy(t)


        depth=0.0
        flag=0
        do i=1,nlevsoi
          if (flag.ne.1) then
            soilm1m(t)=0.0
            depth=depth+clm(t)%dz(i)
            if (depth.ge.1.0) then
              factor=(1.0-((depth-1.0)/clm(t)%dz(i)))
              if (depth.le.1.00001) then
                do j=1,i
                  soilm1m(t)=soilm1m(t)+(clm(T)%h2osoi_liq(j)+
     &                   clm(T)%h2osoi_ice(j))
                enddo
                flag=1
              else
                do j=1,i-1
                  soilm1m(t)=soilm1m(t)+(clm(T)%h2osoi_liq(j)+
     &                       clm(T)%h2osoi_ice(j))
                enddo
                soilm1m(t)=soilm1m(t)+(factor*(clm(T)%h2osoi_liq(i)+
     &                     clm(T)%h2osoi_ice(i)))
                flag=1
              endif
            endif
          endif
        enddo
        if (flag.eq.0) then
          do i=1,nlevsoi
            soilm1m(t)=soilm1m(t)+(clm(T)%h2osoi_liq(i)+
     &  clm(T)%h2osoi_ice(i))
          enddo
          soilm1m(t)=soilm1m(t)/depth
        endif

        SOILM(3,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILM(3,patchvec%ixy(t),patchvec%jxy(t))+
     1        (soilm1m(t))*patchvec%wtxy(t)


        SOILL(1,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(1,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(1)*patchvec%wtxy(t)
        SOILL(2,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(2,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(2)*patchvec%wtxy(t)
        SOILL(3,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(3,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(3)*patchvec%wtxy(t)
        SOILL(4,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(4,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(4)*patchvec%wtxy(t)
        SOILL(5,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(5,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(5)*patchvec%wtxy(t)
        SOILL(6,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(6,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(6)*patchvec%wtxy(t)
        SOILL(7,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(7,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(7)*patchvec%wtxy(t)
        SOILL(8,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(8,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(8)*patchvec%wtxy(t)
        SOILL(9,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(9,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(9)*patchvec%wtxy(t)
        SOILL(10,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SOILL(10,patchvec%ixy(t),patchvec%jxy(t))+
     1        clm(T)%h2osoi_liq(10)*patchvec%wtxy(t)


        swetint(t)=0.
        avgwatsat(t)=0.
	soilwtc=0.

        do m=1,nlevsoi
        avgwatsat(t)=avgwatsat(t)+clm(t)%dz(m)*1000.0*clm(t)%watsat(m)
         swetint(t)=swetint(t)+clm(t)%h2osoi_liq(m)+
     &        clm(t)%h2osoi_ice(m)
        enddo
        soilwtc(t)=100.*swetint(t)/avgwatsat(t)

        MSTAV(1,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   MSTAV(1,patchvec%ixy(t),patchvec%jxy(t))+soilwtc(t)*
     1      patchvec%wtxy(t)

        swetwilt(t)=0.
        totaldepth(t)=0.
        avgwatsat(t)=0.
	soilwr(t)=0.
	rootfr=0.0

        soilwr(t)=0.
        rootfr=0.0
	totaldepth(t)=0.0
        do m=1,nlevsoi
          if (rootfr.le.0.95) then
           totaldepth(t)=totaldepth(t)+clm(t)%dz(m)
           rootfr=rootfr+clm(t)%rootfr(m)
	 endif
	enddo

        rootfr=0.0
	watsat(t)=0.0
	totmst(t)=0.0
        do m=1,nlevsoi
          if (rootfr.le.0.95) then
           swetwilt(t)=swetwilt(t) + (clm(t)%dz(m)*(clm(t)%watsat(m)*
     &    ((-1)*clm(t)%smpmax/clm(t)%sucsat(m))**(-1/clm(t)%bsw(m))))
     &     *(clm(t)%dz(m)/totaldepth(t))/clm(t)%dz(m)
           watsat(t)=watsat(t)+
     &     (clm(t)%watsat(m)*clm(t)%dz(m)/totaldepth(t))
           totmst(t)=totmst(t)+clm(t)%h2osoi_liq(m)+
     &     clm(t)%h2osoi_ice(m)
           rootfr=rootfr+clm(t)%rootfr(m)
          endif
        enddo

        aaa=((totmst(t))/(totaldepth(t)*  
     &       1000.0*watsat(t)))*watsat(t)
        bbb=(watsat(t))*swetwilt(t)
        ccc=watsat(t)-(watsat(t)*swetwilt(t))

        soilwr(t)=100.0*((aaa-bbb)/ccc)

        MSTAV(2,patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   MSTAV(2,patchvec%ixy(t),patchvec%jxy(t))+soilwr(t)*
     1    patchvec%wtxy(t)

!	Time Averaged
        EVCW(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   EVCW(patchvec%ixy(t),patchvec%jxy(t))+clm(t)%totqflx_ecanop*
     2    patchvec%wtxy(t)
        TRANS(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 TRANS(patchvec%ixy(t),patchvec%jxy(t))+clm(t)%totqflx_tran_veg*
     2 patchvec%wtxy(t)/FLOAT(clm(T)%COUNT)
        EVBS(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 EVBS(patchvec%ixy(t),patchvec%jxy(t))+clm(t)%totqflx_evap_grnd*
     2 patchvec%wtxy(t)/FLOAT(clm(T)%COUNT)
        SBSNO(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 SBSNO(patchvec%ixy(t),patchvec%jxy(t))+clm(t)%totqflx_sub_snow*
     2 patchvec%wtxy(t)/FLOAT(clm(T)%COUNT)
	 
!	Instantaneous
        ACOND(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      ACOND(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%acond*patchvec%wtxy(t)
       CCOND(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 ldas%udef
        VEG(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      VEG(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(T)%acond*patchvec%wtxy(t)

!        LAI(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
!     1   LAI(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%tlai*TILE(T)%FGRD
        LAI(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 LAI(patchvec%ixy(t),patchvec%jxy(t))+clm(t)%tlai*patchvec%wtxy(t)
       SNOD(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1     SNOD(patchvec%ixy(t),patchvec%jxy(t))+
     1     clm(T)%snowdp*patchvec%wtxy(t)
       SNOC(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1     SNOC(patchvec%ixy(t),patchvec%jxy(t))+
     1     clm(T)%frac_sno*patchvec%wtxy(t)*
     1                100.0
        SALBD(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1   SALBD(patchvec%ixy(t),patchvec%jxy(t))+clm(T)%snoalb_old*100.0*
     1      patchvec%wtxy(t)
        TMP(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      TMP(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(t)%forc_t*patchvec%wtxy(t)
        SPFH(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      SPFH(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(t)%forc_q*patchvec%wtxy(t)
        UGRD(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      UGRD(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(t)%forc_u*patchvec%wtxy(t)
        VGRD(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      VGRD(patchvec%ixy(t),patchvec%jxy(t))+
     1      clm(t)%forc_v*patchvec%wtxy(t)
        PRES(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1 PRES(patchvec%ixy(t),patchvec%jxy(t))+(clm(t)%forc_pbot/100.0)*
     1 patchvec%wtxy(t)     
! 	Accumulation
        ACPCP(patchvec%ixy(t),patchvec%jxy(t)) =   !Transfer Tile to Grid
     1      ACPCP(patchvec%ixy(t),patchvec%jxy(t))+
     1      (tile(T)%forcing(9)*ldas%ts)*
     1  patchvec%wtxy(t)

       ENDDO

C     DEFINE THE GDS FOR THE LDAS GRID (THIS ARRAY NEVER CHANGES FROM
C     FIELD TO FIELD. ONLY HAS TO BE DEFINED ONCE, UNLIKE THE PDS).

      IF (LDAS%DOMAIN==1) THEN
        IG = 236    !FOR LDAS GRID IN SUBROUTINE w3fi71.f in w3lib

        CALL MAKGDS(IG, KGDS, GDS, LENGDS, IRET)
      ELSE
        DO J=1,22
          KGDS(J)=LDAS%LDAS_KGDS(J)
        ENDDO
      ENDIF

C     read past header info in KPDS.tbl
      OPEN (UNIT = 30, FILE = './SRC/KPDS_completeclm2.tbl')
      DO K = 1, 42
         READ(30,*)
      END DO

C     MAKE OUTPUT FILENAME
 92    FORMAT(80A1)
 93    FORMAT(A256)
 91    FORMAT(256A1)
 94    FORMAT(A80)
 95    format(80A,A4,3A1,A6,4A1,A1,8A1,A1,10A1,11A1)
101    FORMAT(A9)
102    FORMAT (40I)
103    FORMAT (3I)


        WRITE(95,111,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
        READ(95,112,REC=1)FTIMEB
        DO I=1,10
         IF(FTIMEB(I).EQ.(' '))FTIMEB(I)='0'
        ENDDO

        WRITE(95,109,REC=1)LDAS%YR,LDAS%MO,LDAS%DA
        READ(95,110,REC=1)FTIMEC
        DO I=1,8
         IF(FTIMEC(I).EQ.(' '))FTIMEC(I)='0'
        ENDDO

        WRITE(95,113,REC=1)LDAS%YR
        READ(95,114,REC=1)FTIMED
        DO I=1,4
         IF(FTIMED(I).EQ.(' '))FTIMED(I)='0'
        ENDDO


        WRITE(95,101,REC=1)'.clm2.grb'
        READ(95,92,REC=1) (FSUBGB(I),I=1,9)
        C=0
       DO I=1,40
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO



C      IF ( (HR1.LE. 23).AND.(HR1.GE.13) ) THEN
         WRITE (CHOUR, '(i2.2)') HR1
        WRITE(95,103,REC=1)LDAS%EXPCODE
c	WRITE(95,102,REC=1)LDAS%EXPCODE
        READ(95,92,REC=1)EXPARRAY
	COUNTER=1
        DO I=1,40
         IF(EXPARRAY(I).NE.(' ')) THEN
          EXPCODE(COUNTER)=EXPARRAY(I)
          COUNTER=COUNTER+1
	 ENDIF
        ENDDO

c        WRITE(95,95,REC=1)(FBASE(I),I=1,C),FYRMODIR,'/',
c     &  'LDAS.E',EXPCODE,'.',
c     &   TODAY,CHOUR,'.MOSgrib'
        WRITE(95,95,REC=1)(FBASE(I),I=1,C),'/EXP',EXPCODE,'/CLM2/',
     &  FTIMED,'/',FTIMEC,'/',
     &   FTIMEB,(FSUBGB(I),I=1,9)
	READ(95,93,REC=1) GRIBFILE
c        print *,'gribfile=',gribfile,'end',LDAS%EXPCODE,'end'
        CLOSE (95)
C      ELSE 
C         WRITE (CHOUR, '(i2.2)') HR1
C        WRITE(95,95,REC=1)(FBASE(I),I=1,C),FYRMODIR,'/',
C     &   YESTERDAY,CHOUR,'_MOSAIC.LDASGRIB'
C        READ(95,93,REC=1) GRIBFILE
C      END IF

C     OPEN GRIB FILE

      LUGB = 40
      CALL BAOPEN (LUGB, GRIBFILE, IRET)
c	print *,'gribfile=',gribfile,'iret=',iret

C     WRITE EACH FIELD (IN GRIB) TO HR1LY OUTPUT FIELD
 15   FORMAT (29x, 7I6) 

      NLDAS = LDAS%NC * LDAS%NR

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
      IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
c	print *,'rian, hr1 in write=',hr1
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,NSWRS,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = NSWRS'
         STOP
      ELSE
C         PRINT*, 'wrote NSWRS'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
      IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,NLWRS,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = NLWRS'
         STOP
      ELSE
C         PRINT*, 'wrote NLWRS'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,LHTFL,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = LHTFL'
         STOP
      ELSE
C         PRINT*, 'wrote LHTFL'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SHTFL,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SHTFL'
         STOP
      ELSE
C         PRINT*, 'wrote SHTFL'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,GFLUX,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = GFLUX'
         STOP
      ELSE
C         PRINT*, 'wrote GFLUX'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOHF,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOHF'
         STOP
      ELSE
C         PRINT*, 'wrote SNOHF'
      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DSWRF,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DSWRF'
         STOP
      ELSE
C         PRINT*, 'wrote DSWRF'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DLWRF,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = DLWRF'
         STOP
      ELSE
C         PRINT*, 'wrote DLWRF'
      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ASNOW,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ASNOW'
         STOP
      ELSE
C         PRINT*, 'wrote ASNOW'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ARAIN,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ARAIN'
         STOP
      ELSE
C         PRINT*, 'wrote ARAIN'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,EVP,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = EVP'
         STOP
      ELSE
C         PRINT*, 'wrote EVP'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SSRUN,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SSRUN'
         STOP
      ELSE
C         PRINT*, 'wrote SSRUN'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,BGRUN,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = BGRUN'
         STOP
      ELSE
C         PRINT*, 'wrote BGRUN'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOM,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOM'
         STOP
      ELSE
C         PRINT*, 'wrote SNOM'
      END IF




      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOWT,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOWT'
         STOP
      ELSE
C         PRINT*, 'wrote SNOWT'
      END IF



      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,VEGT,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = VEGT'
         STOP
      ELSE
C         PRINT*, 'wrote VEGT'
      END IF



      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,BARET,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = BARET'
         STOP
      ELSE
C         PRINT*, 'wrote BARET'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,AVSFT,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = AVSFT'
         STOP
      ELSE
C         PRINT*, 'wrote AVSFT'
      END IF





      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)      
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,RADT,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = RADT'
         STOP
      ELSE
C         PRINT*, 'wrote RADT'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ALBDO,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ALBDO'
         STOP
      ELSE
C         PRINT*, 'wrote ALBDO'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,WEASD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = WEASD'
         STOP
      ELSE
C         PRINT*, 'wrote WEASD'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,CWAT,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = CWAT'
         STOP
      ELSE
C         PRINT*, 'wrote CWAT'
      END IF




      DO N = 1,10
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILT(N,I,J)
            END DO
         END DO
         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +        KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C        MAKE TIME DEPENDENT PDS PARAMETERS
         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = SOILT',N
            STOP
         ELSE
C         PRINT*, 'wrote SOILT ', N
         END IF
      END DO

c      DO N = 1, 8
c         DO J = 1,LDAS%NR
c            DO I = 1,LDAS%NC
c               DUMMY(I,J) = SOILV(N,I,J)
c            END DO
c         END DO
c         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
c     +        KPDS(15), KPDS(16), KPDS(22)
c        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2
c
C        MAKE TIME DEPENDENT PDS PARAMETERS
c         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
c         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
c         IF (IRET .NE. 0) THEN
c            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = SOILV',N
c            STOP
c         ELSE
C         PRINT*, 'wrote SOILV ',N
c         END IF
c      END DO




      DO N = 1,13
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILM(N,I,J)
            END DO
         END DO
         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +        KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C        MAKE TIME DEPENDENT PDS PARAMETERS
         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = SOILM',N
            STOP
         ELSE
C         PRINT*, 'wrote SOILM ',N
         END IF
      END DO

      DO N = 1,10
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = SOILL(N,I,J)
            END DO
         END DO
         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +        KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C        MAKE TIME DEPENDENT PDS PARAMETERS
         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = SOILL',N
            STOP
         ELSE
C         PRINT*, 'wrote SOILL ',N
         END IF
      END DO

      DO N = 1, 2
         DO J = 1,LDAS%NR
            DO I = 1,LDAS%NC
               DUMMY(I,J) = MSTAV(N,I,J)
            END DO
         END DO
         READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +        KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C        MAKE TIME DEPENDENT PDS PARAMETERS
         CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
         CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DUMMY,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HR1,',','FIELD = MSTAV',N
            STOP
         ELSE
C         PRINT*, 'wrote MSTAV ',N
         END IF
      END DO


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,EVCW,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = EVCW'
         STOP
      ELSE
C         PRINT*, 'wrote EVCW'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,TRANS,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = TRANS'
         STOP
      ELSE
C         PRINT*, 'wrote TRANS'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,EVBS,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = EVBS'
         STOP
      ELSE
C         PRINT*, 'wrote EVBS'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SBSNO,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SBSNO'
         STOP
      ELSE
C         PRINT*, 'wrote SBSNO'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ACOND,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ACOND'
         STOP
      ELSE
C         PRINT*, 'wrote ACOND'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,CCOND,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = CCOND'
         STOP
      ELSE
C         PRINT*, 'wrote CCOND'
      END IF


c      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
c     +              KPDS(15), KPDS(16), KPDS(22)
c        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2
c
cC     MAKE TIME DEPENDENT PDS PARAMETERS
c      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
c      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,VEG,IRET)
c      IF (IRET .NE. 0) THEN
c         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = VEG'
c         STOP
c      ELSE
cC         PRINT*, 'wrote VEG'
c      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,LAI,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = LAI'
         STOP
      ELSE
C         PRINT*, 'wrote LAI'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOD'
         STOP
      ELSE
C         PRINT*, 'wrote SNOD'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SNOC,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SNOC'
         STOP
      ELSE
C         PRINT*, 'wrote SNOC'
      END IF


      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SALBD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SALBD'
         STOP
      ELSE
C         PRINT*, 'wrote SALBD'
      END IF



      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,TMP,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = TMP'
         STOP
      ELSE
C         PRINT*, 'wrote TMP'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SPFH,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = SPFH'
         STOP
      ELSE
C         PRINT*, 'wrote SPFH'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,UGRD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = UGRD'
         STOP
      ELSE
C         PRINT*, 'wrote UGRD'
      END IF
      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,VGRD,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = VGRD'
         STOP
      ELSE
C         PRINT*, 'wrote VGRD'
      END IF
      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14),
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,PRES,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = PRES'
         STOP
      ELSE
C         PRINT*, 'wrote PRES'
      END IF

      READ (30, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
        IF(KPDS(16).ne.0) KPDS(15)=LDAS%writeintc2


C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(LDAS,TODAY, YESTERDAY, KPDS, HR1)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,ACPCP,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HR1,',','FIELD = ACPCP'
      ELSE
C         PRINT*, 'wrote ACPCP'
C         STOP
      END IF

      CALL BACLOSE (LUGB, JRET)
      CLOSE (30)
!=== Write statistical output


      END
