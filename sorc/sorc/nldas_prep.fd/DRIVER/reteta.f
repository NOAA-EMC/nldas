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
! reteta.f: 
!
! DESCRIPTION:
!  Retrieves the name of a file from geteta.f and then gets the ETA-LDAS
!  forcing data using the zterp.f subroutine.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!  The strategy for missing data is to backwards up to 10 days to get
!  forcing at the same time of day.
!  However, the program will also look at >12 hr forecasts
!
!  Also because some of the ETA files are forecast the precipitation
!  is an accumulation over a period, thus you need the previous file
!  so you can get the precip for that period.  This extra file
!  is passed in as well as a flag to say if you're going to need it.
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  20 Nov 1999: Jared Entin a little more revision to mesh with new geteta.f 
!  11 Apr 2000: Brian Cosgrove; Changed code to use Forcing Mask (with
!               inland water filled in)
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!                not a hard-wired undefined value of -999.9 and -999.0
!  22 Aug 2000: Brian Cosgrove; Altered code for US/Mexico/Canada Mask
!  07 Mar 2001: Brian Cosgrove; Replaced errorflag with errorflag1 and errorflag2
!=========================================================================
      Subroutine RetETA(ORDER,LDAS,GRID,name,ferror,ftype,dataflag,
     &    prevname,precflag,try)
	  
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (griddec) GRID(LDAS%NC,LDAS%NR)   

      CHARACTER*80 NAME,prevname
      real*8 ftime
	real varfieldnarr(ldas%ncold,ldas%nrold)
	real prevvarfieldnarr(ldas%ncold,ldas%nrold)
      real varfield(LDAS%NC,LDAS%NR)
      real outvaru(LDAS%NC,LDAS%NR)
      real outvarv(LDAS%NC,LDAS%NR)		
      REAL UTEMP1D(LDAS%NCOLD*LDAS%NROLD)
      character*1 pvanme(67)
      character*1 dirnom(42)
      character*1 infile(24)
      character*1 dirnom1(42)
      character*1 infile1(24)
      character*24 ndummy
	character*22 validtime
      INTEGER ORDER  !NCEPDATA 1 or 2
      integer dataflag   !1=edas/eta, 0=ncep/hourly
	integer precflag   ! flag for precip if previous time must be subtracted out
      integer ferror
      integer ftype      !file type,1 EDAS, 2 ETA3hr,  3 ETA6hr
      integer iv
      integer tnv   !total number of variables to retrieve
      parameter(tnv=21)
      integer pds5(tnv),pds7(tnv)
      integer try
      integer lubi,kf,kpds(200),kgds(200),km,ipopt(20),ibi,neta
      integer  lugb,iret,gbret,jret,jpds(200),jgds(200),lugbpre
      real emptyvf(LDAS%NC,LDAS%NR)
      real prevvarfield(LDAS%NC,LDAS%NR)
      integer errorflag1,errorflag2
      integer i,j,k,li,endloop
      integer ii,jj,kk,lili
      integer badprecip
      integer goodcnt
      real outvar(LDAS%NC,LDAS%NR)
      character*6 evars(LDAS%NF)
      integer eparams(LDAS%NF)
      integer eclev(LDAS%NF)
      integer etimran(LDAS%NF),count
      integer iire,ip,nldas,no,ibo
        parameter (nldas=103936)

	parameter (neta=349*277)
	real etatemp1d(neta)
	logical*1 lb(neta),lo(nldas)

	real xyz1,xyz2,xyz3
        INTEGER MXSIZE          ! Set for size of largest GRIB record (number of gridpoints)
        PARAMETER (MXSIZE=260000)
	real fld(mxsize),rlat(nldas),rlon(nldas)
        INTEGER KGDSLDAS(200)   ! Storage array for GDS elements
	REAL ETALDAS1D(NLDAS)
      data ndummy/'none'/
! standard edas/eta values are on next line
!      data nineoctet/011,051,204,205,033,034,001,061,063,157/
! in narr data, 134 is pressure, not 1 and precip (apcpsfc) is  202 not 61
!      data nineoctet/011,051,204,205,033,034,134,202,063,157/
!UGRD VGRD PRESN APCPN ACPCP CAPE TMP PRES HGT SPFH UGRD VGRD HPBL TMP SPFH FRICV SFEXC CRAIN PEVAP DSWRF DLWRF
      data pds5/33,34,134,202,63,157,11,1,7,51,33,34,
     &  221,11,51,253,208,140,228,204,205/
      data pds7/10,10,0,0,0,46080,1,1,1,1,1,
     &  1,0,2,2,0,0,0,0,0,0/

      badprecip=0
      ferror=1
      iv=0
      goodcnt=0
! If there were a problem then ferror would
! be set to 0

      errorflag1=0    !should turn to 1 if there is a problem with a variable
      errorflag2=0    !should turn to 1 if there is a problem with a variable


      endloop=0
      iv=0

      do while (endloop.ne.1)
       iv=iv+1  
!	 print*,'In retrieve iv is now ',iv,'  errorflag1',errorflag1
c	print *,'NAME=',NAME
C	print *,'nineoctet(iv),twen1oct(iv),tenoctet(iv),sixoctet(iv)',
C     &  nineoctet(iv),twen1oct(iv),tenoctet(iv),sixoctet(iv)
       !=== Set up to open file and retrieve specified field
!      Changed next line to ensure lugb varies when rolling back.
       lugb = iv + try+10
!       lugb = iv ! alternate file unit # -- value of iv not important
       lubi=0
        j=0
        jpds=-1
        jpds(5)=pds5(iv)
        jpds(7)=pds7(iv)
        jgds(1)=3
        jgds(2)=ldas%ncold
        jgds(3)=ldas%nrold
        kpds(3)=221
       call baopenr(lugb,name,iret)
        jgds=-1
       call getgb(lugb,lubi,neta,j,jpds,jgds, 
     &       kf,k,kpds,kgds,lb,varfieldnarr,gbret)

        if (gbret.ne.0) then
        print *,'problem getting field ',iv
        print *,'potential corrupt file ',name
        print *,'iret,gbret=',iret,gbret
        stop
        endif
        call baclose(lugb,jret)
        if (gbret.ne.0) errorflag1=1
c       Move data into different sized array
!         DO I=1,neta
!           ETATEMP1D(I)=FLD(I)
!         ENDDO

        DO I=1,200
         KGDSLDAS(I)=0
        ENDDO

         KGDSLDAS(1)=0
         KGDSLDAS(2)=464
         KGDSLDAS(3)=224
         KGDSLDAS(4)=25063
         KGDSLDAS(5)=-124938
         KGDSLDAS(6)=128
         KGDSLDAS(7)=52938
         KGDSLDAS(8)=-67063
         KGDSLDAS(9)=125
         KGDSLDAS(10)=125
         KGDSLDAS(11)=64
         KGDSLDAS(12)=0
         KGDSLDAS(13)=0
         KGDSLDAS(14)=0
         KGDSLDAS(15)=0
         KGDSLDAS(16)=0
         KGDSLDAS(17)=0
         KGDSLDAS(18)=0
         KGDSLDAS(19)=0
         KGDSLDAS(20)=255
         KGDSLDAS(21)=0
         KGDSLDAS(22)=0

C       set to bilinear interpolation
        IP=0
        DO I=1,20
         IPOPT(I)=0
        ENDDO
C       Set interp to nearest neighbor value 2 if crain variable
        IF (PDS5(iv).EQ.140) THEN
         IP=2
!         IP=3
!         IPOPT(1) = 4
!         IPOPT(2) = -1
!	print *,'warning setting crain interp to new value orig was 2'
        ENDIF

C       Set interp to bilinear budget if precip or potential evap variable
        IF ((PDS5(IV).EQ.202).OR.(PDS5(IV).EQ.63)
     &  .OR.(PDS5(IV).EQ.228)) THEN
         IP=3
         IPOPT(1) = 4
         IPOPT(2) = -1
        ENDIF



        KM=1
        IBI=0
         DO I=1,neta
           LB(I)=.TRUE.
         ENDDO
                CALL IPOLATES (IP,IPOPT,KGDS,KGDSLDAS,neta,NLDAS,
     +                         KM,IBI,LB,varfieldnarr,NO,RLAT,
     &                         RLON,IBO,LO,
     +                         ETALDAS1D,IRET)


	count=1
	do jj=1,LDAS%Nr
        do ii=1,LDAS%Nc
           varfield(ii,jj)=ETALDAS1D(count)
	count=count+1
	enddo
	enddo

!       call ungribeta(LDAS,GRID,NAME,LDAS%NC,LDAS%NR,
!     &  LDAS%ncold,LDAS%nrold,
!     &  nineoctet(iv),twen1oct(iv),tenoctet(iv),
!     &  ftype,varadj(iv),
!     &  varfield,errorflag1,GRID%FMASK,dataflag,
!     &  sixoctet(iv),validtime,outvaru,outvarv,utemp1d)
	if ((iv.eq.20).and.(errorflag1.eq.0)) LDAS%shortflag=2
!        if ((iv.eq.20).and.(errorflag1.eq.0)) print *,'got time avg sw'
        if ((iv.eq.21).and.(errorflag1.eq.0)) LDAS%longflag=2
!        if ((iv.eq.21).and.(errorflag1.eq.0)) print *,'got time avg lw'


C	print *,'after call, errorflag=',errorflag1
	if(iv.eq.4.or.iv.eq.5) then  !it is precip
	   if(precflag.eq.1)  then ! we must also get the previous timestep's precip to subtract out
        j=0
        jpds=-1
        jpds(5)=pds5(iv)
        jpds(7)=pds7(iv)
        lugbpre=lugb+15
       call baopenr(lugbpre,prevname,iret)
       call getgb(lugbpre,lubi,neta,j,jpds,jgds,kf,k,kpds,
     &  kgds,lb,prevvarfieldnarr,gbret)
        if (gbret.ne.0) then
        print *,'problem getting field ',iv
        print *,'potential corrupt file ',name
        stop
        endif

       call baclose(lugbpre,jret)
        if (gbret.ne.0) errorflag2=1
C       set to bilinear interpolation
        IP=0
        DO I=1,20
         IPOPT(I)=0
        ENDDO

C       Set interp to nearest neighbor if crain variable
        IF (PDS5(iv).EQ.140) THEN
         IP=2
        ENDIF

C       Set interp to bilinear budget if precip variable
        IF ((PDS5(IV).EQ.202).OR.(PDS5(IV).EQ.63)) THEN
         IP=3
         IPOPT(1) = 4
         IPOPT(2) = -1
        ENDIF

        KM=1
        IBI=0
         DO I=1,neta
           LB(I)=.TRUE.
         ENDDO
                CALL IPOLATES (IP,IPOPT,KGDS,KGDSLDAS,neta,NLDAS,
     +                         KM,IBI,LB,prevvarfieldnarr,NO,RLAT,
     &                         RLON,IBO,LO,
     +                         ETALDAS1D,IRET)

        count=1
        do jj=1,LDAS%NR
        do ii=1,LDAS%NC
           prevvarfield(ii,jj)=ETALDAS1D(count)
        count=count+1
        enddo
        enddo


!	     Call ungribeta(LDAS,GRID,PrevName,
!     &    LDAS%NC,LDAS%NR,LDAS%ncold,LDAS%nrold,
!     &    nineoctet(iv),twen1oct(iv),tenoctet(iv),	     
!     &    ftype,varadj(iv),
!     &    prevvarfield,errorflag2,GRID%FMASK,dataflag,
!     &    sixoctet(iv),validtime,outvaru,outvarv,utemp1d)
         endif
	endif
!        if(iv.eq.20) then  !it is instantaneous shortwave radiation
!        LDAS%shortflag=1       !Assume instantaneous was successfully retrieved
!         if (errorflag1.eq.1) then !it was not successfully retrieved
!	LDAS%shortflag=0       !No SW at this point, so set to 0
!C	Call ungrib for time averaged shortwave
!        call ungribeta(LDAS,GRID,NAME,LDAS%NC,LDAS%NR,
!     &  LDAS%ncold,LDAS%nrold,
!     &  204,3,001,
!     &  ftype,0,
!     &  varfield,errorflag1,GRID%FMASK,dataflag,
!     &  140,validtime,outvaru,outvarv,utemp1d)
!
!	if (errorflag1.eq.0) LDAS%shortflag=2  !Got time averaged SW, so set to 2
!
!c	if (errorflag1.eq.0) print *,'successful sw avg ret'
!c        if (errorflag1.eq.1) print *,'UNsuccessful sw avg ret'
!	 endif
!	endif

!        if(iv.eq.4) then  !it is instantaneous longwave radiation
!        LDAS%longflag=1       !Assume instantaneous was successfully retrieved
!         if (errorflag1.eq.1) then !it was not successfully retrieved
!        LDAS%longflag=0       !No LW at this point, so set to 0
!C       Call ungrib for time averaged longwave
!        call ungribeta(LDAS,GRID,NAME,LDAS%NC,LDAS%NR,
!     &  LDAS%ncold,LDAS%nrold,
!     &  205,3,001,
!     &  ftype,0,
!     &  varfield,errorflag1,GRID%FMASK,dataflag,
!     &  089,validtime,outvaru,outvarv,utemp1d)
!
!	if (errorflag1.eq.0) LDAS%longflag=2  !Got time averaged LW, so set to 2
!c        if (errorflag1.eq.0) print *,'successful lw avg ret'
!c        if (errorflag1.eq.1) print *,'UNsuccessful lw avg ret'
!         endif
!        endif



       if((errorflag1.eq.0).and.(errorflag2.eq.0)) then
        goodcnt=goodcnt+1
       endif
       if((errorflag1.eq.1).or.(errorflag2.eq.1)) then
         if (dataflag.eq.1) then
           endloop=1
           ferror=0
         endif
       else

       if((errorflag1.eq.0).and.(errorflag2.eq.0)) then
         if(iv.eq.1.and.ORDER.eq.1) 
     &   LDAS%EVT1=validtime
	   if(iv.eq.1.and.ORDER.eq.2) 
     &   LDAS%EVT2=validtime
	 endif
	   
          do jj=1,LDAS%NR
         do ii=1,LDAS%NC
            if((iv.le.3.or.iv.ge.6)) then  !i.e. a non-precip value
              if(ORDER.EQ.1)THEN 
                GRID(ii,jj)%ETADATA1(iv)=varfield(ii,jj)
              else
                GRID(ii,jj)%ETADATA2(iv)=varfield(ii,jj)	
              endif       
            else   !it is precip  (i.e.  4,5)
            
! zero out low precip values to correct for NARR issues
              if (varfield(ii,jj).le.0.01) varfield(ii,jj)=0.0
 
C Check to see if we need to subtract out the previous time's precip

             IF(precflag.eq.0) then   ! Do NOT subtract it out 
              if(varfield(ii,jj).ge.0.00) then
                IF(ORDER.EQ.1)THEN
                  GRID(ii,jj)%ETADATA1(iv)=varfield(ii,jj)
                ELSE
                  GRID(ii,jj)%ETADATA2(iv)=varfield(ii,jj)
                ENDIF
              else
                IF(ORDER.EQ.1)THEN
                  GRID(ii,jj)%ETADATA1(iv)=0.0
                ELSE
                  GRID(ii,jj)%ETADATA2(iv)=0.0
                ENDIF
                badprecip=iv
              endif
		 ELSE    !subtract out previous time's precip
C	 
		   if(varfield(ii,jj).ge.0.00.and.
     &	   prevvarfield(ii,jj).ge.0.00) then
C	   
		     IF(ORDER.EQ.1)THEN
		       GRID(ii,jj)%ETADATA1(iv)=
     &  			 varfield(ii,jj)-prevvarfield(ii,jj)
                 ELSE
		       GRID(ii,jj)%ETADATA2(iv)=
     &                   varfield(ii,jj)-prevvarfield(ii,jj)
                 ENDIF
		   else

                 IF(ORDER.EQ.1)THEN
                   GRID(ii,jj)%ETADATA1(iv)=0.0
                 ELSE
                   GRID(ii,jj)%ETADATA2(iv)=0.0
                 ENDIF		
C
		   endif

		 ENDIF			 

	      endif   !the 1-3/6-21
          enddo
         enddo
C
       endif


      IF ((ERRORFLAG1.EQ.1).or.(errorflag2.eq.1)) then
       if (errorflag1.eq.1) then
c       print *,'errorflag1=1'
       PRINT *,'bc IN NEXT HOUR COULDNT FIND PARAMETER #',IV,NAME
	stop
       endif
       if (errorflag2.eq.1) then
c       print *,'errorflag2=1'
       PRINT *,'bcc IN PREVIOUS HOUR COULDNT FIND PARAMETER #',IV,
     &  PREVNAME
	stop
       endif

      ELSE
       do i=1,LDAS%NC
        do j=1,LDAS%NR
         if(varfield(i,j).eq.LDAS%UDEF.and.GRID(i,j)%FMASK.ge.1.0)then
          print*,'RetEta Error -- Undef Value @',i,j,'IV=',iv,
     &  varfield(i,j)
        open (unit=96,file=
     &  '/SCRATCH2/LDAS-S/rad.bin',
     &  form='unformatted')
        WRITE(96)((varfield(ii,jj),ii=1,LDAS%NC),jj=1,LDAS%NR)
        close(96)

         endif 
        enddo
       enddo
      ENDIF

       if(dataflag.eq.1.and.iv.eq.21) endloop=1
      enddo   !the while loop

C      print*,'ENDING RETETA'
C        print*,'EVT1  ',LDAS%EVT1
C        print*,'EVT2  ',LDAS%EVT2
C        print*,'validtime  ',validtime
      Return
      End

