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
!  SPFH forcing data.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!
! REVISION HISTORY:
!  05 Sep 2001: Brian Cosgrove; Initial code based on reteta.f
!=========================================================================
      Subroutine RetETAspfh(ORDER,LDAS,GRID,name,ferror,ftype,dataflag,
     &    prevname,precflag)
	  
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (griddec) GRID(LDAS%NC,LDAS%NR)   

      CHARACTER*80 NAME,prevname
      real*8 ftime
      real varfield(LDAS%NC,LDAS%NR)
      real outvaru(LDAS%NC,LDAS%NR)
      real outvarv(LDAS%NC,LDAS%NR)
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
      parameter(tnv=3)
      integer nineoctet(tnv),tenoctet(tnv),twen1oct(tnv)
      integer nineoctet2(tnv),tenoctet2(tnv),twen1oct2(tnv)
      integer sixoctet(tnv),sixoctet2(tnv)
      integer varadj(tnv),varadj2(tnv)
	
      real emptyvf(LDAS%NC,LDAS%NR)
      real prevvarfield(LDAS%NC,LDAS%NR)
      REAL UTEMP1D(LDAS%NCOLD*LDAS%NROLD)
      integer errorflag
      integer i,j,k,li,endloop
      integer ii,jj,kk,lili
      integer badprecip
      integer goodcnt
      real outvar(LDAS%NC,LDAS%NR)
      character*6 evars(LDAS%NF)
      integer eparams(LDAS%NF)
      integer eclev(LDAS%NF)
      integer etimran(LDAS%NF)
      integer iire
	real xyz1,xyz2,xyz3
      data ndummy/'none'/


      data nineoctet/011,051,001/
      data tenoctet/116,116,116/
      data sixoctet/089,089,089/
      data twen1oct/000,000,000/
      data varadj/0,0,0/

c	get temp,spechumid,press

C       print*,'Begin RETETA'
C	 print*,'EVT1  ',LDAS%EVT1
C	 print*,'EVT2  ',LDAS%EVT2
C	 print*,'validtime  ',validtime

C For NCEP      data nineoctet2/011,051,204,205,033,034,001,061,063,100,204,61/
C For NCEP      data sixoctet2/089,089,089,089,089,089,089,089,089,089,154,155/
C For NCEP      data tenoctet2/105,105,001,001,105,105,001,001,001,001,001,001/
C For NCEP      data twen1oct2/000,000,000,000,000,000,000,004,004,000,000,004/
C For NCEP      data varadj2/0,0,0,0,0,0,0,0,0,0,0,0/
      badprecip=0
      ferror=1
      iv=0
      goodcnt=0
! If there were a problem then ferror would
! be set to 0

      errorflag=0    !should turn to 1 if there is a problem with a variable

      endloop=0
      iv=0

      do while (endloop.ne.1)
       iv=iv+1  
c	 print*,'In retrieve iv is now ',iv,'  errorflag',errorflag
C	print *,'NAME=',NAME
	print *,'calling ungribeta for',name
!       call ungribeta(LDAS,GRID,NAME,LDAS%NC,LDAS%NR,
!     &  LDAS%ncold,LDAS%nrold,
!     &  nineoctet(iv),twen1oct(iv),tenoctet(iv),
!     &  ftype,varadj(iv),
!     &  varfield,errorflag,GRID%FMASK,dataflag,
!     &  sixoctet(iv),validtime,outvaru,outvarv,utemp1d)



       if(errorflag.eq.0) then
        goodcnt=goodcnt+1
       endif
       if(errorflag.eq.1) then
         if (dataflag.eq.1) then
           endloop=1
           ferror=0
         endif
       else

       if(errorflag.eq.0) then
         if(iv.eq.1.and.ORDER.eq.1) 
     &   LDAS%EVT1=validtime
	   if(iv.eq.1.and.ORDER.eq.2) 
     &   LDAS%EVT2=validtime
	 endif
	  
	IF (IV.EQ.1) THEN
         do ii=1,LDAS%NC
          do jj=1,LDAS%NR
              if(ORDER.EQ.1)THEN 
                GRID(ii,jj)%TEMPDATA1=varfield(ii,jj)
              else
                GRID(ii,jj)%TEMPDATA2=varfield(ii,jj)	
              endif              
          enddo
         enddo
        ELSEIF (IV.EQ.2) THEN
         do ii=1,LDAS%NC
          do jj=1,LDAS%NR
              if(ORDER.EQ.1)THEN
                GRID(ii,jj)%QDATA1=varfield(ii,jj)
              else
                GRID(ii,jj)%QDATA2=varfield(ii,jj)    
              endif
          enddo
         enddo
        ELSEIF (IV.EQ.3) THEN
         do ii=1,LDAS%NC
          do jj=1,LDAS%NR
              if(ORDER.EQ.1)THEN
                GRID(ii,jj)%PRESDATA1=varfield(ii,jj)
              else
                GRID(ii,jj)%PRESDATA2=varfield(ii,jj)       
              endif
          enddo
         enddo
	ENDIF

C
       endif


      IF (ERRORFLAG.EQ.1) then
       PRINT *,'COULDNT FIND PARAMETER #',IV,NAME
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

       if(dataflag.eq.1.and.iv.eq.10) endloop=1
      enddo   !the while loop

C      print*,'ENDING RETETA'
C        print*,'EVT1  ',LDAS%EVT1
C        print*,'EVT2  ',LDAS%EVT2
C        print*,'validtime  ',validtime
      Return
      End

