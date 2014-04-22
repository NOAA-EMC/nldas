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
! retncep.f: 
!
! DESCRIPTION:
!  Retrieves the name of a file from getncep.f and then gets the NCEP-LDAS
!  forcing data using the zterp.f subroutine.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!  The strategy for missing data is to backwards up to 10 days to get
!  forcing at the same time of day.
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  11 Apr 2000: Brian Cosgrove; changed code to use Forcing Mask (With inland
!               water filled in).  Deleteted unused variables.
!  27 Apr 2000: Brian Cosgrove; changed code to use the original 
!               mask again since that is the
!               mask which  NCEP has already applied to the forcing data
!               by the time NASA gets it......not possible to use the 
!               expanded NASA forcing mask
!  1  May 2000: Brian Cosgrove; changed code so that if parameter 11 (sw)
!               is not found in hourly ncep data, it will just use
!               edas-based shortwave from the hourly ncep files
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!                not a hard-wired undefined value of -999.9 and -999.0
!  18 Aug 2000: Brian Cosgrove; changed code so that FMASK and not MASK
!               is used when ungribbing.  NCEP data already has a mask applied
!               to it and so may not be able to supply forcing data to
!               all LDAS land forcing points.  In areas where LDAS
!               forcing mask states that land exists, but where NCEP forcing
!               data is non-existant, assign undefined value to forcing data.
!  22 Aug 2000: Brian Cosgrove; Altered code for US/Mexico/Canada Mask
!  05 Sep 2001: Brian Cosgrove; Removed dirnom and infile variables, changed
!               call to ungribncep to match removal.  Added code to make use
!               of precip weighting mask
!=========================================================================
      Subroutine RetNCEP(ORDER,LDAS,GRID,name,ferror,ftype,dataflag)
	  
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (griddec) GRID(LDAS%NC,LDAS%NR)   

      CHARACTER*80 NAME
      real*8 ftime
      real varfield(LDAS%NC,LDAS%NR)
      INTEGER ORDER  !NCEPDATA 1 or 2
      integer dataflag   !1=edas/eta, 0=ncep/hourly
      integer ferror
      integer ftype      !file type, EDAS,ETA3hr, ETA6hr
      integer iv
      integer tnv   !total number of variables to retrieve
      parameter(tnv=12)
      integer nineoctet2(tnv),tenoctet2(tnv),twen1oct2(tnv)
      integer sixoctet2(tnv)
      integer varadj2(tnv)
	
      real emptyvf(LDAS%NC,LDAS%NR)
      integer errorflag
      integer i,j,k,li,endloop
      integer ii,jj,kk,lili
      integer badprecip
      integer goodcnt
      real outvar(LDAS%NC,LDAS%NR)
	real maskdiff(LDAS%NC,LDAS%NR)
      real xrat
      integer icycold
      character*6 evars(LDAS%NF)
      integer eparams(LDAS%NF)
      integer eclev(LDAS%NF)
      integer etimran(LDAS%NF)
      integer iire


      data nineoctet2/011,051,204,205,033,034,001,061,063,100,204,61/
      data sixoctet2/089,089,089,089,089,089,089,089,089,089,154,155/
      data tenoctet2/105,105,001,001,105,105,001,001,001,001,001,001/
      data twen1oct2/000,000,000,000,000,000,000,004,004,000,000,004/
      data varadj2/0,0,0,0,0,0,0,0,0,0,0,0/
      badprecip=0
 
      ferror=1
      iv=0
      goodcnt=0
! If there were a problem then ferror would
! be set to 0
      endloop=0
      iv=0
      do while (endloop.ne.1)
       iv=iv+1  
!       call ungribncep(LDAS,GRID,NAME,LDAS%NC,LDAS%NR,
!     &  LDAS%ncold,LDAS%nrold,
!     &  nineoctet2(iv),twen1oct2(iv),tenoctet2(iv),
!     &  ftype,varadj2(iv),
!     &  varfield,emptyvf,errorflag,GRID%FMASK,dataflag,
!     &  sixoctet2(iv))

        if(iv.eq.1) then
	 do j=1,LDAS%NR
	  do i=1,LDAS%NC
	   if (varfield(i,j).eq.0) then
	     maskdiff(i,j)=GRID(i,j)%FMASK-varfield(i,j)
	     if(maskdiff(i,j).ge.1) then
                maskdiff(i,j)=1
	     endif
	   else
	     maskdiff(i,j)=0
	   endif
	  enddo
	 enddo
	endif

         do j=1,LDAS%NR
          do i=1,LDAS%NC
            if (maskdiff(i,j).eq.1) then
              varfield(i,j)=LDAS%UDEF
	    endif
	  enddo
 	 enddo

       icycold=0

       if(errorflag.eq.0) then
        goodcnt=goodcnt+1
       endif

       do i=1,LDAS%NC
        do j=1,LDAS%NR

        if(iv.eq.3.or.iv.eq.11) then
         if(varfield(i,j).gt.10000.0) then
          icycold=iv
         endif
        endif
		 
        enddo
       enddo
	
       if(icycold.ne.0) then
         print*,'ICYCOLD ICYCOLD ',icycold,'  ICYCOLD'
       endif	 
		          
       if(errorflag.eq.1) then
        if (dataflag.eq.1) then
         endloop=1
         ferror=0
        endif
       else
        do ii=1,LDAS%NC
         do jj=1,LDAS%NR
          if(iv.le.7.or.iv.eq.10.or.iv.eq.11) then    !i.e. a non-precip value
           if(iv.ne.11) then
            if(ORDER.EQ.1)THEN 
              GRID(ii,jj)%NCEPDATA1(iv)=varfield(ii,jj)	
            else
              GRID(ii,jj)%NCEPDATA2(iv)=varfield(ii,jj)	
            endif              
           else
	   if (varfield(ii,jj).eq.-999) varfield(ii,jj)=LDAS%UDEF
           if(varfield(ii,jj).ne.LDAS%UDEF) then
            if(ORDER.EQ.1)THEN
             GRID(II,JJ)%NCEPDATA1(3)=varfield(ii,jj)    !replace old shortwave
            ELSE 
             GRID(II,JJ)%NCEPDATA2(3)=varfield(ii,jj)    !replace old shortwave
           endif  
          endif
	endif
         else   !it is precip  (i.e.  8,9, or 12)
         if(iv.ne.12) then
          if(varfield(ii,jj).ge.0.00) then
           if(ORDER.EQ.1)THEN
            GRID(ii,jj)%NCEPDATA1(iv)=varfield(ii,jj)
           ELSE
            GRID(ii,jj)%NCEPDATA2(iv)=varfield(ii,jj)
           ENDIF
          else
          if(ORDER.EQ.1)THEN
            GRID(ii,jj)%NCEPDATA1(iv)=0.0
          ELSE
            GRID(ii,jj)%NCEPDATA2(iv)=0.0
          ENDIF
         badprecip=iv
         endif
        else
		     
!  Okay, we're about to replace the total precip with a new total
!  SO we must find out what percentage of the Total
!  the convective is.

        IF(ORDER.EQ.1)THEN
         if(GRID(ii,jj)%NCEPDATA1(8).ne.0.0.and.
     &    GRID(ii,jj)%NCEPDATA1(8).ne.LDAS%UDEF.and.
     &    GRID(ii,jj)%NCEPDATA1(9).ne.LDAS%UDEF) then
          xrat=GRID(II,JJ)%NCEPDATA1(9)/
     &    GRID(II,JJ)%NCEPDATA1(8)	
            if(xrat.gt.1.0) xrat=1.0
            if(xrat.lt.0.0) xrat=0.0
           else
	    xrat=0.0
           endif			     
           if(varfield(ii,jj).ne.LDAS%UDEF.and.
     &      varfield(ii,jj).ge.0.00) then
           IF (LDAS%PRECIPMASK.EQ.0) THEN
            GRID(ii,jj)%NCEPDATA1(8)=varfield(ii,jj)  !replace old precip w/new
           ENDIF
           IF (LDAS%PRECIPMASK.EQ.1) THEN
            GRID(ii,jj)%NCEPDATA1(8)=(
     &       ((grid(ii,jj)%precipweight/16.0)*
     &       varfield(ii,jj))
     &       +( (1.0-(grid(ii,jj)%precipweight/16.0))*
     &       GRID(ii,jj)%NCEPDATA1(8)) )
           ENDIF
           IF(GRID(ii,jj)%FIMASK.EQ.0)
     &      GRID(ii,jj)%NCEPDATA1(8)=LDAS%UDEF
              GRID(ii,jj)%NCEPDATA1(9)=(xrat)* GRID(ii,jj)%NCEPDATA1(8)
c             GRID(ii,jj)%NCEPDATA1(9)=(xrat)*varfield(ii,jj)
           endif 
          ELSE
           if(GRID(ii,jj)%NCEPDATA2(8).ne.0.0.and.
     &      GRID(ii,jj)%NCEPDATA2(8).ne.LDAS%UDEF.and.
     &      GRID(ii,jj)%NCEPDATA2(9).ne.LDAS%UDEF) then
            xrat=GRID(II,JJ)%NCEPDATA2(9)/
     1      GRID(II,JJ)%NCEPDATA2(8)	
              if(xrat.gt.1.0) xrat=1.0
	      if(xrat.lt.0.0) xrat=0.0
	     else
	      xrat=0.0
	     endif			     
 		    
	     if(varfield(ii,jj).ne.LDAS%UDEF.and.
     &	      varfield(ii,jj).ge.0.00) then

           IF (LDAS%PRECIPMASK.EQ.0) THEN
            GRID(ii,jj)%NCEPDATA2(8)=varfield(ii,jj)  !replace old precip w/new
           ENDIF
           IF (LDAS%PRECIPMASK.EQ.1) THEN
            GRID(ii,jj)%NCEPDATA2(8)=(
     &       ((grid(ii,jj)%precipweight/16.0)*
     &       varfield(ii,jj))
     &       +( (1.0-(grid(ii,jj)%precipweight/16.0))*
     &       GRID(ii,jj)%NCEPDATA2(8)) )
           ENDIF
           IF(GRID(ii,jj)%FIMASK.EQ.0)
     &      GRID(ii,jj)%NCEPDATA2(8)=LDAS%UDEF
              GRID(ii,jj)%NCEPDATA2(9)=(xrat)* GRID(ii,jj)%NCEPDATA2(8)
c	      GRID(ii,jj)%NCEPDATA2(9)=(xrat)*varfield(ii,jj)
	     endif
            ENDIF
           endif !the iv.ne.12	 	 
	  endif   !the 1-7/10	 
			 
         enddo
        enddo	
       endif

        IF (IV.EQ.11) THEN
         IF (ERRORFLAG.EQ.1) THEN
          LDAS%FSOURCE(8)=0
         ELSEIF (ERRORFLAG.EQ.0) THEN
          LDAS%FSOURCE(8)=1
         ENDIF
        ENDIF

        IF (IV.EQ.12) THEN
         IF (ERRORFLAG.EQ.1) THEN
          LDAS%FSOURCE(12)=0
          LDAS%FSOURCE(9)=1
         ELSEIF (ERRORFLAG.EQ.0) THEN
          LDAS%FSOURCE(12)=1
          LDAS%FSOURCE(9)=0
         ENDIF
        ENDIF


       IF (ERRORFLAG.EQ.1.AND.IV.NE.10.AND.IV.NE.12.
     &     AND.IV.NE.11) then
	PRINT *,'COULDNT FIND PARAMETER #',IV,NAME
        ferror=0
       ELSE
        do i=1,LDAS%NC
         do j=1,LDAS%NR
c	print *,i,j,varfield(i,j),GRID(i,j)%FMASK
          if(varfield(i,j).eq.LDAS%UDEF.and.GRID(i,j)%FMASK.ge.1.0
     &  .AND.IV.NE.10.AND.IV.NE.12
     &  .AND.IV.NE.11.and.maskdiff(i,j).ne.1) then
c	print *,'in loop'
           print*,'RetNCEP Error -- Undef Value @',i,j,'IV=',iv	
          endif
         ENDDO
        ENDDO
       ENDIF

       if(dataflag.eq.1.and.iv.eq.10) endloop=1
       if(dataflag.eq.0.and.iv.eq.12) endloop=1

      enddo   !the while loop



      Return
      End

