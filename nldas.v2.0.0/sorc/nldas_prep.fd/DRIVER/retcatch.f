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
! retcatch.f: 
!
! DESCRIPTION:
!  Retrieves the name of a file from getcatch.f and then gets the CATCHMENT-LDAS
!  forcing data using the zterp.f subroutine.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
! REVISION HISTORY:
!  27 Feb 2001: Brian Cosgrove; Initial code based heavily on reteta.f
!=========================================================================
      Subroutine Retcatch(ORDER,LDAS,GRID,name,ferror,ftype,
     &    prevname,count)
	  
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (griddec) GRID(LDAS%NC,LDAS%NR)   

      CHARACTER*80 NAME(9),prevname
      real*8 ftime
      real varfield(LDAS%NC,LDAS%NR)
      character*1 pvanme(67)
      character*1 dirnom(42)
      character*1 infile(24)
      character*1 dirnom1(42)
      character*1 infile1(24)
      character*24 ndummy
	character*22 validtime
      INTEGER ORDER  !NCEPDATA 1 or 2
      integer ferror
      integer ftype      !file type,1 EDAS, 2 ETA3hr,  3 ETA6hr
      integer iv
      integer tnv   !total number of variables to retrieve
      parameter(tnv=10)
      integer nineoctet(tnv),tenoctet(tnv),twen1oct(tnv)
      integer nineoctet2(tnv),tenoctet2(tnv),twen1oct2(tnv)
      integer sixoctet(tnv),sixoctet2(tnv)
      integer varadj(tnv),varadj2(tnv)
	
      real emptyvf(LDAS%NC,LDAS%NR)
      real prevvarfield(LDAS%NC,LDAS%NR)
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
        real gridsp(ldas%nc,ldas%nr)
        real *4 catchsp(5018)
	integer count



c       print*,'Begin RETCATCH'
      ferror=1
      iv=0
      goodcnt=0
! If there were a problem then ferror would
! be set to 0


      endloop=0
      iv=0

	do I=1,9
c        print *,'NAME=',NAME(I)
        open (unit=20,file=NAME(I),form='unformatted',status='old')
	do ii=1,count
        read (20) catchsp
	enddo
        close(20)

        CALL TRNSFMAARON(LDAS,GRIDSP,CATCHSP)
         do ii=1,LDAS%NC
          do jj=1,LDAS%NR
              if(ORDER.EQ.1)THEN 
                GRID(ii,jj)%CATCHDATA1(I)=GRIDSP(ii,jj)
              else
                GRID(ii,jj)%CATCHDATA2(I)=GRIDSP(ii,jj)	
	      endif  
          enddo
         enddo

        do jj=1,LDAS%NR
       do ii=1,LDAS%NC
         if(gridsp(ii,jj).eq.LDAS%UDEF.and.GRID(ii,jj)%FIMASK.gt.0)then
	gridsp(ii,jj)=150.0
          print*,'RetEta Error -- Undef Value @',ii,jj,'I=',i,
     &  gridsp(ii,jj)
         endif 
        enddo
       enddo
	ENDDO !The I=1,9 loop

c      print*,'ENDING RETcatch'
C        print*,'EVT1  ',LDAS%EVT1
C        print*,'EVT2  ',LDAS%EVT2
C        print*,'validtime  ',validtime
	
      Return
      End

