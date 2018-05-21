      SUBROUTINE READ_NOAH_GRIB2(nldas,potevap,filnam,lugb)
!
!-------------------------------------------------------------------------
! ABSTRACT: This routine is to read a grib2 data from Noah model output
!   Oct 2013     Y. Xia    
!-------------------------------------------------------------------------
!
      use grib_mod
      implicit none
 
      integer,intent(in) :: lugb
      integer,dimension(200) :: jids,jpdt,jgdt
      integer j, k, jdisc, jpdtn, jgdtn, iret, jret, kf, nldas
      
      type(gribfield) :: gfld
      logical :: unpack=.true.
      real  potevap(nldas)
      character filnam*150

!     open grib2 file
      call baopen (lugb, trim(filnam), jret)
      if (jret.ne.0) then
       print *, 'jret=',jret
        write(6,*) 'GRIB:BAOPEN ERR FOR DATA ',trim(filnam)
         write(6,*) 'PLEASE CHECK DATA AVAILABLE OR NOT'
          endif
! Set GRIB2 field identification values to search for
      j=0         ! search from 0 record
      jdisc=0    ! set discipline, for met field:0 hydro: 1, land:2
!-- set id section
      jids=-9999
!-- set product def template, using template 4.8
      jpdtn=8          ! set template number
!-- set product def template array
      jpdt=-9999
!for pdt, define catogory, parameter and level (e.g., potevap)    
      jpdt(1)=1        
      jpdt(2)=40      
      jpdt(10)=1     ! 
      jpdt(11)=0       
      jpdt(12)=0
!
!-- set grid def template
      jgdtn=-1
!-- set product def array
      jgdt=-9999      
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,unpack,k,gfld,iret)
      if ( iret.ne.0) then
            print *,' getgb2 error = ',iret
            call errexit(iret)
      endif
     
      kf = gfld%ngrdpts
      if (kf.eq.nldas) then
      potevap(1:kf)=gfld%fld(1:kf)
      else
      print *, 'Reading data problem, please check!!!'
      print *, 'nldas= ',nldas,'kf= ', kf
      endif

      call gf_free(gfld)
      call baclose(lugb,jret)
      end     