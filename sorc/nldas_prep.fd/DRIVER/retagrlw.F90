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
! retagrlw.F90:
!
! DESCRIPTION:
!  Opens, reads, interpolates and overlays LW radiation forcing.
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data
!
!	1. Load AGRMET 3 level cloud amount (convertted from RTNEPH)
!	2. Convert 2 m forcing specific humidity to vapor pressure
!	3. Use AGRMET subroutine longwv() to calculate LW DOWN. Arguments are:
!	   cloud amount [%], 2 m temperature [k], and 2 m vapor pressure [pa]
!
! REVISION HISTORY:
!  28 Oct 1999: Brian Cosgrove; Initial code
!  11 Apr 2000: Brian Cosgrove; changed code to use Forcing Mask (With inland
!               water filled in).  Deleted unused variables.
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses LDAS%UDEF and
!                not a hard-wired undefined value of -999.9 and -999.0
!  25 Jun 2001: Urszula Jambor; Modified for AGRMET data use in GLDAS.
!  25 Oct 2001: Jesse Meng; Modified for AGRMET LW scheme implementation
!=========================================================================
        
subroutine retagrlw ( order, ldas, grid, yr, mo, da, hr, ferror, flag )

  use ldas_module      	! LDAS non-model-specIFic 1-D variables
  use grid_module      	! LDAS non-model-specIFic grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc, ldas%nr)
  integer :: yr, mo, da, hr
  
!=== Local Variables =====================================================

  character*4  :: cyr
  character*2  :: cmo, cda, chr
  character*80 :: nameNH, nameSH

  integer :: flag               !data source, 1=AGRMET, all we have for now
  integer :: order              !retrieve data for time1 or time2
  integer :: ferror             !0=no radiation data found
                                !1=found observational data
  
  integer :: fvalid		!0=cloud data undef

  integer :: c, r

  real :: outdata(ldas%nc, ldas%nr)

  real :: cldamtH(ldas%nc, ldas%nr)	!high layer cloud amount [%]
  real :: cldamtM(ldas%nc, ldas%nr)	!mid  layer cloud amount [%]
  real :: cldamtL(ldas%nc, ldas%nr)	!low  layer cloud amount [%]
  
  real :: cldamt1d( 3 )			!single column 3 layet cldamt [%]
  
  real :: tair				!2 m temperature [K]
  real :: vaporP			!2 m vapor pressure [Pa]

  real :: rldown			!LW down [W/m2]
  
!=== End Variable Definition =============================================

!=== If using AGRMET data, open, read in and interpolate AGRMET files
!===    to appropriate GLDAS resolution;
!=== If error reading file, or data completely missing, ferror=0;

  ferror = 0
  do r = 1, ldas%nr
  do c = 1, ldas%nc
     if (order == 1) then
        grid(c,r)%OBLWDATA1 = ldas%udef
     else if (order == 2) then
        grid(c,r)%OBLWDATA2 = ldas%udef
     end if
  end do
  end do

  if (flag == 1) then

!=== 1. Generate AGRMET cloud filenames; =================================
!===    Load 3 layers cldamt data ========================================

!  print*,"Cloud files:"

  write(cyr, '(I4.4)') yr
  write(cmo, '(I2.2)') mo
  write(cda, '(I2.2)') da
  write(chr, '(I2.2)') hr
   
  nameNH = trim(ldas%agrmdir)//'/CloudAGR/'//cyr//cmo// &
           '/NH/cldamtH_'//cyr//cmo//cda//chr//'n'
  nameSH = trim(ldas%agrmdir)//'/CloudAGR/'//cyr//cmo// &
           '/SH/cldamtH_'//cyr//cmo//cda//chr//'s'  
!  print*, trim(nameNH)
!  print*, trim(nameSH)  
  
!  call agrmet2latlon( nameNH, nameSH, cldamtH, ldas, grid, ferror )
!  print*,'LOAD cldamtH ferror =', ferror

  if ( ferror .EQ. 0 )  return

  nameNH = trim(ldas%agrmdir)//'/CloudAGR/'//cyr//cmo// &
           '/NH/cldamtM_'//cyr//cmo//cda//chr//'n'
  nameSH = trim(ldas%agrmdir)//'/CloudAGR/'//cyr//cmo// &
           '/SH/cldamtM_'//cyr//cmo//cda//chr//'s'
!  print*, trim(nameNH)
!  print*, trim(nameSH)  
  
!  call agrmet2latlon( nameNH, nameSH, cldamtM, ldas, grid, ferror )
!  print*,'LOAD cldamtM ferror =', ferror

  if ( ferror .EQ. 0 ) return
  
  nameNH = trim(ldas%agrmdir)//'/CloudAGR/'//cyr//cmo// &
           '/NH/cldamtL_'//cyr//cmo//cda//chr//'n'
  nameSH = trim(ldas%agrmdir)//'/CloudAGR/'//cyr//cmo// &
           '/SH/cldamtL_'//cyr//cmo//cda//chr//'s'
!  print*, trim(nameNH)
!  print*, trim(nameSH)  

!  call agrmet2latlon( nameNH, nameSH, cldamtL, ldas, grid, ferror )
!  print*,'LOAD cldamtL ferror =', ferror
  
  if ( ferror .EQ. 0 ) return

!=== If AGRMET cloud data is not undefined, set ferror=1

     fvalid = 0
     do c=1, ldas%nc
        do r=1, ldas%nr
           if (cldamtH(c,r) >= 0.0) fvalid = 1
           if (cldamtM(c,r) >= 0.0) fvalid = 1
           if (cldamtL(c,r) >= 0.0) fvalid = 1
        end do
     end do

!     print*,'VALIDATE cldamt, fvalid = ', fvalid

   if ( fvalid .EQ. 0 ) return
     
!=== (C,R) loops

        do c=1, ldas%nc
           do r=1, ldas%nr

           rldown = 0.

           tair = grid(c,r)%forcing(1)

!=== 2. CONVERT 2 METER SPECIFIC HUMIDITY TO VAPOR PRESSURE ==============

           vaporP = grid(c,r)%forcing(2) * grid(c,r)%forcing(7) / &
                    ( 0.622 + grid(c,r)%forcing(2) * (1-0.622) )	      
	   
!=== If tair, vaporP, and ALL 3 AGRMET LAYERS cldamt are defined, 
!===    calculate rldown; transfer to proper output array
           
           fvalid = 1
	   if (tair         .LT. 0.0) fvalid = 0
	   if (vaporP       .LT. 0.0) fvalid = 0 
           if (cldamtH(c,r) .LT. 0.0) fvalid = 0
           if (cldamtM(c,r) .LT. 0.0) fvalid = 0
           if (cldamtL(c,r) .LT. 0.0) fvalid = 0
	      
	   if (fvalid == 1) then
	      
!=== Extract the 3 layers cldamt at a single grid

              cldamt1d(1) = cldamtL(c,r)
	      cldamt1d(2) = cldamtM(c,r)
	      cldamt1d(3) = cldamtH(c,r)
     
!=== 3. Calculate lw down ================================================

              call agrlwdn( tair, vaporP, cldamt1d, rldown )
	      
              if (order == 1) then
                  grid(c,r)%OBLWDATA1 = rldown
              else if (order == 2) then
                  grid(c,r)%OBLWDATA2 = rldown
              end if

!=== If tair/vaporP/cldamt undefined
!===    fill output array with undefined values

           else	!(fvalid \= 1)

              if (order == 1) then
                  grid(c,r)%OBLWDATA1 = ldas%udef
              else if (order == 2) then
                  grid(c,r)%OBLWDATA2 = ldas%udef
              end if

           end if !(fvalid=1); tair/vaporP/cldamt defined
	      
           end do
        end do

  end if !(flag)

  return
end subroutine retagrlw

