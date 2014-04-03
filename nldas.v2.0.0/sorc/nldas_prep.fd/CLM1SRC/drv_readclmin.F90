#include <misc.h>

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:	  
!  1-D user defined CLM parameters.
!  NOTE: For 2-D runs, it is recommended that drv_readclmin.f or drv_main.f be 
!  modified to allow the specification of spatially-variable parameters.  
!  This file is still used for all spatially-constant variabiles when
!  running CLM in 2-D. 
!
!  This subroutine works entirely in grid space.  It assigns initial spatially
!  constant values to the entire model domain in grid space based on values
!  contained in the drv_clmin.dat file.  If spatially variable grid fields are
!  to be read in, then this should be done near the end of this routine.
!  if spatially-variable tile space fields are desired, then they should be 
!  read in at ... PRH
!
!  NOTE on INDEXES: There are several index soil and vegetation values.  These
!  are followed by several parameters that are defined by the index.  If a
!  -999.0 is placed in the fields after the index, then the default index
!  value will be used.  However, if it is desired to override this default
!  value, simply replace the -999.0 with the value you desire.
!
! INPUT DATA FORMAT:
!  FORTRAN PARAMETER NAME, VALUE, description (not read in)
!  This is free format, in any order.  readclmin.f skips any comment lines
!
! REVISION HISTORY:
!   6 May 1999: Paul Houser; initial code
!  15 Jan 2000: Paul Houser; significant revision for new CLM version
!   3 Mar 2000: Jon Radakovich; Revision for diagnostic output
!  27 Nov 2001: Jon Gottschalck; Added LAI flag
!=========================================================================

subroutine drv_readclmin(ldas,ldas_tile,ldas_grid,drv,tile,grid,clm1)
  use precision
  use ldas_module      ! LDAS non-model-specific 1-D variables
  use tile_module      ! LDAS non-model-specific tile variables
  use grid_module      ! LDAS non-model-specific grid variables
  use drv_module       ! 1-D Land Model Driver variables
  use drv_tilemodule   ! CLM Tile-space variables
  use drv_gridmodule   ! CLM Grid-space variables
  use clm1type          ! 1-D CLM variable

!=== Arguments ===========================================================

  implicit none
  type (ldasdec)     ::     ldas
  type (tiledec)     ::     ldas_tile(ldas%nch)
  type (griddec)     ::     ldas_grid(ldas%nc,ldas%nr)
  type (drvdec)      ::     drv   
  type (clm_tiledec) ::     tile(ldas%nch)        
  type (clm_griddec) ::     grid(ldas%nc,ldas%nr)
  type (clm11d)       ::     clm1(ldas%nch) 

!=== Local Variables =====================================================

  integer :: c,r,t,ierr    ! Column and Row dimensions [-]
  character(15) :: vname   ! variable name read from clm1_in.dat
  integer :: ioval         ! Read error code

!=== End Variable List ===================================================

!=== Reading of the clm1_in file occurs in 2 steps:
!===  (1) Read in CLM Domain, and allocate module sizes
!===  (2) Read in CLM parameters and assign to grid module variables

!===  (1) Open and read domain info into drv_module from clm1_in input file

!=== Equate ldas_module variables to the corresponding CLM drv_module variables
! CLM domain
     drv%nch=ldas%nch
     drv%nc=ldas%nc
     drv%nr=ldas%nr
     drv%nt=ldas%nt
     drv%maxt=ldas%maxt
     drv%mina=ldas%mina
     drv%udef=ldas%udef
     drv%vclass=ldas%vclass
! CLM files
     drv%vegtf=ldas%clm1_vfile
     drv%vegpf=ldas%clm1_mvfile
     drv%poutf1d=ldas%clm1_pfile
     drv%rstf=ldas%clm1_rfile
 
! Run timing parameters
     drv%startcode=ldas%startcode
     drv%ts=ldas%ts
     drv%sss=ldas%sss
     drv%smn=ldas%smn
     drv%shr=ldas%shr
     drv%sda=ldas%sda
     drv%smo=ldas%smo
     drv%syr=ldas%syr
     drv%ess=ldas%ess
     drv%emn=ldas%emn
     drv%ehr=ldas%ehr
     drv%eda=ldas%eda
     drv%emo=ldas%emo
     drv%eyr=ldas%eyr

! IC Source: (1) restart file, (2) drv_clm1in.dat (this file)
     drv%clm_ic=ldas%clm1_ic 

! CLM initial conditions (Read into 1D drv_module variables)
     drv%t_ini=ldas%clm1_it
     drv%h2osno_ini=ldas%clm1_iscv
     drv%sw_ini=ldas%clm1_ism

!=== Equate tile_module variables to the corresponding CLM drv_tilemodule variables
     tile%col=ldas_tile%col     
     tile%row=ldas_tile%row
     tile%fgrd=ldas_tile%fgrd
     tile%vegt=ldas_tile%vegt
     tile%pveg=ldas_tile%pveg

!=== Equate tile_module variables to the corresponding clm_module variables
     do t=1,drv%nch
      clm1(t)%laiflag=ldas%lai
      clm1(t)%latdeg=ldas_tile(t)%lat
      clm1(t)%londeg=ldas_tile(t)%lon
      clm1(t)%lat=clm1(t)%latdeg*4.*atan(1.)/180.
      clm1(t)%lon=clm1(t)%londeg*4.*atan(1.)/180.
     enddo

     do t=1,drv%nch
        grid(tile(t)%col,tile(t)%row)%latdeg=clm1(t)%latdeg
        grid(tile(t)%col,tile(t)%row)%londeg=clm1(t)%londeg
     enddo !t
!                 clm1(drv%nch)%londeg   = grid(c,r)%londeg    !Longitude of tile (degrees)
!                 clm1(drv%nch)%latdeg   = grid(c,r)%latdeg    !Latitude of tile (degrees
!                 clm1(drv%nch)%lat = clm1(drv%nch)%latdeg*4.*atan(1.)/180. !tile latitude  (radians)
!                 clm1(drv%nch)%lon = clm1(drv%nch)%londeg*4.*atan(1.)/180. !tile longitude (radians)

!=== Equate grid_module variables to the corresponding CLM drv_gridmodule variables
     do r=1,drv%nr  !rows
      do c=1,drv%nc  !columns
       allocate (grid(c,r)%fgrd(drv%nt))
       allocate (grid(c,r)%pveg(drv%nt))
       grid(c,r)%mask=ldas_grid(c,r)%imask
      enddo !R
     enddo !C

     do t=1,drv%nt
      do r=1,drv%nr  !rows
       do c=1,drv%nc  !columns
        grid(c,r)%fgrd(t)=ldas_grid(c,r)%fgrd(t)
        grid(c,r)%pveg(t)=ldas_grid(c,r)%pveg(t)
       enddo !R
      enddo !C
     enddo !t
!=== Open and read 1-D  CLM input file
  open(10, file=ldas%clm1_cfile, form='formatted', status = 'old')
  ioval=0
  do while (ioval == 0)
     vname='!'
     read(10,'(a15)',iostat=ioval)vname
     if (vname == 'metf1d')          call drv_get1dcvar(drv%metf1d)
     if (vname == 'outf1d')          call drv_get1dcvar(drv%outf1d) 

     if (vname == 'surfind')         call drv_get1divar(drv%surfind)
     if (vname == 'soilind')         call drv_get1divar(drv%soilind)
     if (vname == 'snowind')         call drv_get1divar(drv%snowind)
  enddo
  close(10)

  open(10, file=ldas%clm1_cfile, form='formatted', status = 'old')
  ioval=0
  do while (ioval == 0)
     vname='!'
     read(10,'(a15)',iostat=ioval)vname
     c=drv%nc
     r=drv%nr

! CLM Forcing parameters (read into 2-D grid module variables)
     if (vname == 'forc_hgt_u')      call drv_get2drvar(c,r,grid%forc_hgt_u)
     if (vname == 'forc_hgt_t')      call drv_get2drvar(c,r,grid%forc_hgt_t)
     if (vname == 'forc_hgt_q')      call drv_get2drvar(c,r,grid%forc_hgt_q)

! CLM Vegetation parameters (read into 2-D grid module variables)

     if (vname == 'dewmx')           call drv_get2drvar(c,r,grid%dewmx)
     if (vname == 'rootfr')          call drv_get2drvar(c,r,grid%rootfr)

! CLM Soil parameters	(read into 2-D grid module variables)

     if (vname == 'smpmax')          call drv_get2drvar(c,r,grid%smpmax)
     if (vname == 'scalez')          call drv_get2drvar(c,r,grid%scalez)
     if (vname == 'hkdepth')         call drv_get2drvar(c,r,grid%hkdepth)
     if (vname == 'wtfact')          call drv_get2drvar(c,r,grid%wtfact)

! Roughness lengths (read into 2-D grid module variables)

     if (vname == 'zlnd')            call drv_get2drvar(c,r,grid%zlnd)
     if (vname == 'zsno')            call drv_get2drvar(c,r,grid%zsno)
     if (vname == 'csoilc')          call drv_get2drvar(c,r,grid%csoilc)

! Numerical finite-difference parameters (read into 2-D grid module variables)

     if (vname == 'capr')            call drv_get2drvar(c,r,grid%capr)
     if (vname == 'cnfac')           call drv_get2drvar(c,r,grid%cnfac)
     if (vname == 'smpmin')          call drv_get2drvar(c,r,grid%smpmin)
     if (vname == 'ssi')             call drv_get2drvar(c,r,grid%ssi)
     if (vname == 'wimp')            call drv_get2drvar(c,r,grid%wimp)
     if (vname == 'pondmx')          call drv_get2drvar(c,r,grid%pondmx)

  enddo
  close(10)

!=== Open Files (to be read in later)
!  open(11,file=drv%metf1d, form='formatted')  !Meteorological Input

! If restarting from a restart file then assume append to old output file
  if (drv%startcode == 1)then  !Append to old output file
!     open(20,file=drv%outf1d,  form='unformatted',position='append')  !Timeseries output
     open(21,file='leaf_err.dat', form='unformatted',position='append')
  else
!     open(20,file=drv%outf1d,  form='unformatted')  !Timeseries output
     open(21,file='leaf_err.dat', form='unformatted')
!     open(57,file= 'alma_'//drv%outf1d,form='unformatted') ! ALMA
  endif

!=== Set clm diagnostic indices and allocate space

  clm1%surfind = drv%surfind 
  clm1%soilind = drv%soilind
  clm1%snowind = drv%snowind

  do t=1,drv%nch 
     allocate (clm1(t)%diagsurf(1:drv%surfind             ),stat=ierr); call drv_astp(ierr) 
     allocate (clm1(t)%diagsoil(1:drv%soilind,1:nlevsoi   ),stat=ierr); call drv_astp(ierr)
     allocate (clm1(t)%diagsnow(1:drv%snowind,-nlevsno+1:0),stat=ierr); call drv_astp(ierr)
  end do


!=== Read in 2-D and 3-D (GRID SPACE) parameter arrays here (to overwrite 1-D arrays read above)
!===  NOTE TO USER: READ IN YOUR 2-D PARAMETERS & INITIAL CONDITIONS HERE

end subroutine drv_readclmin

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
! The following subroutine simply reads and distributes spatially-constant
!  data from clm1_in.dat into clm arrays.
!
! REVISION HISTORY:
!  6 May 1999: Paul Houser; initial code
!=========================================================================

subroutine drv_get2divar(nc,nr,clmvar)  
  implicit none
  character*15 vname  
  integer nc,nr,x,y
  integer clmvar(nc,nr)
  integer ivar

  backspace(10)
  read(10,*)vname,ivar
  do x=1,nc
     do y=1,nr
        clmvar(x,y)=ivar
     enddo
  enddo
end subroutine drv_get2divar

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
! The following subroutine simply reads and distributes spatially-constant
!  data from clm_in.dat into clm arrays.
!
! REVISION HISTORY:
!  6 May 1999: Paul Houser; initial code
!=========================================================================

subroutine drv_get2drvar(nc,nr,clmvar)  
  use precision
  implicit none  
  character*15 vname  
  integer nc,nr,x,y
  real(r8) clmvar(nc,nr)
  real(r8) rvar

  backspace(10)
  read(10,*)vname,rvar
  do x=1,nc
     do y=1,nr
        clmvar(x,y)=rvar
     enddo
  enddo
end subroutine drv_get2drvar

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
! The following subroutine simply reads data from clm_in.dat drv module.
!
! REVISION HISTORY:
!  6 May 1999: Paul Houser; initial code
!=========================================================================

subroutine drv_get1divar(drvvar)  
  use precision
  implicit none
  character*15 vname  
  integer drvvar

  backspace(10)
  read(10,*)vname,drvvar
end subroutine drv_get1divar

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
! The following subroutine simply reads data from clm_in.dat drv module.
!
! REVISION HISTORY:
!  6 May 1999: Paul Houser; initial code
!=========================================================================

subroutine drv_get1drvar(drvvar)  
  use precision
  implicit none
  character*15 vname  
  real(r8) drvvar

  backspace(10)
  read(10,*)vname,drvvar
end subroutine drv_get1drvar

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
! The following subroutine simply reads data from clm_in.dat drv module.
!
! REVISION HISTORY:
!  6 May 1999: Paul Houser; initial code
!=========================================================================

subroutine drv_get1dcvar(drvvar)  
  use precision
  implicit none
  character*15 vname  
  character*40 drvvar

  backspace(10)
  read(10,*)vname,drvvar
end subroutine drv_get1dcvar





