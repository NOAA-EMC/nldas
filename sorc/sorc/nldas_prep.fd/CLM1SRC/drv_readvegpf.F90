#include <misc.h>

subroutine drv_readvegpf (ldas,drv,grid,tile,clm1)

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
!   Read in vegetation class paramters from input file and assign to
!   CLM variables.
!
! INPUT DATA FORMAT:
!  FORTRAN PARAMETER NAME, description (not read in)
!  values (number of types in vegetation classification)
!  
!  This is free format, in any order.  drv_readvegp.f skips any comment lines
!
! REVISION HISTORY:
!   6 May 1999: Paul Houser; initial code
!  15 Jan 2000: Paul Houser; revised for F90
!  27 Nov 2001: Jon Gottschalck; Added call to LAIREAD (AVHRR LAI)
!=========================================================================
! $Id: drv_readvegpf.F90,v 1.1.1.1 2003/02/06 16:10:47 jgottsch Exp $
!=========================================================================

  use precision
  use ldas_module         ! LDAS run variables
  use drv_module          ! 1-D Land Model Driver variables
  use drv_gridmodule      ! Grid space module
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! 1-D clm1 variables
  use clm1_varcon, only : istwet, istice, istdlak , istslak
  implicit none

!=== Arguments ===========================================================

  type (ldasdec)     :: ldas
  type (drvdec)      :: drv              
  type (clm_griddec) :: grid(drv%nc,drv%nr)
  type (clm_tiledec) :: tile(drv%nch)
  type (clm11d)       :: clm1 (drv%nch)

!=== Local Variables =====================================================

  character(15) :: vname   ! variable name read from clm_in.dat
  integer :: ioval,t       ! Read error code; tile space counter
  integer :: tmp           ! temporary for soil type determination

!=== End Variable List ===================================================

! Open and read 1-D  CLM input file

  open(10, file=drv%vegpf, form='formatted', status = 'old')

  ioval=0
  do while (ioval == 0)

     vname='!'
     read(10,'(a15)',iostat=ioval)vname
     if (vname == 'itypwat'  ) call drv_vpi(drv,tile,clm1%itypwat)
     if (vname == 'lai0')      call drv_vpr(drv,tile,clm1%minlai) 
     if (vname == 'lai')       call drv_vpr(drv,tile,clm1%maxlai) 
     clm1%tlai=clm1%maxlai
     if (vname == 'sai')       call drv_vpr(drv,tile,clm1%tsai  )
     if (vname == 'z0m')       call drv_vpr(drv,tile,clm1%z0m   )
     if (vname == 'displa')    call drv_vpr(drv,tile,clm1%displa)
     if (vname == 'dleaf')     call drv_vpr(drv,tile,clm1%dleaf )
     if (vname == 'roota')     call drv_vpr(drv,tile,tile%roota)
     if (vname == 'rootb')     call drv_vpr(drv,tile,tile%rootb)
     if (vname == 'rhol_vis')  call drv_vpr(drv,tile,clm1%rhol(1))
     if (vname == 'rhol_nir')  call drv_vpr(drv,tile,clm1%rhol(2))
     if (vname == 'rhos_vis')  call drv_vpr(drv,tile,clm1%rhos(1))
     if (vname == 'rhos_nir')  call drv_vpr(drv,tile,clm1%rhos(2))
     if (vname == 'taul_vis')  call drv_vpr(drv,tile,clm1%taul(1))
     if (vname == 'taul_nir')  call drv_vpr(drv,tile,clm1%taul(2))
     if (vname == 'taus_vis')  call drv_vpr(drv,tile,clm1%taus(1))
     if (vname == 'taus_nir')  call drv_vpr(drv,tile,clm1%taus(2))
     if (vname == 'xl')        call drv_vpr(drv,tile,clm1%xl)
     if (vname == 'htop')      call drv_vpr(drv,tile,clm1%htop)
     if (vname == 'hbot')      call drv_vpr(drv,tile,clm1%hbot)
     if (vname == 'qe25')      call drv_vpr(drv,tile,clm1%qe25)
     if (vname == 'ko25')      call drv_vpr(drv,tile,clm1%ko25)
     if (vname == 'kc25')      call drv_vpr(drv,tile,clm1%kc25)
     if (vname == 'vcmx25')    call drv_vpr(drv,tile,clm1%vcmx25)
     if (vname == 'ako')       call drv_vpr(drv,tile,clm1%ako)
     if (vname == 'akc')       call drv_vpr(drv,tile,clm1%akc)
     if (vname == 'avcmx')     call drv_vpr(drv,tile,clm1%avcmx)
     if (vname == 'bp')        call drv_vpr(drv,tile,clm1%bp)
     if (vname == 'mp')        call drv_vpr(drv,tile,clm1%mp)
     if (vname == 'folnmx')    call drv_vpr(drv,tile,clm1%folnmx)
     if (vname == 'folnvt')    call drv_vpr(drv,tile,clm1%folnvt)
     if (vname == 'c3psn')     call drv_vpr(drv,tile,clm1%c3psn)

! initialize lakpoi from itypwat variable

     do t=1,drv%nch 

        if (clm1(t)%itypwat == istdlak .or. clm1(t)%itypwat == istslak) then
           clm1(t)%lakpoi = .true.
        else
           clm1(t)%lakpoi = .false.
        endif

!        if (tile(t)%vegt == 18) then  !bare soil index
! Change for LDAS UMD parameters
        if (tile(t)%vegt == 12) then  !bare soil index
           clm1(t)%baresoil = .true.
        else
           clm1(t)%baresoil = .false.
        endif

        clm1(t)%irrig = .false.  !for now - no irrigation 

     end do

  enddo
  close(10)

!=== If using AVHRR LAI make a call to read in LAI

  IF (LDAS%LAI .EQ. 2) CALL CLM1LAIREAD(LDAS,DRV,TILE,CLM1)

end subroutine drv_readvegpf

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! drv_vp.f:
!
! DESCRIPTION:
! The following subroutine simply reads and distributes spatially-constant
!  data from drv_vegp.dat into clm arrays.
!
! REVISION HISTORY:
!  6 May 1999: Paul Houser; initial code
!=========================================================================

subroutine drv_vpi(drv,tile,clmvar)  

! Declare Modules and data structures
  use drv_module          ! 1-D Land Model Driver variables
  use drv_tilemodule      ! Tile-space variables
  implicit none
  type (drvdec)               :: drv              
  type (clm_tiledec)          :: tile(drv%nch)

  integer t
  integer clmvar(drv%nch)
  integer ivar(drv%nt)

  read(10,*)ivar
  do t=1,drv%nch
     clmvar(t)=ivar(tile(t)%vegt)
  enddo

end subroutine drv_vpi      


subroutine drv_vpr(drv,tile,clmvar)  

! Declare Modules and data structures
  use drv_module          ! 1-D Land Model Driver variables
  use drv_tilemodule      ! Tile-space variables
  implicit none
  type (drvdec)               :: drv              
  type (clm_tiledec)          :: tile(drv%nch)

  integer t
  real(r8) clmvar(drv%nch)
  real(r8) rvar(drv%nt)

  read(10,*)rvar
  do t=1,drv%nch
     clmvar(t)=rvar(tile(t)%vegt)
  enddo

end subroutine drv_vpr      








