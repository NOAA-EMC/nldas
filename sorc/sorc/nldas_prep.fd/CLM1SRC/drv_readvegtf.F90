#include <misc.h>

subroutine drv_readvegtf (drv,grid,ldas_tile,tile,clm1,ldas)

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
!  The primary goal of this routine is to assign soil parameters.
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!   1 Aug 2001: Matt Rodell; Rewrite for new soil parameter definition
!=========================================================================      
! $Id: drv_readvegtf.F90,v 1.1.1.1 2003/02/06 16:10:47 jgottsch Exp $
!=========================================================================

  use precision
  use drv_module          ! 1-D Land Model Driver variables
  use tile_module         ! LDAS non-model-specific tile variables
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! 1-D CLM variables
  use drv_gridmodule      ! Grid-space variables
  use ldas_module	  ! LDAS non-model-specific 1-D variables
  implicit none

!=== Arguments ===========================================================

  type (drvdec)      :: drv              
  type (tiledec)     :: ldas_tile(drv%nch)
  type (clm_tiledec) :: tile(drv%nch)
  type (clm11d)       :: clm1 (drv%nch)
  type (clm_griddec) :: grid(drv%nc,drv%nr)   
  type (ldasdec)     :: ldas

!=== Local Variables =====================================================

  integer  :: c,r,t,i,j,l     !Loop counters
  character(15) :: vname   ! variable name read 
  integer :: ioval         ! Read error code
  real :: sand1(ldas%nc,ldas%nr)	! Layer 1 sand
  real :: clay1(ldas%nc,ldas%nr)	! Layer 1 clay
  integer :: color(ldas%nc,ldas%nr)	! Soil color

!=== End Variable Definition =============================================

! Vegetation-based soil parameterization scheme.
  if (ldas%soil .eq. 1) then

!=== Read in Vegetation Data
   open(10, file=drv%vegtf, form='formatted', status = 'old')

   ioval=0
   do while (ioval == 0)
     vname='!'
     read(10,'(a15)',iostat=ioval)vname
 
     if(vname == 'isc')   call drv_soipi(drv,ldas_tile,clm1%isoicol)        
     if(vname == 'sand')  then
        call drv_soipr(drv,ldas_tile,tile(:)%sand(1))
        do l=2,nlevsoi
            do t=1,drv%nch
               tile(t)%sand(l)=tile(t)%sand(1)
            enddo
        enddo
     endif
     if(vname == 'clay')  then
        call drv_soipr(drv,ldas_tile,tile(:)%clay(1))
        do l=2,nlevsoi
            do t=1,drv%nch
               tile(t)%clay(l)=tile(t)%clay(1)
            enddo
        enddo
     endif
   enddo     
   close(10)

  end if	!soil=1

! Reynolds soil parameterization scheme.
  if (ldas%soil .eq. 2) then

!+++ Read in soil data maps to gridded arrays
   open(11, file=ldas%safile, form='unformatted', status='old')
   open(12, file=ldas%clfile, form='unformatted', status='old')
   open(13, file=ldas%iscfile, form='unformatted', status='old')
   read(11) sand1
   read(12) clay1
   read(13) color
   close(11)
   close(12)
   close(13)

!+++ Transfer soil data to tile space
   do t=1,ldas%nch
     tile(t)%sand(1) = sand1(tile(t)%col,tile(t)%row)
     tile(t)%clay(1) = clay1(tile(t)%col,tile(t)%row)
     clm1(t)%isoicol = color(tile(t)%col,tile(t)%row)
   end do

!+++ Set soil properties in all layers equal to the top layer
!+++ NOTE THIS MUST BE IMPROVED IN THE FUTURE
   do l=2,nlevsoi
     tile%sand(l)=tile%sand(1)
     tile%clay(l)=tile%clay(1)
   enddo

  end if	!soil=2

!=== Transfer tile vegetation data to grid space variables

  call drv_t2gi(clm1%isoicol,grid%isoicol,drv%nc,drv%nr,drv%nch,tile%fgrd,tile%col,tile%row)

  do l=1,nlevsoi
     call drv_t2gr(tile(:)%sand(l),grid(:,:)%sand(l),drv%nc,drv%nr,drv%nch,tile%fgrd,tile%col,tile%row)
     call drv_t2gr(tile(:)%clay(l),grid(:,:)%clay(l),drv%nc,drv%nr,drv%nch,tile%fgrd,tile%col,tile%row)
  enddo

  return
end subroutine drv_readvegtf

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely  
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! drv_soipi.f:
!
! DESCRIPTION:
! The following subroutine simply reads and distributes spatially-constant
!  data from drv_vegp.dat into clm arrays.
!
! REVISION HISTORY:
!  6 May 1999: Paul Houser; initial code
!=========================================================================

  subroutine drv_soipi(drv,ldas_tile,clmvar)  

! Declare Modules and data structures
  use drv_module          ! 1-D Land Model Driver variables
  use tile_module         ! LDAS non-model-specific tile variables
  implicit none
  type (drvdec)               :: drv 
  type (tiledec)              :: ldas_tile(drv%nch)

  integer t
  integer clmvar(drv%nch)
  integer ivar(drv%nt)

  read(10,*)ivar
  do t=1,drv%nch
     clmvar(t)=ivar(ldas_tile(t)%soilt)     
  enddo 
  end subroutine drv_soipi      


  subroutine drv_soipr(drv,ldas_tile,clmvar)  

! Declare Modules and data structures
  use drv_module          ! 1-D Land Model Driver variables
  use tile_module         ! LDAS non-model-specific tile variables
  implicit none
  type (drvdec)               :: drv 
  type (tiledec)              :: ldas_tile(drv%nch)

  integer t
  real(r8) clmvar(drv%nch)
  real(r8) rvar(drv%nt)

  read(10,*)rvar
  do t=1,drv%nch
     clmvar(t)=rvar(ldas_tile(t)%soilt)
  enddo 
end subroutine drv_soipr      
