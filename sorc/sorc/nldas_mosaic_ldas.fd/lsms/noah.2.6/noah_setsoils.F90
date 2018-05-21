!-------------------------------------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
! Released October 2005
!
! See SOFTWARE DISTRIBUTION POLICY for software distribution policies
!
! The LIS source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See COPYRIGHT.TXT for copyright details.
!
!-------------------------------------------------------------------------
#include "misc.h"
!BOP
!
! !ROUTINE: noahsetsoils
!
! !DESCRIPTION:
!  This subroutine retrieves NOAH parameters - Significant F90 revisions
!   below this subroutine will be required in the future.  
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added NOAH LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!
! !INTERFACE:
subroutine noah_setsoils
! !USES:
  use lisdrv_module, only : grid,tile,lis
  use noah_varder      ! NOAH tile variables
  use lis_openfileMod
  use lis_indices_module
!EOP      
  implicit none

!+++ Local Parameters for new soil definition ++++++++++++++++++++++++++++
!+++  Layer depths correspond to soil property map files - do not change!
!      Thicknesses of soil layers for any simulations
!      These are compatible with the FIRST set of Matt Rodell's
!      soil files, i.e. sand_nlis3.5.bfsa

  real,allocatable :: stype(:,:)
  integer :: cindex, rindex
  
  integer :: n,i,j,k,jj,c,r                  !loop counters
  integer :: line,line1,line2,glnc,glnr
  integer :: ios1
  integer, allocatable :: soiltyp(:,:)
  integer, allocatable :: placesltyp(:)
  real :: basicset(noahdrv%noah_zst,noahdrv%noah_nsoilp)
  real, allocatable :: sand1(:,:)
  real, allocatable :: clay1(:,:)
!BOC  
  call lis_log_msg('MSG: noah_setsoils -- reading soil and clay files')
!-----------------------------------------------------------------------
! Open soil files (sand, clay, and porosity).
!-----------------------------------------------------------------------

  allocate(sand1(lis_nc_data,lis_nr_data))
  allocate(clay1(lis_nc_data,lis_nr_data))

!-----------------------------------------------------------------------
!     Read soil properties and convert to Zobler soil classes.
!     Note that the 3 files each contain 3 global records, one
!     for each layer depth (0-2, 2-150, 150-350 cm).  At this
!     time only the top layer data are used to map soil 
!     parameters.  Since NOAH has 4 layers, the third layer of 
!     porosity is temporarily used for layer 4 as well.
!-----------------------------------------------------------------------
  call readsand(lis%d%soil, sand1)
  call readclay(lis%d%soil, clay1)

  call lis_log_msg('MSG: noah_setsoils -- read sand and clay files')
!-----------------------------------------------------------------------
! Determine Zobler-equivalent soil classes derived from
! percentages of sand and clay, used in NOAH.
!-----------------------------------------------------------------------
  allocate(soiltyp(lis_nc_data, lis_nr_data))
  CALL SOILTYPE(lis_nc_data, lis_nr_data,SAND1,CLAY1,SOILTYP)
  deallocate(sand1)
  deallocate(clay1)
!-----------------------------------------------------------------------
! Read in the NOAH Soil Parameter File
!-----------------------------------------------------------------------
  open(unit=18,file=noahdrv%noah_sfile,status='old', & 
       access='sequential')
  
  do i=1,noahdrv%noah_nsoilp
     read(18,*)(basicset(jj,i),jj=1,noahdrv%noah_zst)
  enddo
  close(18)
  call lis_log_msg('MSG: noah_setsoils -- read sfile: '//&
                   trim(noahdrv%noah_sfile))
!-----------------------------------------------------------------------
! Convert grid space to tile space for soil type values.
!-----------------------------------------------------------------------
  allocate(placesltyp(lis%d%nch))
  do i=1,lis%d%nch
     placesltyp(i) = soiltyp(tile(i)%col, tile(i)%row-lis_tnroffset)
     noah(i)%zobsoil = placesltyp(i)
  end do

  if(lis%o%wparam.eq.1) then 
     allocate(stype(lis%d%lnc,lis%d%lnr))
     stype = -9999.0
#if ( defined OPENDAP )
     do i=1,lis%d%nch      
        stype(tile(i)%col,tile(i)%row) = noah(i)%zobsoil(1)*1.0
     enddo
#else
     do i=1,lis%d%nch      
        if(grid(i)%lat.ge.lis%d%gridDesc(4).and. & 
             grid(i)%lat.le.lis%d%gridDesc(7).and. & 
             grid(i)%lon.ge.lis%d%gridDesc(5).and. & 
             grid(i)%lon.le.lis%d%gridDesc(8)) then
           rindex = tile(i)%row - nint((lis%d%gridDesc(4)-lis%d%gridDesc(44)) &
                   /lis%d%gridDesc(9))
           cindex = tile(i)%col - nint((lis%d%gridDesc(5)-lis%d%gridDesc(45)) &
                   /lis%d%gridDesc(10))
           stype(cindex,rindex) = noah(i)%zobsoil(1)*1.0
        endif
     enddo
#endif
     
     open(32,file="soiltype.bin",form='unformatted')
     write(32) stype
     close(32)  
     deallocate(stype)
  endif
  
  deallocate(soiltyp)
!-----------------------------------------------------------------------
! Assign SOIL Parameters to each tile based on the
! type of Zobler soil class present in that tile.
!-----------------------------------------------------------------------
  do i=1,lis%d%nch            !tile loop
     k=placesltyp(i)           !soil type
     do j=1,noahdrv%noah_nsoilp   !soil parameter loop
        noah(i)%soilp(j)=basicset(k,j)
     enddo !j
  enddo !i
  deallocate(placesltyp)
  return
!EOC
end subroutine noah_setsoils
