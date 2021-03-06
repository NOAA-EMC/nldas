#include <misc.h>
#include <preproc.h>      

module clm_varpar

  use precision
  implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: 
! land surface model array dimensions
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: clm_varpar.F90,v 1.1.1.1 2003/02/06 16:10:50 jgottsch Exp $
!-----------------------------------------------------------------------

! Define land surface 2-d grid. This sets the model resolution according
! to cpp directives LSMLON and LSMLAT in preproc.h. 

  integer, parameter :: lsmlon = LSMLON  !maximum number of longitude points on lsm grid
  integer, parameter :: lsmlat = LSMLAT  !number of latitude points on lsm grid

! Define maximum number of PFT patches per grid cell and set
! patch number for urban, lake, wetland, and glacier patches

#if (defined DGVM)
  integer, parameter :: maxpatch_pft = 13
#else
  integer, parameter :: maxpatch_pft = 13                !maximum number of PFT subgrid patches per grid cell
#endif
  integer, parameter :: npatch_urban = maxpatch_pft + 1 !urban   patch number: 1 to maxpatch
  integer, parameter :: npatch_lake  = npatch_urban + 1 !lake    patch number: 1 to maxpatch
  integer, parameter :: npatch_wet   = npatch_lake  + 1 !wetland patch number: 1 to maxpatch
  integer, parameter :: npatch_gla   = npatch_wet   + 1 !glacier patch number: 1 to maxpatch
  integer, parameter :: maxpatch     = npatch_gla       !maximum number of subgrid patches per grid cell

! Define history file parameters

  integer , parameter :: maxhist      =   3             !max number of history files
  integer , parameter :: maxflds      = 200             !max number of fields in list
  integer , parameter :: max_slevflds =  75             !max number of active single-level fields
  integer , parameter :: max_mlevflds =  10             !max number of active multi-level fields (either snow or soil)
  integer , parameter :: maxalflds = max_slevflds + max_mlevflds !max number of active fields (all levels)

! Define number of level parameters

  integer, parameter :: nlevsoi     =  10   !number of soil layers
  integer, parameter :: nlevlak     =  10   !number of lake layers
  integer, parameter :: nlevsno     =   5   !maximum number of snow layers

! Define miscellaneous parameters

  integer, parameter :: numwat      =   5   !number of water types (soil, ice, 2 lakes, wetland)
  integer, parameter :: numpft      =  16   !number of plant types

! next variable used with DGVM
  integer, parameter :: npftpar     =  32   !number of pft parameters (in LPJ)
  integer, parameter :: numcol      =   8   !number of soil color types
  integer, parameter :: numrad      =   2   !number of solar radiation bands: vis, nir
  integer, parameter :: ndst        =   4   !number of dust size classes
  integer, parameter :: dst_src_nbr =   3   !number of size distns in src soil
  integer, parameter :: nvoc        =   4   !number of voc categories

! Define parameters for RTM river routing model

  integer, parameter :: rtmlon = 720  !# of rtm longitudes
  integer, parameter :: rtmlat = 360  !# of rtm latitudes

end module clm_varpar

