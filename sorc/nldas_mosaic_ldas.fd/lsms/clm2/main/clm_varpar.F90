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
! $Id: clm_varpar.F90,v 1.6 2004/11/24 22:57:13 jim Exp $
!-----------------------------------------------------------------------

! Define land surface 2-d grid. This sets the model resolution according
! to cpp directives LSMLON and LSMLAT in preproc.h. 

!  integer, parameter :: lsmlon = LSMLON  !maximum number of longitude points on lsm grid
!  integer, parameter :: lsmlat = LSMLAT  !number of latitude points on lsm grid

! Define maximum number of PFT patches per grid cell and set
! patch number for urban, lake, wetland, and glacier patches

#if (defined DGVM)
  integer, parameter :: maxpatch_pft = 13
#else
  integer, parameter :: maxpatch_pft = 13                !maximum number of PFT subgrid patches per grid cell
#endif
!<TODO: clean up>
#if ( defined CLM_GSWP_SUPPORT )
  integer, parameter :: maxpft = 16
#else
  integer, parameter :: maxpft = 13
#endif
  !integer, parameter :: npatch_urban = maxpatch_pft + 1 !urban   patch number: 1 to maxpatch
  integer, parameter :: npatch_urban = maxpft + 1       !urban   patch number: 1 to maxpatch
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
!<TODO: clean up>
#if ( defined CLM_GSWP_SUPPORT )
  integer, parameter :: numpft      =  19   !number of plant types
#else
  integer, parameter :: numpft      =  16   !number of plant types
#endif

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

