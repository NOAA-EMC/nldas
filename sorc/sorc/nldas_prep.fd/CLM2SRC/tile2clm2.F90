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
! tile2clm2.F90:
!
! DESCRIPTION:
!  This subroutine equates tile space variables between standard LDAS and that used in CLM2
!
! REVISION HISTORY:
!  7 Nov 2002: Jon Gottschalck; Initial code
! 12 Feb 2003: Jon Gottschalck; Modified to include tilespace fix for defined 
!                               land but no vegetation set (similar to maketiles.f fix)
!=========================================================================

  SUBROUTINE TILE2CLM2 (LDAS,TILE,GRID,WTXY,VEGXY)

  use clm_varder          ! CLM2 tile variables
  use clm_varctl, only    : nsrest, fpftcon
  use clm_varpar, only    : lsmlon, lsmlat, maxpatch, maxpatch_pft, npatch_urban, npatch_lake, npatch_wet, npatch_gla
  use clm_varsur, only    : numlon, landmask
  use ldas_module         ! LDAS variables
  use tile_module         ! Tile variables
  use grid_module         ! Grid variables  

  IMPLICIT NONE

!=== Arguments ===========================================================

  TYPE (LDASDEC)     :: LDAS
  TYPE (GRIDDEC)     :: GRID(LDAS%NC,LDAS%NR)
  TYPE (TILEDEC)     :: TILE(LDAS%NCH)

!=== Local variables

  INTEGER :: T,I,C,R,LAKEFLAG,M
  REAL(R8), INTENT(OUT) :: WTXY(LSMLON,LSMLAT,MAXPATCH)
  INTEGER,  INTENT(OUT) :: VEGXY(LSMLON,LSMLAT,MAXPATCH)
  REAL :: tsum

!=== End Local variable list

  do r=1,ldas%nr
   do c=1,ldas%nc
    if (grid(c,r)%mask == 1.0) landmask(c,r) = 1
    if (landmask(c,r) .eq. 1) then
    tsum=0.0
    do m=1,maxpatch_pft
     tsum = tsum + grid(c,r)%fgrd(m)
    enddo
    if (tsum .eq. 0.0) then
     wtxy (c,r,5) = 1.0
     vegxy(c,r,5) = 5	   
    else
     do m=1,maxpatch_pft
      wtxy (c,r,m) = grid(c,r)%fgrd(m)
      vegxy(c,r,m) = m
     enddo
    endif
    endif
    wtxy (c,r,npatch_urban) = 0.0
    vegxy(c,r,npatch_urban) = 12
!    wtxy (c,r,npatch_lake)  = grid(c,r)%fgrd(ldas%nt)
    wtxy (c,r,npatch_lake)  = 0.0
    vegxy(c,r,npatch_lake)  = 12
    wtxy (c,r,npatch_wet)   = 0.0
    vegxy(c,r,npatch_wet)   = 12
    wtxy (c,r,npatch_gla)   = 0.0
    vegxy(c,r,npatch_gla)   = 12
    
!    write(*,23) r,c,grid(c,r)%lat,grid(c,r)%lon,maxpatch,maxpatch_pft
! 23   format(i3,1x,i3,1x,f9.4,1x,f9.4,1x,i2,1x,i2)
!    do m=1,maxpatch
!     if (r .eq. 169 .and. c .eq. 1431) then
!     write(*,24) r,c,grid(c,r)%mask,landmask(c,r),(grid(c,r)%fgrd(m),m=1,13),(wtxy(c,r,m),m=1,13)
! 24  format(i3,1x,i4,1x,f3.1,1x,i2,1x,26(f3.1,1x))
!     endif
!    enddo

   enddo
  enddo

  return

end subroutine tile2clm2
