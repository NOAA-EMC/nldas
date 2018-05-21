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
! canhtset.F90:
!
! DESCRIPTION:
!  This subroutine reads in canopy height information into CLM2 for now
!
! REVISION HISTORY:
!  15 Nov 2002: Jon Gottschalck; Initial code
!=========================================================================
 
  SUBROUTINE canhtset (pthtop,pthbot)

  use clm_varder                          ! CLM2 tile variables
  use clm_varpar     , only : maxpatch
  use clm_varmap     , only : numpatch, begpatch, endpatch

  implicit none

!=== Arguments ===========================================================

!=== Local variables

  integer :: i,t
  real    :: pthtop(numpatch),pthbot(numpatch)  ! Canopy top and bottom height for the 13 UMD vegetation types

!=== End Local variable list

!=== Open canopy heights file and read into temporary variables
  open(57,file="BCS/CLM2/clm2_ptcanhts.txt",status='old')
  read(57,*) (pthtop(i),i=1,maxpatch)
  read(57,*) (pthbot(i),i=1,maxpatch)
  close(57)
  
!=== Assign canopy top and height values for each vegetation type to CLM2 tiles
!  do t=1,numpatch
!   clm(t)%htop = pthtop(clm(t)%itypveg)
!   clm(t)%hbot = pthbot(clm(t)%itypveg)
!   print*, t, clm(t)%itypveg, clm(t)%htop, clm(t)%hbot   
!  enddo

  return

end subroutine canhtset
