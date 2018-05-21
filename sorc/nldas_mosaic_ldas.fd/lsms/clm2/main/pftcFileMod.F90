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

module pftcFileMod

!=======================================================================
CONTAINS
!=======================================================================

  subroutine pftconrd

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read and initialize vegetation (PFT) constants 
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: pftcFileMod.F90,v 1.7 2004/11/24 22:57:21 jim Exp $
!-----------------------------------------------------------------------

    use precision
    use clm_varpar       !lsm parameters
    use clm_varctl       !run control variables
    use pft_varcon       !vegetation type constants
    use fileutils, only : opnfil, getfil, relavu, getavu
    use spmdMod          !spmd variables and routines
    implicit none

! ------------------------ local variables ------------------------
    integer :: i,n              !loop indices
    character(len=256) :: locfn !local file name
    integer :: ier              !error code
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Set specific vegetation type values
! -----------------------------------------------------------------

!=== LDAS modification: Vegetation parameters used for the PFT are mapped
!=== offline to UMD types (13) and read in 
!    ncorn  = 15
!    nwheat = 16

! Set value for last type of tree

!    ntree = 8  !value for last type of tree

! Set value for non-vegetated

!    noveg = 0  !value

! ----------------------------------------------------------------------
! Assign unit number to file. Get local file. Open file and read PFT's.
! Close and release file.
! ----------------------------------------------------------------------

    if (masterproc) then
       write (6,*) 'Attempting to read PFT physiological data .....'
       n = getavu()
       call getfil (fpftcon, locfn, 0)
       call opnfil (locfn, n, 'f')
!=== LDAS modification: Using 13 UMD types not 16 PFTs       
       do i = 1, numpft-3
!          read (n,*)  pftname(i),              &
           read (n,60)  pftname(i),              &
                      z0mr(i)   , displar(i), dleaf(i)  , c3psn(i)  , &
                      vcmx25(i) , mp(i)     , qe25(i)   , rhol(i,1) , &
                      rhol(i,2) , rhos(i,1) , rhos(i,2) , taul(i,1) , &
                      taul(i,2) , taus(i,1) , taus(i,2) , xl(i)     , &
                      roota_par(i), rootb_par(i)
       end do
       call relavu (n)
    endif
    
60     format(A27,13x,f5.3,1x,f4.2,1x,f4.2,1x,f2.0,1x,f3.0,1x,f2.0,1x,&
              7(f4.2,1x),2(f5.3,1x),f5.2,1x,f4.1,1x,f3.1)


!=== LDAS modification: Bare ground already set offline
! ----------------------------------------------------------------------
! Define PFT zero to be bare ground
! ----------------------------------------------------------------------

!    pftname(noveg) = 'not_vegetated'
!    z0mr(noveg) = 0.
!    displar(noveg) = 0.
!    dleaf(noveg) = 0.
!    c3psn(noveg) = 1.
!    vcmx25(noveg) = 0.
!    mp(noveg) = 9.
!    qe25(noveg) = 0.
!    rhol(noveg,1) = 0.
!    rhol(noveg,2) = 0.
!    rhos(noveg,1) = 0.
!    rhos(noveg,2) = 0.
!    taul(noveg,1) = 0.
!    taul(noveg,2) = 0.
!    taus(noveg,1) = 0.
!    taus(noveg,2) = 0.
!    xl(noveg) = 0.
!    roota_par(noveg) = 0.
!    rootb_par(noveg) = 0.

#if ( defined SPMD )
! ----------------------------------------------------------------------
! pass surface data to all processors
! ----------------------------------------------------------------------

#if ( defined OPENDAP )
    call mpi_bcast (z0mr, size(z0mr), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (displar, size(displar), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (dleaf, size(dleaf), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (c3psn, size(c3psn), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (vcmx25, size(vcmx25), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (mp, size(mp), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (qe25, size(qe25), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (rhol, size(rhol), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (rhos, size(rhos), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (taul, size(taul), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (taus, size(taus), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (xl, size(xl), MPI_REAL, 0, MPI_COMM_WORLD, ier)
    call mpi_bcast (roota_par, size(roota_par), MPI_REAL,0, MPI_COMM_WORLD, ier)
    call mpi_bcast (rootb_par, size(rootb_par), MPI_REAL,0, MPI_COMM_WORLD, ier)
!
!    call mpi_bcast (z0mr, size(z0mr), mpir8, 0, mpicom, ier)
!    call mpi_bcast (displar, size(displar), mpir8, 0, mpicom, ier)
!    call mpi_bcast (dleaf, size(dleaf), mpir8, 0, mpicom, ier)
!    call mpi_bcast (c3psn, size(c3psn), mpir8, 0, mpicom, ier)
!    call mpi_bcast (vcmx25, size(vcmx25), mpir8, 0, mpicom, ier)
!    call mpi_bcast (mp, size(mp), mpir8, 0, mpicom, ier)
!    call mpi_bcast (qe25, size(qe25), mpir8, 0, mpicom, ier)
!    call mpi_bcast (rhol, size(rhol), mpir8, 0, mpicom, ier)
!    call mpi_bcast (rhos, size(rhos), mpir8, 0, mpicom, ier)
!    call mpi_bcast (taul, size(taul), mpir8, 0, mpicom, ier)
!    call mpi_bcast (taus, size(taus), mpir8, 0, mpicom, ier)
!    call mpi_bcast (xl, size(xl), mpir8, 0, mpicom, ier)
!    call mpi_bcast (roota_par, size(roota_par), mpir8, 0, mpicom, ier)
!    call mpi_bcast (rootb_par, size(rootb_par), mpir8, 0, mpicom, ier)
#endif
#endif

! ----------------------------------------------------------------------
! Return
! ----------------------------------------------------------------------

    if (masterproc) then
       write (6,*) 'Successfully read PFT physiological data'
       write (6,*)
    endif

    return 
  end subroutine pftconrd

end module pftcFileMod
