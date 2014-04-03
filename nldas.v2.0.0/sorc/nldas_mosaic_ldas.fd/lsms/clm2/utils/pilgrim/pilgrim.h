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

!-------------------------------------------------------------------------
! MPI Context:
!-----------------------------------------------------------------------
!
!  Now only contains cpp defines.  This file should really be called
!     mpi_defines.h
!

#define  ROOT_PE      0
#define  MAX_PE    1024
#define  MAX_TRF     10
#define  MAX_PAX MAX_PE

!
! Max buffer size for Shared Memory Arenas
!
#define  MAX_BUF  10000000

#if defined(CRAY)
#define CPP_MPI_INTEGER MPI_INTEGER
#define CPP_INTEGER i8
#else
#define CPP_MPI_INTEGER MPI_INTEGER
#define CPP_INTEGER i4
#endif

#if defined(CRAY)
#define CPP_MPI_REAL MPI_REAL
#define CPP_REAL     r8
#else
#define CPP_MPI_REAL MPI_DOUBLE_PRECISION
#define CPP_REAL     r8
#endif
#define CPP_REAL_WIDTH 8
