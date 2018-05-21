#include <misc.h>
   subroutine endrun
!-----------------------------------------------------------------------
! Purpose:
!
! Abort the model for abnormal termination
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
! $Id: endrun.F90,v 1.1.1.1 2003/02/06 16:10:59 jgottsch Exp $
!-----------------------------------------------------------------------
#if (defined SPMD || defined COUP_CSM)
   use mpishorthand, only: MPI_COMM_WORLD
#endif
   use shr_sys_mod, only: shr_sys_flush
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
   write(6,*)'ENDRUN IS BEING CALLED'
   call shr_sys_flush( 6 )   ! Flush all output to standard output

#if (defined SPMD) || (defined COUP_CSM) 
! passing an argument of 1 to mpi_abort will lead to a STOPALL output 
! error code of 257
   call mpi_abort (MPI_COMM_WORLD, 1)  
#else
   call abort
#endif

   end
