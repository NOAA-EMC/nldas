//-------------------------------------------------------------------------
// NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
// Released October 2005
//
// See SOFTWARE DISTRIBUTION POLICY for software distribution policies
//
// The LIS source code and documentation are in the public domain,
// available without fee for educational, research, non-commercial and
// commercial purposes.  Users may distribute the binary or source
// code to third parties provided this statement appears on all copies and
// that no charge is made for such copies.
//
// NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
// SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
// IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
// LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
//
// See COPYRIGHT.TXT for copyright details.
//
//-------------------------------------------------------------------------
#include "misc.h"
#if (defined SPMD)
#include "mpi.h"
#endif
#include "vicNl.h"
#include "ftn.h"
//BOP
// !ROUTINE: vic_bcast.c
//
// !DESCRIPTION: 
// Broadcasts VIC variables
//
// !INTERFACE:
void FTN(vic_bcast)(int * ntype){
  //EOP  
  extern par_struct par; 
  extern veg_lib_struct *veg_lib; 

  //BOC 
  if ( par.rank != 0 )
  {
    veg_lib = (veg_lib_struct *)lis_calloc(*ntype,sizeof(veg_lib_struct),
                                           "vic_bcast");
  }
#if (defined SPMD)
  MPI_Bcast(veg_lib,*ntype,par.MPI_VEG_LIB,0,MPI_COMM_WORLD);
#endif
  //EOC	     
}
