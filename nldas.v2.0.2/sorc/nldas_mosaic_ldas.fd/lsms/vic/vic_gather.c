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
#include <mpi.h>
#endif
#include "vicNl.h"
#include "ftn.h"
//BOP
//
// !ROUTINE: vic_gather.c
//
// !DESCRIPTION: 
//  Gathers VIC tiles before writing output
//
// !INTERFACE:
void FTN(vic_gather)()
//EOP  
{
  extern par_struct par; 
  extern ctile_spmd tspmd; 
  extern out_data_struct *outdata1;
  extern atmos_data_struct *atmos;
#if (defined SPMD)  
  MPI_Gatherv(outdata1,tspmd.cdi_array[par.rank],
	      par.MPI_OUT_STRUCT,outdata1,tspmd.cdi_array,tspmd.cdispls,
	      par.MPI_OUT_STRUCT,0,MPI_COMM_WORLD);

  MPI_Gatherv(atmos,tspmd.cdi_array[par.rank],
	      par.MPI_ATMOS_STRUCT,atmos,tspmd.cdi_array,tspmd.cdispls,
	      par.MPI_ATMOS_STRUCT,0,MPI_COMM_WORLD);
#endif
  
}
