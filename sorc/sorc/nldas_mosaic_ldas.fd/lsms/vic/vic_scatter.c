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
//
// !ROUTINE: vic_scatter.c
//
// !DESCRIPTION: 
//  Scatters VIC tiles onto compute nodes before LSM runs
//
// !INTERFACE:
void FTN(vic_scatter)()
  //EOP
{
  extern par_struct par; 
  extern ctile_spmd tspmd; 
  extern soil_con_struct *soil_con; 
#if (defined SPMD)
  MPI_Scatterv(soil_con,tspmd.cdi_array,tspmd.cdispls,
               par.MPI_SOIL_CON,soil_con,tspmd.cdi_array[par.rank],
               par.MPI_SOIL_CON,0,MPI_COMM_WORLD);
#endif 

}
