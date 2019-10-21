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
#include "vicNl.h"
#include "ftn.h"
#include "misc.h"
//BOP
//
// !ROUTINE: gather_model_state.c
//
// !DESCRIPTION: 
//  Gathers model state data for restart file
//
// !INTERFACE:
void FTN(gather_model_state)(int *snowband,
                             int *nlayer,
                             int *nnode)
{
//EOP
   extern par_struct par;
   extern ctile_spmd tspmd;
   extern dist_prcp_struct *prcp; 
   extern model_state_struct *model_state; 
 
   cell_data_struct     ***cell;
   snow_data_struct      **snow;
   energy_bal_struct     **energy; 
   veg_var_struct       ***veg_var;
 
   int i, veg, band, dist, lidx, nidx; 

   int Ndist = 1; 
   int Nveg = 1; 
 
//BOC
 
   for ( i = 0; i < tspmd.cdi_array[par.rank]; i++ )
   {
      energy  = prcp[i].energy;
      cell    = prcp[i].cell;
      veg_var = prcp[i].veg_var;
      snow    = prcp[i].snow;
  
      for ( veg = 0; veg <= Nveg; ++veg )
//      for ( veg = 0; veg < Nveg; ++veg )
      {
         for ( band = 0; band < *snowband; ++band )
         {
            for ( dist = 0; dist < Ndist; ++dist )
            {
               for ( lidx = 0; lidx < *nlayer; ++lidx )
               {
                  model_state[i].moist[dist][veg][band][lidx] = 
                     cell[dist][veg][band].layer[lidx].moist; 
                  model_state[i].ice[dist][veg][band][lidx]   = 
                    cell[dist][veg][band].layer[lidx].ice; 
               }
  
               if ( veg < Nveg )
               {
                  model_state[i].Wdew[dist][veg][band]  = 
                     veg_var[dist][veg][band].Wdew;
               }
            }

            for ( nidx = 0; nidx < *nnode; ++nidx )
            {
               model_state[i].T[veg][band][nidx] = energy[veg][band].T[nidx];
            }
  
            model_state[i].last_snow[veg][band]   = snow[veg][band].last_snow;
            model_state[i].swq[veg][band]         = snow[veg][band].swq;
            model_state[i].surf_temp[veg][band]   = snow[veg][band].surf_temp;
            model_state[i].pack_temp[veg][band]   = snow[veg][band].pack_temp;
            model_state[i].density[veg][band]     = snow[veg][band].density;
            model_state[i].snow_canopy[veg][band] = snow[veg][band].snow_canopy;
            model_state[i].pack_water[veg][band]  = snow[veg][band].pack_water;
            model_state[i].surf_water[veg][band]  = snow[veg][band].surf_water;

            model_state[i].MELTING[veg][band]     = snow[veg][band].MELTING;
            model_state[i].coverage[veg][band]    = snow[veg][band].coverage;
            model_state[i].coldcontent[veg][band] = snow[veg][band].coldcontent;
  
         }
      }
   }
#if ( ( defined SPMD ) && ( ! defined OPENDAP ) )
   if ( par.npes > 1 )
   {
      MPI_Gatherv(model_state,tspmd.cdi_array[par.rank],
                  par.MPI_MS_STRUCT,model_state,tspmd.cdi_array,tspmd.cdispls,
                  par.MPI_MS_STRUCT,0,MPI_COMM_WORLD);
   }
#endif  
//EOC
}
