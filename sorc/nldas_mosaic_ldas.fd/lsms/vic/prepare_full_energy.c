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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vicNl.h"

void prepare_full_energy(int               iveg,
                         int               Nveg,
                         int               Nnodes,
                         dist_prcp_struct *prcp,
                         soil_con_struct  *soil_con,
                         float            *moist,
                         float            *ice0,
                         int               Nlayer,
                         int               snowband,
                         int               frozen_soil_flag)
{
/*******************************************************************
  prepare_full_energy.c      Keith Cherkauer       January 20, 2000

  This subroutine returns the soil thermal properties, moisture 
  and ice contents for the top two layers for use with the QUICK_FLUX
  ground heat flux solution.

  Modifications:
  01-20-00 split into separate file, formerly at the end of 
           full_energy.c                                      KAC

*******************************************************************/


  int                i, band;
  layer_data_struct *layer;

  layer = (layer_data_struct *)lis_calloc(Nlayer,
                                          sizeof(layer_data_struct),
                                          "prepare_full_energy");

  for ( band = 0; band < snowband; band++)
  {

    /* Compute average soil moisture values for distributed precipitation */

    for ( i = 0; i < Nlayer; i++) 
    {
      layer[i] = find_average_layer(&(prcp->cell[WET][iveg][band].layer[i]),
                                    &(prcp->cell[DRY][iveg][band].layer[i]),
                                    soil_con->depth[i], 1);
    }
    
    /* Compute top soil layer moisture content (mm/mm) */

    (*moist) = layer[0].moist / ( soil_con->depth[0] * 1000. );

    /* Compute top soil layer ice content (mm/mm) */

    if ( frozen_soil_flag && soil_con->FS_ACTIVE)
    {
      if ( (prcp->energy[iveg][band].T[0] + prcp->energy[iveg][band].T[1]) / 2.
           < 0. )
      {
        (*ice0) = (*moist) - 
                  maximum_unfrozen_water((prcp->energy[iveg][band].T[0] + 
                                          prcp->energy[iveg][band].T[1]) / 2.,
                                         soil_con->max_moist[0] / 
                                          (soil_con->depth[0] * 1000.),
                                         soil_con->bubble[0], 
                                         soil_con->expt[0]);
        if ( (*ice0) < 0. )
        {
          (*ice0) = 0.;
        }
      }
      else
      {
        (*ice0) = 0.;
      }
    }
    else
    {
      (*ice0) = 0.;
    }
    //printf("SOIL PARAM..%f %f\n",soil_con->quartz[0],soil_con->quartz[1]);
    /** Compute Soil Thermal Properties **/
    compute_soil_layer_thermal_properties(layer,soil_con->depth,
                                          soil_con->bulk_density,
                                          soil_con->soil_density,
                                          soil_con->quartz,Nlayer);
    
    /** Save Thermal Conductivities for Energy Balance **/
    prcp->energy[iveg][band].kappa[0] = layer[0].kappa; 
    prcp->energy[iveg][band].Cs[0]    = layer[0].Cs; 
    prcp->energy[iveg][band].kappa[1] = layer[1].kappa; 
    prcp->energy[iveg][band].Cs[1]    = layer[1].Cs; 
  }

  free((void *)layer);

}
