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
#include "ftn.h"
#include "vicNl.h"
void  FTN(initialize_model_state)(int *nlayer, 
                                  int *snowband,
                                  int *Nnodes,
                                  int *quick_flux,
                                  int *grnd_flux,
                                  int *frozen_soil,
                                  float *init_surf_temp)
/**********************************************************************
  initialize_model_state       Keith Cherkauer	    April 17, 2000

  This routine initializes the model state (energy balance, water balance,
  and snow components).  If a state file is provided to the model than its
  contents are checked to see if it agrees with the current simulation
  set-up, if so it is used to initialize the model state.  If no state
  file is provided the model initializes all variables with defaults and
  the user should expect to throw out the beginning of the simulation 
  period as model start-up.

  UNITS: (m, s, kg, C, moisture in mm) unless otherwise specified

  Modifications:
  4-17-00 Modified from initialize_energy_bal.c and initialize_snow.c
          to provide a single controlling routine for initializing the
          model state.
**********************************************************************/
{
  extern par_struct par; 
  extern ctile_spmd tspmd; 
  extern soil_con_struct *soil_con; 
  extern dist_prcp_struct *prcp; 
  /* extern veg_con_struct **veg_con;  */
  /* extern atmos_data_struct *atmos; */

  int      i, veg, index;
  int      nidx,lidx,dry;
  int      Nveg,Ndist; 
  int      band;
  float   sum, dp, Ltotal;
  float   moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];
  float   ice[MAX_VEG][MAX_BANDS][MAX_LAYERS];
  
  float Zsum,tmpadj,tmpdp;
  char FIRST_VEG;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;
  /* air temp.. needs to be properly assigned.. */
  float surf_temp; 
  Ndist = 1; 

  for ( i = 0; i < tspmd.cdi_array[par.rank]; i++)
  {

     /* Justin: FS_ACTIVE needs to be set if frozen_soil flag is true */
     if ( *frozen_soil == 1 )
     {
        soil_con[i].FS_ACTIVE = 1;
     }
     /* Justin: end of changes */
     /* surf_temp = atmos[i].air_temp; */
     surf_temp = *init_surf_temp;
     cell    = prcp[i].cell;
     veg_var = prcp[i].veg_var;
     snow    = prcp[i].snow;
     energy  = prcp[i].energy;
    //    Nveg   = veg_con[0][0].vegetat_type_num; 
     Nveg   = 1;

     dp = soil_con[i].dp;
     Ltotal = 0;

     for ( index = 0; index < *nlayer; index++)
     {
        Ltotal += soil_con[i].depth[index];
     }

     FIRST_VEG = TRUE;

     /********************************************
      Initialize all snow pack variables 
      - some may be reset if state file present
      ********************************************/
     initialize_snow(snow, Nveg, i, *snowband);

     /********************************************
      Initialize all soil layer variables 
      - some may be reset if state file present
      ********************************************/
      
     initialize_soil(cell[WET], &(soil_con[i]), Nveg, *snowband, *nlayer);

     /********************************************
      Initialize all vegetation variables 
      - some may be reset if state file present
      ********************************************/

     initialize_veg(veg_var[WET],Nveg, *snowband);


  
     /************************************************************************
      CASE 2: Initialize soil if using quick heat flux, and no initial
      soil properties file given
      ************************************************************************/
    
     if ( *quick_flux )
     {
        for ( veg = 0 ; veg <= Nveg ; veg++)
        {
           for ( band = 0; band < *snowband; band++ )
           {
    
              /* Initialize soil node temperatures and thicknesses */
    
              soil_con[i].dz_node[0] = soil_con[i].depth[0];
              soil_con[i].dz_node[1] = soil_con[i].depth[0];
              soil_con[i].dz_node[2] = 2. * (dp - 1.5 * soil_con[i].depth[0]);
              energy[veg][band].T[0] = surf_temp;
              energy[veg][band].T[1] = surf_temp;
              energy[veg][band].T[2] = soil_con[i].avg_temp;
    
              for ( lidx = 0; lidx < *nlayer; lidx++) {
                moist[veg][band][lidx] = soil_con[i].init_moist[lidx];
                ice[veg][band][lidx] = 0.;
              }
             
           }
        }
      }
      else if ( !(*quick_flux) )
      {
         for ( veg = 0; veg <= Nveg; veg++)
         {
            for ( band = 0; band < *snowband; band++)
            {
    
               /* Initialize soil node temperatures and thicknesses 
                  Nodes set at surface, the depth of the first layer,
                  twice the depth of the first layer, and at the
                  damping depth.  Extra nodes are placed equal distance
                  between the damping depth and twice the depth of the
                  first layer. */

               energy[veg][band].T[0] = surf_temp;
               soil_con[i].dz_node[0] = soil_con[i].depth[0];
               soil_con[i].dz_node[1] = soil_con[i].depth[0];
               soil_con[i].dz_node[2] = soil_con[i].depth[0];
               energy[veg][band].T[*Nnodes-1] = soil_con[i].avg_temp;
               energy[veg][band].T[1] = exp_interp(soil_con[i].depth[0], 
                                                   0.,
                                                   dp, 
                                                   surf_temp, 
                                                   soil_con[i].avg_temp);
               energy[veg][band].T[2] = exp_interp(2. * soil_con[i].depth[0], 
                                                   0., 
                                                   dp,
                                                   surf_temp, 
                                                   soil_con[i].avg_temp);
    
               Zsum   = 2. * soil_con[0].depth[0];
               tmpdp  = dp - soil_con[0].depth[0] * 2.5;
               tmpadj = 3.5;

               for ( index = 3; index < *Nnodes-1; index++)
               {
                  if ( veg == 0 && band == 0 )
                  {
                     soil_con[i].dz_node[index]=tmpdp/(((float)*Nnodes-tmpadj));
                  }
                  Zsum += (soil_con[i].dz_node[index] +
                           soil_con[i].dz_node[index-1])/2.;
                  energy[veg][band].T[index] = exp_interp(Zsum,
                                                          0.,
                                                          soil_con[0].dp,
                                                          surf_temp,
                                                          soil_con[0].avg_temp);
               }
               if ( veg == 0 && band == 0)
               {
                  soil_con[i].dz_node[*Nnodes-1] = 
                     (dp - Zsum - soil_con[i].dz_node[*Nnodes-2] / 2. ) * 2.;
                  Zsum += (soil_con[i].dz_node[*Nnodes-2] +
                           soil_con[i].dz_node[*Nnodes-1])/2.;
                  if ( (int)(Zsum*1000+0.5) != (int)(dp*1000+0.5) )
                  {
                     printf("Sum of thermal node thicknesses (%f) in initialize_model_state do not equal dp (%f), check initialization procedure\n",Zsum,dp);
                     //nrerror(ErrStr);
                  }
               }
    
               for ( lidx = 0; lidx < *nlayer; lidx++) { 
                 moist[veg][band][lidx] = soil_con[i].init_moist[lidx];
                 ice[veg][band][lidx] = 0.;
               }
            }
         }
      }
      else
      {
         for ( veg = 0 ; veg <= Nveg ; veg++)
         {
            for ( band = 0; band < *snowband; band++)
            {
               for ( index = 0; index < *nlayer; index++)
               {
                  soil_con->dz_node[index] = 1.;
               }
            }
         }
      }

      /* dz_node[Nnodes] is accessed later despite not being set. This can
      cause run-time errors on some platforms. Therefore, set it to zero */
      soil_con->dz_node[*Nnodes]=0.0;

      if ( *grnd_flux )
      {
         for ( veg = 0 ; veg <= Nveg ; veg++)
         {
            for( band = 0; band < *snowband; band++ )
            {

               /** Set soil properties for all soil nodes **/
               if ( FIRST_VEG )
               {
                  FIRST_VEG = FALSE;
                  set_node_parameters(soil_con[i].dz_node,
                                      soil_con[i].max_moist_node,
                                      soil_con[i].expt_node,
                                      soil_con[i].bubble_node,
                                      soil_con[i].alpha,
                                      soil_con[i].beta,
                                      soil_con[i].gamma,
                                      soil_con[i].depth,
                                      soil_con[i].max_moist,
                                      soil_con[i].expt, 
                                      soil_con[i].bubble,
                                      soil_con[i].quartz, 
                                      soil_con[i].layer_node_fract,
                                      *Nnodes, 
                                      *nlayer, 
                                      soil_con[i].FS_ACTIVE);

                  sum = soil_con[i].dz_node[0]/2. + 
                        soil_con[i].dz_node[*Nnodes-1]/2.;
                  for ( nidx = 1; nidx < *Nnodes-1; nidx++)
                  {
                     sum += soil_con[i].dz_node[nidx];
                  }
               }

               /* set soil moisture properties for all soil thermal nodes */
               distribute_node_moisture_properties(energy[veg][band].moist,
                                                   energy[veg][band].ice,
                                                   energy[veg][band].kappa_node,
                                                   energy[veg][band].Cs_node,
                                                   soil_con[i].dz_node,
                                                   energy[veg][band].T,
                                                   soil_con[i].max_moist_node,
                                                   soil_con[i].expt_node,
                                                   soil_con[i].bubble_node,
                                                   moist[veg][band],
                                                   soil_con[i].depth,
                                                   soil_con[i].soil_density,
                                                   soil_con[i].bulk_density,
                                                   soil_con[i].quartz,
                                                   *Nnodes,
                                                   *nlayer,
                                                   soil_con[i].FS_ACTIVE,
                                                   *frozen_soil);

               /* initialize layer moistures and ice contents */
               for ( dry = 0; dry < Ndist; dry++)
               {
                  for ( lidx = 0; lidx < *nlayer; lidx++)
                  {
                     cell[dry][veg][band].layer[lidx].moist =
                        moist[veg][band][lidx];
                     cell[dry][veg][band].layer[lidx].ice = 
                        ice[veg][band][lidx];
                  }
                  if(soil_con->FS_ACTIVE && *frozen_soil)
                    estimate_layer_ice_content(cell[dry][veg][band].layer,
                                             soil_con[i].dz_node,
                                             energy[veg][band].T,
                                             soil_con[i].max_moist_node,
                                             soil_con[i].expt_node,
                                             soil_con[i].bubble_node,
                                             soil_con[i].depth,
                                             soil_con[i].max_moist,
                                             soil_con[i].expt,
                                             soil_con[i].bubble,
                                             soil_con[i].bulk_density,
                                             soil_con[i].soil_density,
                                             soil_con[i].quartz, 
                                             soil_con[i].layer_node_fract,
                                             *Nnodes,
                                             *nlayer, 
                                             soil_con[i].FS_ACTIVE,
                                             *frozen_soil);

               }

               /* Find freezing and thawing front depths */
               if ( !(*quick_flux) && soil_con[i].FS_ACTIVE)
               { 
                  find_0_degree_fronts(&energy[veg][band],
                                       soil_con[i].dz_node,
                                       energy[veg][band].T,
                                       *Nnodes);
               }

            }
         }

      }
   }
}
