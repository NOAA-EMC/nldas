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
#include "ftn.h"
float calc_moist_by_depth(float, layer_data_struct *, float *, int );
float calc_max_moist_by_depth(float, float *, float *, int);
//BOP
//
// !ROUTINE: vic_run.c
//
// !DESCRIPTION: 
//  Executes VIC model runs on each tile
//
// !INTERFACE:
void FTN(vic_run)(int *nlayer,
                  int *nnode,
                  int *snowband,
                  int *day, 
                  int *day_in_year,
                  int *hour, 
                  int *month,
                  int *year,
                  int *full_energy_flag,
                  int *frozen_soil_flag,
                  int *grnd_flux_flag,
                  int *quick_flux_flag,
                  int *outputflag)
{
  //EOP
  extern par_struct par;
  extern ctile_spmd tspmd;
  extern dmy_struct dmy;
  extern atmos_data_struct *atmos;
  extern soil_con_struct *soil_con; 
  extern veg_con_struct **veg_con; 
  extern dist_prcp_struct *prcp; 
  extern out_data_struct *outdata1; 

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy; 
  veg_var_struct       ***veg_var;

  int i, veg, band, dist,index; 
  int lidx,nidx;
  float tmp_rainf = 0.0;
  float tmp_depth = 0.0;
  float tmp_evap = 0.0;
  float t_veg = 0.0;
  float tmp_ice = 0.0;
  float tmp_moist = 0.0;
  float surf_temp = 0.0;
  float rad_temp = 0.0;
  float albedo = 0.0;
  float swe = 0.0;
  float tmp_moist_sum = 0.0; 
  float tmp_max_moist_sum = 0.0;
  int Ndist = 1; 
  int Nveg = 1; 
  //BOC
  dmy.day = *day; 
  dmy.day_in_year = *day_in_year; 
  dmy.hour = *hour;
  dmy.month = *month; 
  dmy.year = *year; 
  for ( i = 0; i < tspmd.cdi_array[par.rank]; i++)
  {
//<Justin's frozen soil fix>
    if ( *frozen_soil_flag == 1 )
    {
       soil_con[i].FS_ACTIVE = 1;
    }
//</Justin's frozen soil fix>
    full_energy(i,&atmos[i], &soil_con[i], veg_con[i], &prcp[i], 
                dmy, *nlayer, *nnode, *snowband, *full_energy_flag, 
                *frozen_soil_flag,
                *grnd_flux_flag, *quick_flux_flag); 
  }
#if 1
  for ( i = 0; i < tspmd.cdi_array[par.rank]; i++ )
  {
    outdata1[i].count += 1;
    energy  = prcp[i].energy;
    cell    = prcp[i].cell;
    veg_var = prcp[i].veg_var;
    snow    = prcp[i].snow;

    for ( veg = 0; veg <= Nveg; veg++ )
    {
      for ( dist = 0; dist < Ndist; dist++ )
      {
        for ( band = 0; band < *snowband; band++ )
        {
          for ( lidx = 0; lidx < *nlayer; lidx++ )
          {
            cell[dist][veg][band].layer[lidx].T = 0.;
            for ( nidx = 0; nidx < *nnode; nidx++ )
            {
              if ( soil_con[i].layer_node_fract[lidx][nidx] > 0 )
              {
                cell[dist][veg][band].layer[lidx].T += 
                  energy[veg][band].T[nidx] * 
                  soil_con[i].layer_node_fract[lidx][nidx] * 
                  soil_con[i].dz_node[nidx];
              }
            }
            cell[dist][veg][band].layer[lidx].T /= soil_con[i].depth[lidx];
          }
        }
      }
    }

    outdata1[i].acond = 0.0;
    for ( veg = 0; veg <= Nveg; veg++ )
    {
      if ( cell[WET][veg][0].aero_resist[0] != 0.0 )
         outdata1[i].acond += 1.0/cell[WET][veg][0].aero_resist[0];
      for ( band = 0; band < *snowband; band++ )
      {
        outdata1[i].swnet += energy[veg][band].shortwave;
        outdata1[i].lwnet += energy[veg][band].longwave;
        outdata1[i].qle   += energy[veg][band].latent; 
        outdata1[i].qh    += energy[veg][band].sensible; 
        outdata1[i].qg    += energy[veg][band].grnd_flux +
                             energy[veg][band].deltaH; 
        outdata1[i].qfz   += energy[veg][band].refreeze_energy/3.33e5; 
      }
/* Justin: rainfall and snowfall should not be updated within the veg loop */
//      tmp_rainf =calc_rainonly(atmos[i].air_temp,atmos[i].prec,
//                               1.0);
//      outdata1[i].rainf += tmp_rainf; 
//      outdata1[i].snowf +=(atmos[i].prec-tmp_rainf);
/* Justin: end of changes */
    }
/* Justin: rainfall and snowfall updates moved outside of veg loop */
    tmp_rainf =calc_rainonly(atmos[i].air_temp,atmos[i].prec,1.0);
      outdata1[i].rainf += tmp_rainf/dmy.dt; 
      outdata1[i].snowf +=( (atmos[i].prec-tmp_rainf)/dmy.dt );
/* Justin: end of changes */
/* start of changes -- JS */
    outdata1[i].soilt[0] = 0.0;
    outdata1[i].soilt[1] = 0.0;
    outdata1[i].soilt[2] = 0.0;
/* end of changes -- JS */
    outdata1[i].soilwet   = 0.0;
    outdata1[i].rootmoist = 0.0;
    for ( veg = 0; veg <= Nveg; veg++ )
    {
      for (dist = 0; dist < Ndist; dist++ )
      {
        for ( band = 0; band < *snowband; band++ )
        {
          outdata1[i].soilt[0] += cell[dist][veg][band].layer[0].T;
          outdata1[i].soilt[1] += cell[dist][veg][band].layer[1].T;
          outdata1[i].soilt[2] += cell[dist][veg][band].layer[2].T;
          tmp_depth = 0.0;
          for ( index = 0; index < *nlayer; index++ )
          {
            tmp_depth +=soil_con[i].depth[index]; 
          }
          tmp_moist_sum = calc_moist_by_depth(tmp_depth,
                          cell[dist][veg][band].layer,soil_con[i].depth, *nlayer);
          tmp_max_moist_sum = calc_max_moist_by_depth(tmp_depth,
                              soil_con[i].porosity,soil_con[i].depth, *nlayer);
          outdata1[i].soilwet +=(tmp_max_moist_sum==0.0?0.0:
                                 tmp_moist_sum/tmp_max_moist_sum); 
          tmp_depth = 0.0; 
          if ( veg < Nveg )
          {
            for ( index = 0; index < 2; index ++)
            {
              tmp_depth +=veg_con[i][veg].zone_depth[index];
            }
            tmp_moist_sum = calc_moist_by_depth(tmp_depth, 
                            cell[dist][veg][band].layer, soil_con[i].depth, 
                            *nlayer);
            tmp_max_moist_sum = calc_max_moist_by_depth(tmp_depth, 
                                soil_con[i].porosity, soil_con[i].depth, 
                                *nlayer);
          }
          else 
          {
            tmp_moist_sum = 0.0; 
          }
          outdata1[i].rootmoist += tmp_moist_sum;
        }
      }
      outdata1[i].soilwet /=2; 
    }

    outdata1[i].snowt = 0.0;
    for ( veg = 0; veg < Nveg; veg++ )
    {
      tmp_evap  = 0.0; 
      surf_temp = 0.0;
      rad_temp  = 0.0;
      albedo    = 0.0; 
      t_veg     = 0.0;
      tmp_moist = 0.0;
      tmp_ice   = 0.0;
      swe       = 0.0; 
      for ( band = 0; band < *snowband; band++ )
      {
        tmp_evap  += (snow[veg][band].vapor_flux*1000 +
                      snow[veg][band].canopy_vapor_flux*1000);
        surf_temp += energy[veg][band].T[0]+KELVIN;
        if ( snow[veg][band].swq > 0 ) 
        {
          rad_temp += (snow[veg][band].surf_temp+KELVIN);
        }
        else
        {
          rad_temp += (energy[veg][band].T[0]+KELVIN);
        }
        albedo += energy[veg][band].albedo;
        swe    += snow[veg][band].swq;
        outdata1[i].qsm += snow[veg][band].melt/(dmy.dt); 
        for ( dist = 0; dist < Ndist; dist++ )
        {
          for ( index = 0; index < *nlayer; index++ )
          {
            t_veg += cell[dist][veg][band].layer[index].evap;
            if ( veg < Nveg ) 
            {
              outdata1[i].tveg = t_veg/(dmy.dt); 
            }
            else 
            {
              outdata1[i].esoil = t_veg/(dmy.dt); 
            }

            tmp_evap  += cell[dist][veg][band].layer[index].evap;
            tmp_moist  = cell[dist][veg][band].layer[index].moist; 
            tmp_ice    = cell[dist][veg][band].layer[index].ice; 
            tmp_moist -=tmp_ice; 
            outdata1[i].moist[index] = tmp_moist; 
/* Justin: qs and qsb updates should not be within layer loop */
//            outdata1[i].qs+=cell[dist][veg][band].runoff/(dmy.dt);
//            outdata1[i].qsb+=cell[dist][veg][band].baseflow/(dmy.dt);
/* Justin: end of changes */
          }
/* Justin: qs and qsb updates moved outside of layer loop */
          outdata1[i].qs  += cell[dist][veg][band].runoff/(dmy.dt);
          outdata1[i].qsb += cell[dist][veg][band].baseflow/(dmy.dt);
/* Justin: end of changes */
          tmp_evap += veg_var[dist][veg][band].canopyevap;
        }
        outdata1[i].snowt += snow[veg][band].surf_temp; 
      }
      outdata1[i].avgsurft = surf_temp;
      outdata1[i].radt     = rad_temp;
      outdata1[i].evap    += (tmp_evap/dmy.dt); 
      outdata1[i].albedo   = albedo; 
      outdata1[i].swe      = swe; 
      if ( *outputflag )
      {
        outdata1[i].moist_prev = outdata1[i].moist[0] +
                                 outdata1[i].moist[1] +
                                 outdata1[i].moist[2]; 
        outdata1[i].swe_prev = outdata1[i].swe; 
      }
    }
  }
#endif
  //EOC
}

//BOP
//
// !ROUTINE: calc_moist_by_depth
//
// !DESCRIPTION:
//  This functions calculates the soil moisture depth from the
//  surface to a given depth. It assumes that the soil moisture
//  is evenly distributed over a soil layer (which it is, I think).
//
// !INTERFACE:
float calc_moist_by_depth(float the_depth,
                           layer_data_struct *layer,
                           float *depth, int Nlayer)
//EOP
{
  int    index;
  float moist;
  float depth_total;
  
  
  //if(the_depth <= 0.0)
  //   printf("calc_moist_by_depth: depth must be > 0.0\n");  
  //BOC
  moist = 0.0;
  depth_total = 0.0;
  //-------------------------------------------------------------------
  // loop through the soil layers from top to bottom
  //-------------------------------------------------------------------
  for(index=0;index<Nlayer;index++) {
    //-------------------------------------------------------------------
    // calculate the depth to the bottom of the current layer 
    //-------------------------------------------------------------------
     depth_total += depth[index];
     //-------------------------------------------------------------------
     // compare with the depth in question 
    //-------------------------------------------------------------------
     if(the_depth > depth_total) {
       //-------------------------------------------------------------------
       // the depth is past the current layer 
       //-------------------------------------------------------------------
        moist += layer[index].moist;
     }
        
     else {
       //-------------------------------------------------------------------
       // the depth is within the current layer 
       //-------------------------------------------------------------------
        moist += layer[index].moist * (the_depth - 
			   (depth_total - depth[index]))/depth[index];  
        break; 
     }  
    
  }

  return(moist);
  //EOC
}

//BOP
//
// !ROUTINE: calc_max_moist_by_depth
// 
// !DESCRIPTION:
//  This functions calculates the max soil moisture depth from the
//  surface to a given depth. It assumes that the soil moisture
//  is evenly distributed over a soil layer (which it is, I think).
//  Max allowable soil moisture = depth * porosity.
// 
// !INTERFACE:
float calc_max_moist_by_depth(float the_depth,
                           float *porosity,
                           float *depth, int Nlayer) 
     //EOP
{
  int    index;
  float moist;
  float depth_total;

  
  //  if(the_depth <= 0.0)
  //  printf("calc_max_moist_by_depth: depth must be > 0.0\n");  
  //BOC
  moist = 0.0;
  depth_total = 0.0;
  //---------------------------------------------------------------
  // loop through the soil layers from top to bottom 
  //---------------------------------------------------------------
  for(index=0;index<Nlayer;index++) {
    //---------------------------------------------------------------
    // calcaulate the depth to the bottom of the current layer 
    //---------------------------------------------------------------
     depth_total += depth[index];
     //---------------------------------------------------------------
     // compare with the depth in question 
     //---------------------------------------------------------------
     if(the_depth > depth_total) {
       moist += porosity[index] * 1000. * depth[index];
     }
        
     else {
       moist += porosity[index] * 1000. * 
	 depth[index] * (the_depth - (depth_total 
				      - depth[index]))/depth[index];  
       break; 
     }  
     
  }
  return(moist);
  //EOC
}







