#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

/** Start of changes - NGMOD **/
#include "vicNl_ldas.h"
/** End of changes - NGMOD **/

static char vcid[] = "$Id: put_data.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

/** Start of changes -- JS **/
double calc_moist_by_depth(double, layer_data_struct *, double *);
double calc_max_moist_by_depth(double, double *, double *);
/** End of changes -- JS **/


void put_data_ldas(dist_prcp_struct  *prcp,
		   	atmos_data_struct *atmos,
		   	veg_con_struct    *veg_con,
/** Start of changes - NGMOD **/
		   /* outfiles_struct   *outfiles, */
		   	ldas_outfiles_struct   *outfiles,
/** End of changes - NGMOD **/
		   	double            *depth,
		   	double            *dz,
		   	double             dp,
		   	double            *AreaFract,
/** Start of changes - JS **/
		   	double            *porosity,
/** End of changes - JS **/
		   	dmy_struct        *dmy,
		   	int                rec,
		   	int                dt,
		   	int                Nnodes,
		   	int                skipyear,
/** Start of changes -- JS **/
			double               MAX_SNOW_TEMP,
			double               MIN_RAIN_TEMP)
/** End of changes -- JS **/
/**********************************************************************
	put_data.c	Dag Lohmann		January 1996

  This routine converts data units, and stores finalized values
  in an array for later output to the output files.

  modifications:
  06-24-98  modified for new distributed presipitation data structures KAC
  01-20-00 modified to deal with simplified frozen soil moisture layers
           and frost depth / thaw depth accounting                 KAC
  03-08-00 modified to eliminate extra lines for storing bare
           soil variables.                                         KAC
  10-09-00 stored incoming shortwave in output_data struct			JS
  											

**********************************************************************/
{
  extern veg_lib_struct  *veg_lib;
  extern option_struct    options;
#if LINK_DEBUG
  extern debug_struct     debug;
#endif

#if OPTIMIZE
  static int              prtdt;
  static double           runoff;
  static double           baseflow;
#endif

  out_data_struct        *out_data;

  int                     veg;
  int                     index;
  int                     Ndist;
  int                     dist;
  int                     band;
  int                     Nbands;
  double                  Cv;
  double                  mu;
  double                  tmp_evap;
  double                  tmp_moist;
/** Start of changes -- JS **/
  double                  tmp_liquid;
  double                  tmp_max_moist;
  double                  tmp_moist_sum;
  double                  tmp_ice_sum;
  double                  tmp_liquid_sum;
  double                  tmp_max_moist_sum;
  double                  tmp_depth;
  double                  tmp_soil_temp;
  /** End of changes -- JS **/
  double                  tmp_ice;
  double                  rad_temp;
  double                  surf_temp;
  double                  inflow;
  double                  outflow;
  double                  storage;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;
  Nbands = options.SNOW_BAND;

  /** Initialize Output Data Array **/
  out_data = (out_data_struct *) calloc(1,sizeof(out_data_struct)); 
  
  out_data->surf_temp      = 0.;
  out_data->snow_canopy[0] = 0.0;
  
  /*************************************************
    Store Output for General Types
    *************************************************/
    
  /** record precipitation [mm] **/
  out_data->prec           = atmos->prec[NR];
  
  /** record wind speed [m/s] **/
  out_data->wind           = atmos->wind[NR];
  
  /** record air temperature [C] **/
  out_data->air_temp       = atmos->air_temp[NR];
  
  /** record relative humidity [%] **/
  out_data->rel_humid      = 100.*atmos->vp[NR]/(atmos->vp[NR]+atmos->vpd[NR]);
  
/** Start of changes -- JS **/
  /** record incoming shortwave radiation [W/m2] **/
  out_data->shortwave = atmos->shortwave[NR];
  
  /** record rainfall [mm] **/
  out_data->rainfall = calc_rainonly(atmos->air_temp[NR], atmos->prec[NR], 
				MAX_SNOW_TEMP, MIN_RAIN_TEMP, 1.0);
				
  /** record snowfall [mm] **/
  out_data->snowfall = atmos->prec[NR] - out_data->rainfall;
/** End of changes -- JS **/

 
  /*************************************************
    Store Output for Precipitation Distribution Type
    *************************************************/

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;
  
/*printf("init: %f %f %f %f %f %f\n", out_data->evap_veg, out_data->evap_bare, out_data->evap_canop, out_data->sub_snow, out_data->sub_canop, out_data->evap);*/
  /****************************************
    Store Output for all Vegetation Types
  ****************************************/
  for ( veg = 0 ; veg <= veg_con[0].vegetat_type_num ; veg++) {

    if ( veg < veg_con[0].vegetat_type_num ) 
      Cv = veg_con[veg].Cv;
    else
      Cv = (1.0 - veg_con[0].Cv_sum);

    if ( Cv > 0 ) {

/** Start of changes -- JS **/
      /** record LAI [-] **/
      if ( veg < veg_con[0].vegetat_type_num ) 
        out_data->lai += veg_lib[veg_con[veg].veg_class].LAI[dmy->month-1] * Cv; 
/** End of changes -- JS **/

      /*******************************************************
        Compute Average Variables from Wet and Dry Fractions
      *******************************************************/
      for ( dist = 0; dist < Ndist; dist++ ) {
	if(dist==0) 
	  mu = prcp[0].mu[veg];
	else 
	  mu = 1. - prcp[0].mu[veg];
	      
	for(band=0;band<Nbands;band++) {
	  if(AreaFract[band] > 0.) {
	    
	  /*********************************
            Record Water Balance Variables 
	  *********************************/
#if !OPTIMIZE
	    
	    /** evaporation components **/
	    tmp_evap = 0.0;
	    
	    /** record veg/soil evaporation [mm] **/
	    for(index=0;index<options.Nlayer;index++)
	      tmp_evap += cell[dist][veg][band].layer[index].evap;
	    if ( veg < veg_con[0].vegetat_type_num )
	      out_data->evap_veg += tmp_evap * Cv * mu * AreaFract[band];
	    else
	      out_data->evap_bare += tmp_evap * Cv * mu * AreaFract[band];
	      
	    /** record snow pack sublimation [mm] **/
	    tmp_evap += snow[veg][band].vapor_flux * 1000.;
	    out_data->sub_snow += snow[veg][band].vapor_flux * 1000. 
	      * Cv * mu * AreaFract[band]; 
	      
	    /** record canopy snow sublimation [mm] **/
	    if ( veg <= veg_con[0].vegetat_type_num ) {
	      tmp_evap += snow[veg][band].canopy_vapor_flux * 1000.;
	      out_data->sub_canop += snow[veg][band].canopy_vapor_flux 
		* 1000. * Cv * mu * AreaFract[band]; 
	    }
	    
	    /** record canopy evaporation [mm] **/
	    if ( veg < veg_con[0].vegetat_type_num ) {
	      tmp_evap += veg_var[dist][veg][band].canopyevap;
	      out_data->evap_canop += veg_var[dist][veg][band].canopyevap 
		* Cv * mu * AreaFract[band]; 
	    }
	    
	    /** record total evaporation [mm] **/
	    out_data->evap += tmp_evap * Cv * mu * AreaFract[band]; 
/*printf("%d %d %d: %f %f %f %f %f %f\n", veg, dist, band, out_data->evap_veg, out_data->evap_bare, out_data->evap_canop, out_data->sub_snow, out_data->sub_canop, out_data->evap);*/

#endif
	  
	    /** record runoff [mm] **/
	    out_data->runoff   += cell[dist][veg][band].runoff 
	      * Cv * mu * AreaFract[band];
	    
	    /** record baseflow [mm] **/
	    out_data->baseflow += cell[dist][veg][band].baseflow 
	      * Cv * mu * AreaFract[band]; 
	    
#if ! OPTIMIZE
	    
	    /** record inflow [mm] **/
	    if ( veg < veg_con[0].vegetat_type_num ) 
	      out_data->inflow += (cell[dist][veg][band].inflow 
				   + veg_var[dist][veg][band].canopyevap) 
		* Cv * mu * AreaFract[band];
	    else 
	      out_data->inflow += (cell[dist][veg][band].inflow) 
		* Cv * mu * AreaFract[band];
	    
	    /** record canopy interception [mm] **/
	    if ( veg < veg_con[0].vegetat_type_num ) 
	      out_data->Wdew += veg_var[dist][veg][band].Wdew 
		* Cv * mu * AreaFract[band];
	  
	    /** record aerodynamic resistance [s/m] **/
	    out_data->aero_resist += cell[WET][veg][0].aero_resist[0] 
	      * Cv * mu * AreaFract[band];
	  
/** Start of changes -- JS **/
	    /** record aerodynamic conductance [m/s] (=1/resistance) **/
	    out_data->aero_conduct += (1./cell[WET][veg][0].aero_resist[0]) 
	      * Cv * mu * AreaFract[band];
/** End of changes -- JS **/
	      
	    /** recored layer moistures [mm] and temperatures [C] **/
/** Start of changes -- JS **/
            tmp_depth = 0.;
	    for(index=0;index<options.Nlayer;index++) {

              /* get the total moisture, liquid moisture, ice content and max allowable moisture*/
	      tmp_moist     = cell[dist][veg][band].layer[index].moist;
	      tmp_ice       = cell[dist][veg][band].layer[index].ice;
	      tmp_liquid    = tmp_moist - tmp_ice;
	      tmp_max_moist = porosity[index] * depth[index] * 1000.;
	      	      
	      /* convert to fractions if needed */
	      if(options.MOISTFRACT) {
		tmp_moist  /= depth[index] * 1000.;
		tmp_liquid /= depth[index] * 1000.;
		tmp_ice    /= depth[index] * 1000.;
	      }
	      	      
	      /* update the totals for the whole grid cell */
	      out_data->moist[index]  += tmp_moist * Cv * mu * AreaFract[band];
	      out_data->ice[index]    += tmp_ice * Cv * mu * AreaFract[band];
	      out_data->liquid[index] += tmp_liquid * Cv * mu * AreaFract[band];
	      
              /** record soil layer temperatures [C] **/
              tmp_soil_temp = cell[dist][veg][band].layer[index].T;
	      out_data->soil_temp[index] += tmp_soil_temp * Cv * mu * AreaFract[band];

	      /* calculate total column depth */
	      tmp_depth         += depth[index];

	    }
	    
            /** record the total column soil moisture and wetness [mm] **/
            tmp_moist_sum = calc_moist_by_depth(tmp_depth, cell[dist][veg][band].layer, depth);
            tmp_max_moist_sum = calc_max_moist_by_depth(tmp_depth, porosity, depth);
            tmp_moist_sum *= Cv * mu * AreaFract[band];
            tmp_max_moist_sum *= Cv * mu * AreaFract[band];
	    if(options.MOISTFRACT) {
              tmp_moist_sum /= tmp_depth * 1000.;
              tmp_max_moist_sum /= tmp_depth * 1000.;
	    }
            out_data->moist_total += tmp_moist_sum;
            out_data->moist_max_total += tmp_max_moist_sum;
            out_data->wetness_total += (tmp_max_moist_sum==0.0?0.0:tmp_moist_sum/tmp_max_moist_sum) * Cv * mu * AreaFract[band];

            /** record the root zone soil moisture and wetness [mm] **/
            tmp_depth = 0.;
            if ( veg < veg_con[0].vegetat_type_num ) {
              for(index=0;index<options.ROOT_ZONES;index++)
                 tmp_depth += veg_con[veg].zone_depth[index];
              tmp_moist_sum = calc_moist_by_depth(tmp_depth, cell[dist][veg][band].layer, depth);
              tmp_max_moist_sum = calc_max_moist_by_depth(tmp_depth, porosity, depth);
              tmp_moist_sum *= Cv * mu * AreaFract[band];
              tmp_max_moist_sum *= Cv * mu * AreaFract[band];
	      if(options.MOISTFRACT) {
                tmp_moist_sum /= tmp_depth * 1000.;
                tmp_max_moist_sum /= tmp_depth * 1000.;
	      }
            }
            else
               tmp_moist_sum = 0.;
            out_data->moist_root += tmp_moist_sum;
            out_data->wetness_root += (tmp_max_moist_sum==0.0?0.0:tmp_moist_sum/tmp_max_moist_sum) * Cv * mu * AreaFract[band];

            /** record the 1m top layer soil moisture [mm] **/
            tmp_depth = 1.;
            tmp_moist_sum = calc_moist_by_depth(tmp_depth, cell[dist][veg][band].layer, depth);
            tmp_moist_sum *= Cv * mu * AreaFract[band];
	    if(options.MOISTFRACT) {
              tmp_moist_sum /= tmp_depth * 1000.;
	    }
            out_data->moist_1m += tmp_moist_sum;
            
/** End of changes -- JS **/

#endif
	  }
	}
      }

#if !OPTIMIZE
      for(band=0;band<Nbands;band++) {
	if(AreaFract[band] > 0.) {

	  /**********************************
	    Record Frozen Soil Variables
	  **********************************/

	  /** record freezing and thawing front depths [cm] **/
	  if(options.FROZEN_SOIL) {
	    for(index = 0; index < MAX_FRONTS; index++) {
	      if(energy[veg][band].fdepth[index] != MISSING)
		out_data->fdepth[index] += energy[veg][band].fdepth[index] 
		  * Cv * 100. * AreaFract[band];
	      if(energy[veg][band].tdepth[index] != MISSING)
		out_data->tdepth[index] += energy[veg][band].tdepth[index] 
		  * Cv * 100. * AreaFract[band];
	    }
	  }

	  /**********************************
            Record Energy Balance Variables
	  **********************************/

	  /** calculate the surface radiative temperature [K] **/
	  if(snow[veg][band].swq>0) 
	    rad_temp = snow[veg][band].surf_temp + KELVIN;
	  else
	    rad_temp = energy[veg][band].T[0] + KELVIN;
	  
	  /** calculate soil surface temperature [C] **/
	  surf_temp = energy[veg][band].T[0];
	  
	  /** record net shortwave radiation [W/m2] **/
	  out_data->net_short += energy[veg][band].shortwave
	    * Cv * AreaFract[band];
	  	    
	  /** record net longwave radiation [W/m2] **/
	  out_data->net_long  += energy[veg][band].longwave
	    * Cv * AreaFract[band];
	  
	  /** Start of changes -- JS **/
	  /* what is this incoming longwave? - why can't we just use the atmos value of longwave? unless it is the same? TODO: check this out */
	  /** End of changes -- JS **/
	  /** record incoming longwave radiation [W/m2] **/
	  out_data->in_long  += ((energy[veg][band].longwave + STEFAN_B 
				  * (rad_temp) * (rad_temp)
				  * (rad_temp) * (rad_temp))
				 * Cv * AreaFract[band]);
	  
	  /** record albedo [-] **/
	  out_data->albedo    += energy[veg][band].albedo
	    * Cv * AreaFract[band];
	  
	  /** record latent heat flux [W/m2] **/
	  out_data->latent    -= energy[veg][band].latent
	    * Cv * AreaFract[band];
	  
	  /** record sensible heat flux [W/m2] **/
	  out_data->sensible  -= energy[veg][band].sensible
	    * Cv * AreaFract[band];
	  
	  /** record ground heat flux (+ heat storage) [W/m2] **/
	  out_data->grnd_flux -= (energy[veg][band].grnd_flux
				  + energy[veg][band].deltaH)
	    * Cv * AreaFract[band];
	  
	  /** record heat storage [W/m2] **/
	  out_data->deltaH    -= energy[veg][band].deltaH
	    * Cv * AreaFract[band];
	  
	  /** record energy balance error [W/m2] **/
	  out_data->energy_error += energy[veg][band].error
	    * Cv * AreaFract[band];
	  
	  /** record radiative effective temperature [K], 
	      emissivities set = 1.0  **/
	  out_data->rad_temp += ((rad_temp) * (rad_temp) 
				 * (rad_temp) * (rad_temp)) 
	    * Cv * AreaFract[band];
	  
	  /** record mean surface temperature [C]  **/
	  out_data->surf_temp += surf_temp * Cv * AreaFract[band];
	  
	  /*****************************
	    Record Snow Pack Variables 
	  *****************************/
	  
	  /** record snow water equivalence (convert from m to mm) [mm] **/
	  out_data->swq[0]         
	    += snow[veg][band].swq * 1000. * Cv * AreaFract[band];
	  
	  /** record snowpack depth (convert from m to cm) [cm] **/
	  out_data->snow_depth[0]  
	    += snow[veg][band].depth * 100. * Cv * AreaFract[band];
	  
/** Start of changes -- JS **/
          /** record snow melt [mm] **/
          out_data->snow_melt[0]
	    += snow[veg][band].melt * Cv * AreaFract[band];
	    
          /** record snow pack temperature [C] **/
          out_data->snow_pack_temp[0]
	    += snow[veg][band].pack_temp * Cv * AreaFract[band];

          /** record snow albedo [-] **/
          out_data->snow_albedo[0]
	    += snow[veg][band].albedo * Cv * AreaFract[band];
/** End of changes -- JS **/
	  
	  /** record canopy intercepted snow (convert from m to mm) [mm] **/
	  if ( veg < veg_con[0].vegetat_type_num )
#if LDAS_OUTPUT
	    out_data->swq[0] 
#else
	      out_data->snow_canopy[0] 
#endif
	      += (snow[veg][band].snow_canopy) 
	      * Cv * 1000. * AreaFract[band];
	  
	  /** record snow cover fraction [-] **/
	  out_data->coverage[0]    
	    += snow[veg][band].coverage * Cv * AreaFract[band];
	  
	  /** record snowpack cold content [W/m2] **/
	  out_data->deltaCC[0]           
	    += energy[veg][band].deltaCC * Cv * AreaFract[band];
	  
	  /** record snowpack advection [W/m2] **/
	  out_data->advection[0]         
	    += energy[veg][band].advection * Cv * AreaFract[band];
	  
	  /** record snow energy flux [W/m2] **/
	  out_data->snow_flux[0]         
	    += energy[veg][band].snow_flux * Cv * AreaFract[band];
	  
	  /** record refreeze energy [W/m2] **/
	  out_data->refreeze_energy[0]   
	    += energy[veg][band].refreeze_energy 
	    * Cv * AreaFract[band];
	  
/** Start of changes -- JS **/
	  /** record snow energy fluxes [W/m2] **/
	  out_data->snow_fluxes[0]         
	    += (energy[veg][band].advection 
	      - energy[veg][band].deltaCC
	      - energy[veg][band].snow_flux
	      + energy[veg][band].refreeze_energy) * Cv * AreaFract[band];
/** End of changes -- JS **/
	  
	  /*****************************
	    Record Elevation Band Snow Pack Variables 
	  *****************************/
	  if(options.PRT_SNOW_BAND) {
	    
	    /** record band snow water equivalent (convert from m to mm) [mm] **/
	    out_data->swq[band+1]         
	      += snow[veg][band].swq * Cv  * 1000.;
	    
	    /** record band snowpack depth (convert from m to cm) [cm] **/
	    out_data->snow_depth[band+1]  
	      += snow[veg][band].depth * Cv * 100.;
	    
/** Start of changes -- JS **/
            /** record snow melt [mm] **/
            out_data->snow_melt[band+1]
	      += snow[veg][band].melt * Cv;
	      
            /** record snow pack temperature [C] **/
            out_data->snow_pack_temp[band+1]
	      += snow[veg][band].pack_temp * Cv;
	      
            /** record snow albedo [-] **/
            out_data->snow_albedo[band+1]
	      += snow[veg][band].albedo * Cv;
/** End of changes -- JS **/

	    /** record band canopy intercepted snow (convert from m to mm) [mm] **/
	    if ( veg < veg_con[0].vegetat_type_num )
#if LDAS_OUTPUT
	      out_data->swq[band+1]
#else
		out_data->snow_canopy[band+1]  
#endif
		+= (snow[veg][band].snow_canopy) * Cv * 1000.;
	    
	    /** record band snow coverage [-] **/
	    out_data->coverage[band+1]    
	      += snow[veg][band].coverage * Cv;
	    
	    /** record band cold content [W/m2] **/
	    out_data->deltaCC[band+1]           
	      += energy[veg][band].deltaCC * Cv;
	    
	    /** record band advection [W/m2] **/
	    out_data->advection[band+1]         
	      += energy[veg][band].advection * Cv;
	    
	    /** record band snow flux [W/m2] **/
	    out_data->snow_flux[band+1]         
	      += energy[veg][band].snow_flux * Cv;
	    
	    /** record band refreeze energy [W/m2] **/
	    out_data->refreeze_energy[band+1]   
	      += energy[veg][band].refreeze_energy * Cv;
	    
/** Start of changes -- JS **/
	    /** record snow energy fluxes [W/m2] **/
	    out_data->snow_fluxes[band+1]         
	      += (energy[veg][band].advection 
	        - energy[veg][band].deltaCC
	        - energy[veg][band].snow_flux
	        + energy[veg][band].refreeze_energy) * Cv;
/** End of changes -- JS **/
	  }
	}
      }

#endif /* not OPTIMIZE */

    }
  }
  
#if !OPTIMIZE

  /** record radiative temperature [K] **/
  out_data->rad_temp = pow(out_data->rad_temp,0.25);

  /** record net radiation [W/m2] **/
  out_data->r_net    = out_data->net_short + out_data->net_long;

  /********************
    Check Water Balance 
    ********************/
  inflow  = out_data->prec;
  outflow = out_data->evap+out_data->runoff+out_data->baseflow;
  storage = 0.;
  for(index=0;index<options.Nlayer;index++)
    if(options.MOISTFRACT)
      storage += (out_data->moist[index] + out_data->ice[index]) 
	* depth[index] * 1000;
    else
      storage += out_data->moist[index] + out_data->ice[index];
  storage += out_data->swq[0] + out_data->snow_canopy[0] + out_data->Wdew;
  calc_water_balance_error(rec,inflow,outflow,storage);
  if(options.FULL_ENERGY)
    calc_energy_balance_error(rec, out_data->net_short + out_data->net_long,
			      out_data->latent, out_data->sensible,
			      out_data->grnd_flux, out_data->advection[0] 
			      - out_data->deltaCC[0] - out_data->snow_flux[0]
			      + out_data->refreeze_energy[0]);

/*  printf("end: %f %f %f %f %f %f\n", out_data->evap_veg, out_data->evap_bare, out_data->evap_canop, out_data->sub_snow, out_data->sub_canop, out_data->evap);*/
/*************
    Write Data
  *************/

  if(rec >= skipyear)
/** Start of changes - NGMOD **/
    /* write_data(out_data, outfiles, dmy, dt); */
    write_data_ldas(out_data, outfiles, dmy, dt);
/** End of changes - NGMOD **/

#else

  if ( rec == 0 ) prtdt = 0;
  if ( prtdt == 0 ) {
    runoff = out_data->runoff;
    baseflow = out_data->baseflow;
    prtdt ++;
  }
  else {
    runoff += out_data->runoff;
    baseflow += out_data->baseflow;
    prtdt ++;
  }
  if ( prtdt == 24 / dt ) {
    out_data->runoff = runoff;
    out_data->baseflow = baseflow;
/** Start of changes - NGMOD **/
    /* write_data(out_data, outfiles, dmy, dt); */
    write_data_ldas(out_data, outfiles, dmy, dt);
/** End of changes - NGMOD **/
    prtdt = 0;
  } 

#endif

  free((char *)out_data); 

}

/* Start of changes -- JS **/
double calc_moist_by_depth(double the_depth,
                           layer_data_struct *layer,
                           double *depth) {
/*************************************************************
  This functions calculates the soil moisture depth from the
  surface to a given depth. It assumes that the soil moisture
  is evenly distributed over a soil layer (which it is, I think).
**************************************************************/

  int    index;
  double moist;
  double depth_total;
  
  extern option_struct    options;
  
  if(the_depth <= 0.0)
     vicerror("calc_moist_by_depth: depth must be > 0.0");  
  
  moist = 0.0;
  depth_total = 0.0;
  
  /* loop through the soil layers from top to bottom */
  for(index=0;index<options.Nlayer;index++) {
     
     /* calcaulate the depth to the bottom of the current layer */
     depth_total += depth[index];
     
     /* compare with the depth in question */
     if(the_depth > depth_total) {
        /* the depth is past the current layer */
        moist += layer[index].moist;
     }
        
     else {
        /* the depth is within the current layer */
        moist += layer[index].moist * (the_depth - (depth_total - depth[index]))/depth[index];  
        break; 
     }  
    
  }

  return(moist);

}


double calc_max_moist_by_depth(double the_depth,
                           double *porosity,
                           double *depth) {
/*************************************************************
  This functions calculates the max soil moisture depth from the
  surface to a given depth. It assumes that the soil moisture
  is evenly distributed over a soil layer (which it is, I think).
  Max allowable soil moisture = depth * porosity.
**************************************************************/

  int    index;
  double moist;
  double depth_total;

  extern option_struct    options;
  
  if(the_depth <= 0.0)
     vicerror("calc_max_moist_by_depth: depth must be > 0.0");  

  moist = 0.0;
  depth_total = 0.0;
  
  /* loop through the soil layers from top to bottom */
  for(index=0;index<options.Nlayer;index++) {
     
     /* calcaulate the depth to the bottom of the current layer */
     depth_total += depth[index];
     
     /* compare with the depth in question */
     if(the_depth > depth_total) {
        /* the depth is past the current layer */
        moist += porosity[index] * 1000. * depth[index];
     }
        
     else {
        /* the depth is within the current layer */
        moist += porosity[index] * 1000. * depth[index] * (the_depth - (depth_total - depth[index]))/depth[index];  
        break; 
     }  
    
  }

  return(moist);

}/* End of changes -- JS **/

