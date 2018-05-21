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

void surface_fluxes(int                  rec,
		    char                 overstory,
		    int                  band,
		    int                  veg_class,
		    int                  iveg,
		    int                  Nveg,
		    int                  Ndist,
		    int                  Nbands,
		    int                  Nlayers,
		    int                  Nnode,
		    int                  dp,
		    //		    float               mu,
		    float               ice0,
		    float               moist,
		    float               surf_atten,
		    float               height,
		    float               displacement,
		    float               roughness,
		    float               ref_height,
		    float               bare_albedo,
		    float              *aero_resist,
		    float              *baseflow_wet,
		    float              *baseflow_dry,
		    float              *runoff_wet,
		    float              *runoff_dry,
		    float              *out_prec,
		    float              *wind,
		    float              *Le,
		    float              *Ls,
		    float              *Melt,
		    float              *inflow_wet,
		    float              *inflow_dry,
		    float              *snow_inflow,
		    float              *gauge_correction,
		    float               *root,
		    atmos_data_struct   *atmos,
		    soil_con_struct     *soil_con,
		    dmy_struct          dmy,
		    energy_bal_struct   *energy,
		    snow_data_struct    *snow,
		    layer_data_struct   *layer_wet,
		    layer_data_struct   *layer_dry,
		    veg_var_struct      *veg_var_wet,
		    veg_var_struct      *veg_var_dry,
		    int                 full_energy_flag,
		    int                 frozen_soil_flag,
		    int                 grnd_flux_flag,
		    int                 quick_flux_flag)
/**********************************************************************
	surface_fluxes	Keith Cherkauer		February 29, 2000

  Formerally a part of full_energy.c this routine computes all surface
  fluxes, and solves the snow accumulation and ablation algorithm.
  Solutions are for the current snow band and vegetation type (these
  are defined in full_energy before the routine is called).

  modifications:
  02-07-03 fixed indexing problem for sub-daily snow model within
           daily water balance VIC: hour (now hidx) is incremented
           by 1 rather than the sub-daily time step, so the atmospheric
           forcing data is now properly indexed.                KAC
  04-23-03 Indexing fix sent SNOW_STEP to calc_surf_energy_bal rather
           than the model time step, meaning that without snow the
           evaporation was computed for SNOW_STEP hours rather than a
           full day.  This was fixed by introducing step_inc to
           index the arrays, while step_dt keeps track of the correct
           time step.                                            KAC
  07-May-04 Fixed initialization of canopyevap to initialize for every
            value of dist, rather than just dist 0.             TJB

**********************************************************************/
{
  int                    dist;
  int                    lidx;
  /* int                    hour; */
  /* int                    endhour; */
  /* int			endhidx; */
  /* int			hidx; */
  /* int 			step_inc; */
#if TIMESTEPSECS
  float                  step_dt;
#else
  int                    step_dt;
#endif
  int                    N_steps;
  float                 Tsurf;
  float                 air_temp;
  float                 T0;
  float                 step_rad;
  float                 step_net_short;
  float                 ppt[2];
  float                 rainfall[2];
  float                 step_ppt[2];
  float                 step_snow_energy;
  float                 step_out_prec;
  float                 step_out_short;
  float                 step_Evap;
  float                 step_melt;
  float                 step_prec[2];
  float                 store_throughfall[2];
  float                 store_melt;
  float                 store_vapor_flux;
  float                 store_canopy_vapor_flux;
  float                 store_canopyevap[2];
  float                 store_layerevap[2][MAX_LAYERS];
  float                 store_ppt[2];
  float                 store_shortwave;
  float                 store_longwave;
  float                 store_sensible;
  float                 store_latent;
  float                 store_grnd_flux;
  float                 store_deltaH;
  float                 store_advection; 
  float                 store_deltaCC; 
  float                 store_snow_flux; 
  float                 store_refreeze_energy; 
  float                 store_albedo;
  layer_data_struct      step_layer[2][MAX_LAYERS];
  energy_bal_struct      snow_energy;
  energy_bal_struct      bare_energy;
  snow_data_struct       step_snow;
  veg_var_struct         snow_veg_var[2];
  veg_var_struct         bare_veg_var[2];
  float mu = 1; 

  /***********************************************************************
    Set temporary variables - preserves original values until iterations
    are completed
  ***********************************************************************/
  energy->advection       = 0; 
  energy->deltaCC         = 0; 
  energy->snow_flux       = 0; 
  energy->refreeze_energy = 0; 
  snow_energy             = (*energy);
  bare_energy             = (*energy);
  snow_veg_var[WET]       = (*veg_var_wet);
  snow_veg_var[DRY]       = (*veg_var_dry);
  bare_veg_var[WET]       = (*veg_var_wet);
  bare_veg_var[DRY]       = (*veg_var_dry);
  step_snow               = (*snow);
  for(lidx=0;lidx<Nlayers;lidx++) {
    step_layer[WET][lidx] = layer_wet[lidx];
    step_layer[DRY][lidx] = layer_dry[lidx];
  }

  /********************************
    Set-up sub-time step controls
    (May eventually want to set this up so that it is also true 
    if frozen soils are present)
  ********************************/
  /**
  if(snow->swq > 0 || snow->snow_canopy > 0 || atmos->snowflag[NR]) {
    hidx      = 0;
    endhidx   = hidx + NF;
    step_inc  = 1;
    step_dt   = options.SNOW_STEP;
  }
  else {
    hidx      = NR;
    endhidx   = hidx + dmy.dt/3600.0;
    step_inc = step_dt = dmy.dt;
    step_dt   = dmy.dt/3600.0;
    }*/
  step_dt = dmy.dt;

  /*******************************************
    Initialize sub-model time step variables
  *******************************************/

  for(dist = 0; dist < Ndist; dist++) {
    store_throughfall[dist] = 0.;
    store_canopyevap[dist]  = 0.;
    for(lidx = 0; lidx < Nlayers; lidx++) {
      store_layerevap[dist][lidx] = 0.;
    }
  }
  store_canopy_vapor_flux = 0;
  store_vapor_flux        = 0;
  store_ppt[WET]          = 0;
  store_ppt[DRY]          = 0;
  store_melt              = 0;
  step_prec[DRY]          = 0;
  (*snow_inflow)          = 0;
  store_shortwave         = 0;
  store_longwave          = 0;
  store_sensible          = 0;
  store_latent            = 0;
  store_grnd_flux         = 0;
  store_deltaH            = 0;
  store_advection         = 0; 
  store_deltaCC           = 0; 
  store_snow_flux         = 0; 
  store_refreeze_energy   = 0; 
  store_albedo            = 0;
  N_steps                 = 0;
      
  /*************************
    Compute surface fluxes 
  *************************/

  //  do {

    /*********************************************************
      Solve for snow interception, accumulation and ablation
    *********************************************************/

    air_temp = atmos->air_temp;
    //    printf("here... %d %d %f %f %f\n",rec,band,atmos->prec,mu,soil_con->Pfactor[band]);
    step_prec[WET] = atmos->prec / mu;
    //    printf("TEMP,, %f %f \n",atmos->air_temp,soil_con->Tfactor[band]);
    for(dist = 0; dist < Ndist; dist ++) {
      snow_veg_var[dist].canopyevap = 0;
      bare_veg_var[dist].canopyevap = 0;
      for(lidx = 0; lidx < Nlayers; lidx ++) 
	step_layer[dist][lidx].evap = 0;
    }
    step_snow.canopy_vapor_flux = 0;
    step_snow.vapor_flux = 0;

    if ( grnd_flux_flag) T0 = energy->T[0];
    else T0 = air_temp;
    
    step_melt = solve_snow( rec,&(step_snow), step_layer[WET], step_layer[DRY],
			    &(snow_veg_var[WET]), &(snow_veg_var[DRY]), 
			    dmy.month, dmy.day_in_year,
			    &(snow_energy), soil_con, overstory, step_dt, 
			    veg_class, iveg, Nveg, band, dmy.hour, 
			    Nnode, atmos->shortwave, 
			    atmos->longwave, air_temp, 
			    step_prec[WET], atmos->density, 
			    atmos->vp, atmos->vpd, 
			    atmos->pressure, mu, roughness, 
			    displacement, ref_height, surf_atten, 
			    moist, ice0, dp, bare_albedo, 
			    rainfall, &step_out_prec, Le, Ls, aero_resist, 
			    wind, &step_net_short, &step_out_short, &step_rad, 
			    &step_Evap, &step_snow_energy, snow_inflow, 
			    step_ppt, gauge_correction, root, Nlayers,
			    full_energy_flag,frozen_soil_flag);

    /**********************************************************
      Solve Energy Balance Components for Ground Free of Snow
    **********************************************************/
    Tsurf = calc_surf_energy_bal(iveg, Nlayers, Nveg, step_dt, 
				 Nnode, veg_class, band, dmy.hour, 
				 ice0, moist, dp, surf_atten, T0, 
				 atmos->shortwave, 
				 snow_energy.longwave, air_temp, (*Le), 
				 (*Ls), mu, displacement, roughness, 
				 ref_height, step_snow_energy, 
				 step_snow.vapor_flux,  bare_albedo, 
				 aero_resist, wind, 
				 rainfall, step_ppt, root, 
				 atmos, &bare_veg_var[WET], 
				 &bare_veg_var[DRY], &bare_energy, 
				 &(step_snow), step_layer[WET], 
				 step_layer[DRY], soil_con, &(dmy),
				 full_energy_flag, quick_flux_flag,
				 grnd_flux_flag,frozen_soil_flag);  

    /**************************************
      Store sub-model time step variables 
    **************************************/

    for(dist = 0; dist < Ndist; dist++) {

      if(iveg < Nveg) {
	if(step_snow.snow) {
	  store_throughfall[dist] += snow_veg_var[dist].throughfall;
	  store_canopyevap[dist]  += snow_veg_var[dist].canopyevap;
	  bare_veg_var[dist].Wdew  = snow_veg_var[dist].Wdew;
	}
	else {
	  store_throughfall[dist] += bare_veg_var[dist].throughfall;
	  store_canopyevap[dist]  += bare_veg_var[dist].canopyevap;
	  snow_veg_var[dist].Wdew  = bare_veg_var[dist].Wdew;
	}
      }

      for(lidx = 0; lidx < Nlayers; lidx++)
	store_layerevap[dist][lidx] += step_layer[dist][lidx].evap;

      store_ppt[dist] += step_ppt[dist];
    }

    if(iveg < Nveg) 
      store_canopy_vapor_flux += step_snow.canopy_vapor_flux;
    store_vapor_flux += step_snow.vapor_flux;
      
    store_melt  += step_melt;
    out_prec[0] += step_out_prec * mu;
   
    if(overstory){
      store_shortwave += step_snow.coverage * snow_energy.shortwave 
	* surf_atten + (1. - step_snow.coverage) * bare_energy.shortwave; 
      //   printf("sw %d %f %f %f %f %f\n",rec,step_snow.coverage,snow_energy.shortwave,
      //	     surf_atten, bare_energy.shortwave,atmos->shortwave);
    }
    else
      store_shortwave += step_snow.coverage * snow_energy.shortwave 
	+ (1. - step_snow.coverage) * bare_energy.shortwave; 
    store_longwave  += step_snow.coverage * snow_energy.longwave 
      + (1. - step_snow.coverage) * bare_energy.longwave; 
    store_latent    += step_snow.coverage * snow_energy.latent 
      + (1. - step_snow.coverage) * bare_energy.latent; 
    store_sensible  += step_snow.coverage * snow_energy.sensible 
      + (1. - step_snow.coverage) * bare_energy.sensible; 
    store_grnd_flux += bare_energy.grnd_flux; 
    store_deltaH    += bare_energy.deltaH; 
    store_albedo    += step_snow.coverage * snow_energy.albedo
      + (1. - step_snow.coverage) * bare_energy.albedo;
    if (step_snow.snow) {
      store_advection       += snow_energy.advection; 
      store_deltaCC         += snow_energy.deltaCC; 
      store_snow_flux       += bare_energy.snow_flux; 
      store_refreeze_energy += snow_energy.refreeze_energy; 
    }
    
    /* increment time step */
    //    N_steps ++;
    //hour += step_dt;

    //  } while (hour < endhour);

  /************************************************
    Store snow variables for sub-model time steps 
  ************************************************/

  (*snow) = step_snow;
  snow->vapor_flux = store_vapor_flux;
  snow->canopy_vapor_flux = store_canopy_vapor_flux;
  (*Melt) = store_melt;
  snow->melt = store_melt; 
  for(dist = 0; dist < 2; dist++) ppt[dist] = store_ppt[dist];

  /******************************************************
    Store energy flux averages for sub-model time steps 
  ******************************************************/

  (*energy) = bare_energy;
  if(overstory)
    energy->shortwave = store_shortwave ; 
  else
    energy->shortwave = store_shortwave ;
  energy->longwave    = store_longwave ;
  energy->latent      = store_latent ;
  energy->sensible    = store_sensible;
  energy->grnd_flux   = store_grnd_flux;
  energy->deltaH      = store_deltaH ;
  energy->albedo      = store_albedo ;
  if (snow->snow) {
    energy->advection       = store_advection ;
    energy->deltaCC         = store_deltaCC ;
    energy->snow_flux       = store_snow_flux;
    energy->refreeze_energy = store_refreeze_energy;
  }
	  
  /**********************************************************
    Store vegetation variable sums for sub-model time steps 
  **********************************************************/

  if(iveg < Nveg) {
    veg_var_wet->throughfall = store_throughfall[WET];
    veg_var_dry->throughfall = store_throughfall[DRY];
    veg_var_wet->canopyevap  = store_canopyevap[WET];
    veg_var_dry->canopyevap  = store_canopyevap[DRY];
    if(snow->snow) {
      veg_var_wet->Wdew        = snow_veg_var[WET].Wdew;
      veg_var_dry->Wdew        = snow_veg_var[DRY].Wdew;
    }
    else {
      veg_var_wet->Wdew        = bare_veg_var[WET].Wdew;
      veg_var_dry->Wdew        = bare_veg_var[DRY].Wdew;
    }
  }

  /**********************************************************
    Store soil layer variables for sub-model time steps 
  **********************************************************/

  for(lidx=0;lidx<Nlayers;lidx++) {
    layer_wet[lidx]      = step_layer[WET][lidx];
    layer_dry[lidx]      = step_layer[DRY][lidx];
    layer_wet[lidx].evap = store_layerevap[WET][lidx];
    layer_dry[lidx].evap = store_layerevap[DRY][lidx];
  }


  /********************************************************
    Compute Runoff, Baseflow, and Soil Moisture Transport
  ********************************************************/

  (*inflow_wet) = ppt[WET];
  (*inflow_dry) = ppt[DRY];

  runoff(layer_wet, layer_dry, energy, soil_con, runoff_wet, runoff_dry, 
	 baseflow_wet, baseflow_dry, ppt, mu, dmy.dt, Nnode, 
         band, iveg, Nlayers,full_energy_flag,frozen_soil_flag);
	    
  /**  for(lidx=0;lidx<Nlayers;lidx++) {
    printf("SOILM %d %d %f \n",rec,lidx,layer_wet[lidx].moist);
    }*/
}

