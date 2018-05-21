#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

/** Start of changes - NGMOD **/
#include "vicNl_ldas.h"
/** End of changes - NGMOD **/


static char vcid[] = "$Id: initialize_model_state.c,v 4.2 2000/05/16 21:57:54 vicadmin Exp vicadmin $";

void initialize_model_state_ldas(dist_prcp_struct    *prcp,
				 dmy_struct           dmy,
				 double               surf_temp,
				 global_param_struct *global_param,
				 infiles_struct       infiles, 
				 int                  cellnum,
				 int                  Nveg,
				 int                  Nnodes,
				 int                  Ndist,
				 soil_con_struct     *soil_con,
				 veg_con_struct      *veg_con,
				 int                 *dry_time,
/** Start of changes -- JS **/
                         int                *still_storm
/** End of changes -- JS **/
				 )
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

 Modified to read in WET and DRY variable for distributed precip
                                               June 2000, ngmod
**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif
#if QUICK_FS
  extern double temps[];
#endif

  char     tmpstr[MAXSTRING];
  char     ErrStr[MAXSTRING];
  char     FIRST_VEG;
  int      i, j, ii, veg, index;
  int      nidx, lidx;
  int      tmpint;
  int      dry;
  int      band;
  int      zindex;
  double   sum, Lsum, Zsum, dp, Ltotal;
  double   tmpdp, tmpadj;
  double  *kappa, *Cs, *M;
  double   moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];
  double   ice[MAX_VEG][MAX_BANDS][MAX_LAYERS];
  double   unfrozen, frozen;
  double **layer_ice;
  double **layer_tmp;
  double  *EMPTY;
#if QUICK_FS
  double   Aufwc, Bufwc;
#endif
  char    *EMPTY_C;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;
  
  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;

  dp = soil_con->dp;
  Ltotal = 0;
  for(index=0;index<options.Nlayer;index++) Ltotal += soil_con->depth[index];
  FIRST_VEG = TRUE;
  
  /********************************************
    Initialize all snow pack variables 
    - some may be reset if state file present
  ********************************************/
  initialize_snow(snow, Nveg, infiles.init_snow, cellnum);

  /********************************************
    Initialize all soil layer variables 
    - some may be reset if state file present
  ********************************************/

  initialize_soil(cell[WET], soil_con, Nveg);

  /********************************************
    Initialize all vegetation variables 
    - some may be reset if state file present
  ********************************************/

  initialize_veg(veg_var[WET], veg_con, global_param);

#if QUICK_FS
  if(options.FROZEN_SOIL) {

    /***********************************************************
      Prepare table of maximum unfrozen water content values
      - This linearizes the equation for maximum unfrozen water
        content, reducing computation time for the frozen soil
        model.
    ***********************************************************/

    for(lidx=0;lidx<options.Nlayer;lidx++) { 
      for(ii=0;ii<QUICK_FS_TEMPS;ii++) {
	Aufwc = maximum_unfrozen_water(temps[ii], 1.0, 
				       soil_con->bubble[lidx], 
				       soil_con->expt[lidx]);
	Bufwc = maximum_unfrozen_water(temps[ii+1], 1.0, 
				       soil_con->bubble[lidx], 
				       soil_con->expt[lidx]);
	soil_con->ufwc_table_layer[lidx][ii][0] 
	  = linear_interp(0., temps[ii], temps[ii+1], Aufwc, Bufwc);
	soil_con->ufwc_table_layer[lidx][ii][1] 
	  = (Bufwc - Aufwc) / (temps[ii+1] - temps[ii]);
      }
    }
  }  
#endif

  /************************************************************************
    CASE 1: Not using quick ground heat flux, and initial conditions files 
    provided
  ************************************************************************/

  if(options.INIT_STATE) {
    /** Start of changes - NGMOD **/
    /* code change to deal with wet and dry fractions */
    read_initial_model_state_ldas(infiles.statefile, prcp, global_param,  
			     Nveg, options.SNOW_BAND, cellnum, soil_con,
			     veg_con, dry_time
/** Start of changes -- JS **/
                         , still_storm
/** End of changes -- JS **/
			      );

    for( veg = 0; veg <= Nveg; veg++ ) 
      for( band = 0; band < options.SNOW_BAND; band++ )
	for( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	  if(options.DIST_PRCP == TRUE){
	    moist[veg][band][lidx] = prcp->mu[veg]*cell[WET][veg][band].layer[lidx].moist
	      +(1-prcp->mu[veg])*cell[DRY][veg][band].layer[lidx].moist;
	    ice[veg][band][lidx] = prcp->mu[veg]*cell[WET][veg][band].layer[lidx].ice
	      +(1-prcp->mu[veg])*cell[DRY][veg][band].layer[lidx].ice;
	  }
	  else {
	    moist[veg][band][lidx] = cell[0][veg][band].layer[lidx].moist;
	    ice[veg][band][lidx] = cell[0][veg][band].layer[lidx].ice; 
	  }
	}
    /** End of changes - NGMOD **/
  }
  
  /************************************************************************
    CASE 2: Initialize soil if using quick heat flux, and no initial
    soil properties file given
  ************************************************************************/
    
  else if(options.QUICK_FLUX) {

    for ( veg = 0 ; veg <= Nveg ; veg++) {
      for ( band = 0; band < options.SNOW_BAND; band++ ) {

	/* Initialize soil node temperatures and thicknesses */

	soil_con->dz_node[0] = soil_con->depth[0];
	soil_con->dz_node[1] = soil_con->depth[0];
	soil_con->dz_node[2] = 2. * (dp - 1.5 * soil_con->depth[0]);
	energy[veg][band].T[0] = surf_temp;
	energy[veg][band].T[1] = surf_temp;
	energy[veg][band].T[2] = soil_con->avg_temp;

	for( veg = 0; veg <= Nveg; veg++ ) 
	  for( band = 0; band < options.SNOW_BAND; band++ )
	    for(lidx=0;lidx<options.Nlayer;lidx++) {
	      moist[veg][band][lidx] = soil_con->init_moist[lidx];
	      ice[veg][band][lidx] = 0.;
	    }

      }
    }
  }

  /*****************************************************************
    CASE 3: Initialize Energy Balance Variables if not using quick
    ground heat flux, and no Initial Condition File Given 
  *****************************************************************/
  else if(!options.QUICK_FLUX) {
    for ( veg = 0 ; veg <= Nveg ; veg++) {
      for(band=0;band<options.SNOW_BAND;band++) {

	/* Initialize soil node temperatures and thicknesses 
	 Nodes set at surface, the depth of the first layer,
	 twice the depth of the first layer, and at the
	 damping depth.  Extra nodes are placed equal distance
	 between the damping depth and twice the depth of the
	 first layer. */

	energy[veg][band].T[0] = surf_temp;
	soil_con->dz_node[0] = soil_con->depth[0];
	soil_con->dz_node[1] = soil_con->depth[0];
	soil_con->dz_node[2] = soil_con->depth[0];
	energy[veg][band].T[Nnodes-1] = soil_con->avg_temp;
	energy[veg][band].T[1] = exp_interp(soil_con->depth[0], 0., dp, 
					    surf_temp, soil_con->avg_temp);
	energy[veg][band].T[2] = exp_interp(2. * soil_con->depth[0], 0., dp, 
					    surf_temp, soil_con->avg_temp);
	
        Zsum   = 2. * soil_con[0].depth[0];
        tmpdp  = dp - soil_con[0].depth[0] * 2.5;
        tmpadj = 3.5;
        for(index=3;index<Nnodes-1;index++) {
          if(veg==0 && band==0)
	    soil_con->dz_node[index] = tmpdp/(((double)Nnodes-tmpadj));
          Zsum += (soil_con->dz_node[index]
                   +soil_con->dz_node[index-1])/2.;
          energy[veg][band].T[index] = exp_interp(Zsum,0.,soil_con[0].dp,
                                                  surf_temp,
                                                  soil_con[0].avg_temp);
        }
	if(veg==0 && band==0) {
	  soil_con->dz_node[Nnodes-1] = (dp - Zsum 
					 - soil_con->dz_node[Nnodes-2] 
					 / 2. ) * 2.;
	  Zsum += (soil_con->dz_node[Nnodes-2]
		   +soil_con->dz_node[Nnodes-1])/2.;
	  if((int)(Zsum*1000+0.5) != (int)(dp*1000+0.5)) {
	    sprintf(ErrStr,"Sum of thermal node thicknesses (%f) in initialize_model_state do not equal dp (%f), check initialization procedure",Zsum,dp);
	    nrerror(ErrStr);
	  }
        }

/** Start of changes -- JS **/
/* Why are we looping over veg and band within a veg and band loop? */
/*	for( veg = 0; veg <= Nveg; veg++ ) 
	  for( band = 0; band < options.SNOW_BAND; band++ )*/
/** End of changes -- JS **/
	    for(lidx=0;lidx<options.Nlayer;lidx++) {
	      moist[veg][band][lidx] = soil_con->init_moist[lidx];
	      ice[veg][band][lidx] = 0.;
	    }

      }
    }
  }

  /*********************************
    CASE 4: Unknown option
  *********************************/
  else {
    for ( veg = 0 ; veg <= Nveg ; veg++) {
      for(band=0;band<options.SNOW_BAND;band++) {
	for(index=0;index<options.Nlayer;index++) {
	  soil_con->dz_node[index] = 1.;
	}
      }
    }
  }

  /******************************************
    Initialize soil thermal node properties 
  ******************************************/

/** Start of changes -- JS **/
/* dz_node[Nnodes] is accessed later despite not being set. This can 
  cause run-time errors on some platforms (linux/alpha). Therefore, 
  set it to zero */
  soil_con->dz_node[Nnodes]=0.0;
/** End of changes -- JS **/
    
  if ( options.GRND_FLUX ) {

    for ( veg = 0 ; veg <= Nveg ; veg++) {
      for( band = 0; band < options.SNOW_BAND; band++ ) {
	
	/** Set soil properties for all soil nodes **/
	if(FIRST_VEG) {
	  FIRST_VEG = FALSE;
	  set_node_parameters(soil_con->dz_node, soil_con->max_moist_node,
			      soil_con->expt_node, soil_con->bubble_node,
			      soil_con->alpha, soil_con->beta,
			      soil_con->gamma, soil_con->depth,
			      soil_con->max_moist, soil_con->expt, 
			      soil_con->bubble, soil_con->quartz, 
			      soil_con->layer_node_fract,
#if QUICK_FS
			      soil_con->ufwc_table_node,
#endif
			      Nnodes, options.Nlayer, soil_con->FS_ACTIVE);
	  
	  sum = soil_con->dz_node[0]/2. + soil_con->dz_node[Nnodes-1]/2.;
	  for(nidx=1;nidx<Nnodes-1;nidx++) sum += soil_con->dz_node[nidx];
	}
	
	/* set soil moisture properties for all soil thermal nodes */
	distribute_node_moisture_properties(energy[veg][band].moist,
					    energy[veg][band].ice,
					    energy[veg][band].kappa_node,
					    energy[veg][band].Cs_node,
					    soil_con->dz_node,
					    energy[veg][band].T,
					    soil_con->max_moist_node,
#if QUICK_FS
					    soil_con->ufwc_table_node,
#else
					    soil_con->expt_node,
					    soil_con->bubble_node,
#endif
					    moist[veg][band], soil_con->depth,
					    soil_con->soil_density,
					    soil_con->bulk_density,
					    soil_con->quartz,
					    Nnodes, options.Nlayer,
					    soil_con->FS_ACTIVE);
	
	/* initialize layer moistures and ice contents */
	for(dry = 0; dry < Ndist; dry++) {
	  for(lidx=0;lidx<options.Nlayer;lidx++) {
	    if(options.INIT_STATE == FALSE){
	      cell[dry][veg][band].layer[lidx].moist = moist[veg][band][lidx];
	    }
	    cell[dry][veg][band].layer[lidx].ice = ice[veg][band][lidx];
	  }
	  estimate_layer_ice_content(cell[dry][veg][band].layer,
				     soil_con->dz_node,
				     energy[veg][band].T,
				     soil_con->max_moist_node,
#if QUICK_FS
				     soil_con->ufwc_table_node,
#else
				     soil_con->expt_node,
				     soil_con->bubble_node,
#endif
				     soil_con->depth,
				     soil_con->max_moist,
#if QUICK_FS
				     soil_con->ufwc_table_layer,
#else
				     soil_con->expt,
				     soil_con->bubble,
#endif
				     soil_con->bulk_density,
				     soil_con->soil_density,
				     soil_con->quartz, 
				     soil_con->layer_node_fract,
				     Nnodes, options.Nlayer, 
				     soil_con->FS_ACTIVE);
	  
	}
	
	/* Find freezing and thawing front depths */
	if(!options.QUICK_FLUX && soil_con->FS_ACTIVE) 
	  find_0_degree_fronts(&energy[veg][band], soil_con->dz_node,
			       energy[veg][band].T, Nnodes);
	
      }
    }	
  }  
}
