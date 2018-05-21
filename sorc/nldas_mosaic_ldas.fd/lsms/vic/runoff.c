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


void runoff(layer_data_struct *layer_wet,
	    layer_data_struct *layer_dry,
            energy_bal_struct *energy,
            soil_con_struct   *soil_con,
            float            *runoff_wet, 
	    float            *runoff_dry, 
	    float            *baseflow_wet,
	    float            *baseflow_dry,
	    float            *ppt, 
	    float             mu,
	    int                dt,
            int                Nnodes,
	    int                band,
	    int                iveg,
	    int                Nlayer,
	    int                full_energy_flag,
	    int                frozen_soil_flag)
/**********************************************************************
	runoff.c	Keith Cherkauer		May 18, 1996

  This subroutine calculates infiltration and runoff from the surface,
  gravity driven drainage between all soil layers, and generates 
  baseflow from the bottom layer..
  
  sublayer indecies are always [layer number][sublayer number]
  [layer number] is the current VIC model moisture layer
  [sublayer number] is the current sublayer number where: 
         0 = thawed sublayer, 1 = frozen sublayer, and 2 = unfrozen sublayer.
	 when the model is run withoputfrozen soils, the sublayer number
	 is always = 2 (unfrozen).

  UNITS:	Ksat (mm/day)
		Q12  (mm/time step)
		moist (mm)
		inflow (mm)
                runoff (mm)

  Variables:
	ppt	incoming precipitation and snow melt
	mu	fraction of area that receives precipitation
	inflow	incoming water corrected for fractional area of precip (mu)

  MODIFICATIONS:
    5/22/96	Routine modified to account for spatially varying
		precipitation, and it's effects on runoff.	KAC
    11/96		Code modified to account for extra model layers
  		needed for frozen soils modeling.		KAC
    1/9/97	Infiltration and other rate parameters modified
		for time scales of less than 1 day.		KAC
    4-1-98 Soil moisture transport is now done on an hourly time
           step, irregardless to the model time step, to prevent
           numerical stabilities in the solution	Dag and KAC
    01-24-00    simplified handling of soil moisture for the
                frozen soil algorithm.  all option selection
		now use the same soil moisture transport method   KAC
    06-Sep-03   Changed calculation of dt_baseflow to go to zero when
                soil liquid moisture <= residual moisture.  Changed
                block that handles case of total soil moisture < residual
                moisture to not allow dt_baseflow to go negative.  TJB
    17-May-04   Changed block that handles baseflow when soil moisture
                drops below residual moisture.  Now, the block is only
                entered if baseflow > 0 and soil moisture < residual,
                and the amount of water taken out of baseflow and given
                to the soil cannot exceed baseflow.  In addition, error
                messages are no longer printed, since it isn't an error
                to be in that block.                            TJB

**********************************************************************/
{  

  char               ErrStr[MAXSTRING];
  int                firstlayer, lindex;
  int                i;
  int                last_layer[MAX_LAYERS*3];
  int                last_index;
  int                last_cnt;
  int                time_step;
  int                Ndist;
  int                dist;
  int                tmplayer;
  float             ex, A, i_0, basis, frac;
  float             inflow;
  float             resid_moist[MAX_LAYERS];
  float             moist[MAX_LAYERS];
  float             ice[MAX_LAYERS];
  float             max_moist[MAX_LAYERS];
  float             max_infil;
  float             Ksat[MAX_LAYERS];
  float             Q12[MAX_LAYERS-1];
  float             Dsmax;
  float             top_moist;
  float             top_max_moist;
  float             tmp_inflow;
  float             tmp_moist;
  float             dt_inflow, dt_outflow;
  float             dt_runoff;
  float            *runoff;
  float            *baseflow;
  float             tmp_mu;
  float             dt_baseflow;
  layer_data_struct *layer;
  layer_data_struct  tmp_layer;

#if TIMESTEPSECS
  int nsteps;
#endif

  /** Set Residual Moisture **/
  if(soil_con->resid_moist[0] > SMALL) 
    for(i=0;i<Nlayer;i++) resid_moist[i] = soil_con->resid_moist[i] 
				    * soil_con->depth[i] * 1000.;
  else for(i=0;i<Nlayer;i++) resid_moist[i] = 0.;

  /** Initialize Other Parameters **/
  //  if(options.DIST_PRCP) Ndist = 2;
  //else Ndist = 1;
  Ndist = 1; 
  tmp_mu = mu;
  
  /** Allocate and Set Values for Soil Sublayers **/
  for(dist=0;dist<Ndist;dist++) {

    if(dist>0) {
      layer    = layer_dry;
      runoff   = runoff_dry;
      baseflow = baseflow_dry;
      mu       = (1. - mu);
    }
    else {
      layer    = layer_wet;
      runoff   = runoff_wet;
      baseflow = baseflow_wet;
    }
    *baseflow = 0;

    if(mu>0.) {
      
      /** ppt = amount of liquid water coming to the surface **/
      inflow = ppt[dist];
      
      /**************************************************
	Initialize Variables
      **************************************************/
      for(lindex=0;lindex<Nlayer;lindex++) {
#if TIMESTEPSECS
        /* convert ksat in mm/day to mm/timestep */
        if(dt > 3600.0)
          Ksat[lindex]         = soil_con->Ksat[lindex] / 24;
        else
          Ksat[lindex]         = soil_con->Ksat[lindex] * dt / (24. * 3600);
#else
          Ksat[lindex]         = soil_con->Ksat[lindex] / 24.;
#endif
#if LOW_RES_MOIST
	b[lindex]            = (soil_con->expt[lindex] - 3.) / 2.;
#endif

	/** Set Layer Unfrozen Moisture Content **/
	moist[lindex] = layer[lindex].moist - layer[lindex].ice;
	if(moist[lindex]<0) {
	  if(fabs(moist[lindex]) < 1e-6)
	    moist[lindex] = 0;
	  else if (layer[lindex].moist < layer[lindex].ice)
	    moist[lindex] = 0;
	  else {
	    sprintf(ErrStr,
		    "Layer %i has negative soil moisture, %f", 
		    lindex, moist[lindex]);
	    //vicerror(ErrStr);
	  }
	}

	/** Set Layer Ice Content **/
	ice[lindex]       = layer[lindex].ice;

	/** Set Layer Maximum Moisture Content **/
	max_moist[lindex] = soil_con->max_moist[lindex];
	
      }
      
      /******************************************************
        Runoff Based on Soil Moisture Level of Upper Layers
      ******************************************************/

      top_moist = 0.;
      top_max_moist=0.;
      for(lindex=0;lindex<Nlayer-1;lindex++) {
	top_moist += (moist[lindex] + ice[lindex]);
	top_max_moist += max_moist[lindex];
      }
      if(top_moist>top_max_moist) top_moist = top_max_moist;
      
      /**************************************************
        Calculate Runoff from Surface
      **************************************************/
      
      /** Runoff Calculations for Top Layer Only **/
      /** A and i_0 as in Wood et al. in JGR 97, D3, 1992 equation (1) **/
      
      max_infil = (1.0+soil_con->b_infilt) * top_max_moist;
      
      ex        = soil_con->b_infilt / (1.0 + soil_con->b_infilt);
      A         = 1.0 - pow((1.0 - top_moist / top_max_moist),ex);
      i_0       = max_infil * (1.0 - pow((1.0 - A),(1.0 
						    / soil_con->b_infilt))); 
      /* Maximum Inflow */
      
      /** equation (3a) Wood et al. **/
      
      if (inflow == 0.0) *runoff = 0.0;
      else if (max_infil == 0.0) *runoff = inflow;
      else if ((i_0 + inflow) > max_infil) 
	*runoff = inflow - top_max_moist + top_moist;
      
      /** equation (3b) Wood et al. (wrong in paper) **/
      else {
	basis = 1.0 - (i_0 + inflow) / max_infil;
	*runoff = (inflow - top_max_moist + top_moist + top_max_moist
		   * pow(basis,1.0*(1.0+soil_con->b_infilt)));
      }
      if(*runoff<0.) *runoff=0.;
      
      /**************************************************
        Compute Flow Between Soil Layers (using an hourly time step)
      **************************************************/

#if TIMESTEPSECS
      if(dt > 3600.0)
        dt_inflow = inflow / (dt / 3600.0);
      else
        dt_inflow = inflow;
#else
      dt_inflow  =  inflow / (float) dt;
#endif
      dt_outflow =  0.0;
      
#if TIMESTEPSECS
      /* calculate number of hourly runoff timesteps if the global timestep is > 1 hour */
      /* this assumes that the global timestep is a multiple of hours */
      nsteps = 1;
      if(dt > 3600.0) nsteps = dt/3600.0;
      for (time_step = 1; time_step <= nsteps; time_step++) {
#else
      for (time_step = 0; time_step < dt; time_step++) {
#endif
	inflow   = dt_inflow;
	last_cnt = 0;
	
	/*************************************
          Compute Drainage between Sublayers 
        *************************************/
	
	for( lindex = 0; lindex < Nlayer-1; lindex++ ) {
	      
	  /** Brooks & Corey relation for hydraulic conductivity **/

#if TIMESTEPSECS
          if( (tmp_moist = moist[lindex] - layer[lindex].evap / (float)nsteps)
#else
	  if((tmp_moist = moist[lindex] - layer[lindex].evap / (float)dt) 
#endif
	     < resid_moist[lindex])
	    tmp_moist = resid_moist[lindex];
	  
	  if(moist[lindex] > resid_moist[lindex]) {

	    Q12[lindex] 
	      = Ksat[lindex] * pow(((tmp_moist - resid_moist[lindex]) 
				    / (soil_con->max_moist[lindex] 
				       - resid_moist[lindex])),
				   soil_con->expt[lindex]); 
	  }
	  else Q12[lindex] = 0.;
	  last_layer[last_cnt] = lindex;
	}

	/**************************************************
          Solve for Current Soil Layer Moisture, and
          Check Versus Maximum and Minimum Moisture
          Contents.  
	**************************************************/
	
	firstlayer = TRUE;
	last_index = 0;
	for(lindex=0;lindex<Nlayer-1;lindex++) {

#if TIMESTEPSECS
          if(lindex==0) dt_runoff = *runoff / (float)nsteps;
#else
	  if(lindex==0) dt_runoff = *runoff / (float) dt;
#endif
	  else dt_runoff = 0;

	  tmp_inflow = 0.;

	  /** Update soil layer moisture content **/
#if TIMESTEPSECS
          moist[lindex] = moist[lindex] + (inflow - dt_runoff)
            - (Q12[lindex] + layer[lindex].evap/(float)nsteps);
#else
	  moist[lindex] = moist[lindex] + (inflow - dt_runoff) 
	    - (Q12[lindex] + layer[lindex].evap/(float)dt);
#endif

	  /** Verify that soil layer moisture is less than maximum **/
	  if((moist[lindex]+ice[lindex]) > max_moist[lindex]) {
	    tmp_inflow = (moist[lindex]+ice[lindex]) - max_moist[lindex];
	    moist[lindex] = max_moist[lindex] - ice[lindex];

	    if(lindex==0) {
	      Q12[lindex] += tmp_inflow;
	      tmp_inflow = 0;
	    }
	    else {
	      tmplayer = lindex;
	      while(tmp_inflow > 0) {
		tmplayer--;

		if(tmplayer<0) {
		  /** If top layer saturated, add to top layer **/
		  *runoff += tmp_inflow;
		  tmp_inflow = 0;
		}
		else {
		  /** else add excess soil moisture to next higher layer **/
		  moist[tmplayer] += tmp_inflow;
		  if((moist[tmplayer]+ice[tmplayer]) > max_moist[tmplayer]) {
		    tmp_inflow = ((moist[tmplayer] + ice[tmplayer])
				  - max_moist[tmplayer]);
		    moist[tmplayer] = max_moist[tmplayer] - ice[tmplayer];
		  }
		  else tmp_inflow=0;
		}
	      }
	    } /** end trapped excess moisture **/
	  } /** end check if excess moisture in top layer **/
	      
	  firstlayer=FALSE;
	  
	  /** verify that current layer moisture is greater than minimum **/
	  if ((moist[lindex]+ice[lindex]) < resid_moist[lindex]) {
	    /** moisture cannot fall below residual moisture content **/
	    Q12[lindex] += moist[lindex] - resid_moist[lindex];
	    moist[lindex] = resid_moist[lindex];
	  }
	      
	  inflow = (Q12[lindex]+tmp_inflow);
	  Q12[lindex] += tmp_inflow;
	      
	  last_index++;
	      
	} /* end loop through soil layers */
	  
	/**************************************************
	   Compute Baseflow
        **************************************************/
      
	/** ARNO model for the bottom soil layer (based on bottom
	 soil layer moisture from previous time step) **/
      
	lindex = Nlayer-1;
#if TIMESTEPSECS
        if(dt > 3600.0)
          Dsmax = soil_con->Dsmax / 24.;
        else
          Dsmax = soil_con->Dsmax * dt / (24. * 3600.0);
#else
        Dsmax = soil_con->Dsmax / 24.;
#endif

	frac = soil_con->Ds * Dsmax 
	  / (soil_con->Ws * soil_con->max_moist[lindex]);
	dt_baseflow = frac * ( moist[lindex] - resid_moist[lindex]);
	if (moist[lindex] > soil_con->Ws * soil_con->max_moist[lindex]) {
	  frac = (moist[lindex] - soil_con->Ws * soil_con->max_moist[lindex]) 
	    / (soil_con->max_moist[lindex] - soil_con->Ws 
	       * soil_con->max_moist[lindex]);
	  dt_baseflow += (Dsmax - soil_con->Ds * Dsmax / soil_con->Ws) 
	    * pow(frac,soil_con->c);
	}

        /** Make sure baseflow isn't negative **/
	if(dt_baseflow < 0) dt_baseflow = 0;
      
	/** Extract baseflow from the bottom soil layer **/ 
	
#if TIMESTEPSECS
	moist[lindex] += Q12[lindex-1] - (layer[lindex].evap/(float)nsteps 
					  + dt_baseflow);
#else
        moist[lindex] += Q12[lindex-1] - (layer[lindex].evap/(float)dt
                                          + dt_baseflow);		
#endif
      
	/** Check Lower Sub-Layer Moistures **/
	tmp_moist = 0;
	
        /* If baseflow > 0 and soil moisture has gone below residual,
         * take water out of baseflow and add back to soil to make up the difference
         * Note: if baseflow is small, soil moisture may still be < residual after this */
        if(dt_baseflow > 0 && (moist[lindex]+ice[lindex]) < resid_moist[lindex]) {
          if ( dt_baseflow > resid_moist[lindex] - (moist[lindex]+ice[lindex]) ) {
            dt_baseflow -= resid_moist[lindex] - (moist[lindex]+ice[lindex]);
            moist[lindex] += resid_moist[lindex] - (moist[lindex]+ice[lindex]);
          }
          else {
            moist[lindex] += dt_baseflow;
            dt_baseflow = 0.0;
          }
        }

	if((moist[lindex]+ice[lindex]) > max_moist[lindex]) {
	  /* soil moisture above maximum */
	  tmp_moist = ((moist[lindex]+ice[lindex]) - max_moist[lindex]);
	  moist[lindex] = max_moist[lindex] - ice[lindex];
	  tmplayer = lindex;
	  while(tmp_moist > 0) {
	    tmplayer--;

	    if(tmplayer<0) {
	      /** If top layer saturated, add to top layer **/
	      *runoff += tmp_moist;
	      tmp_moist = 0;
	    }
	    else {
	      /** else if sublayer exists, add excess soil moisture **/
	      moist[tmplayer] += tmp_moist ;
	      if((moist[tmplayer]+ice[tmplayer]) > max_moist[tmplayer]) {
		tmp_moist = ((moist[tmplayer] + ice[tmplayer])
				 - max_moist[tmplayer]);
		moist[tmplayer] = max_moist[tmplayer] - ice[tmplayer];
	      }
	      else tmp_moist=0;
	    }
	  }
	}

	for(lindex=0;lindex<Nlayer;lindex++) 
	  layer[lindex].moist      = moist[lindex] + ice[lindex];      
	*baseflow += dt_baseflow;

      } /* end of hourly time step loop */
      
    }
  } /** Loop over wet and dry fractions **/

  /** Recompute Thermal Parameters Based on New Moisture Distribution **/
  if((full_energy_flag==1) || (frozen_soil_flag==1)) {

    for(lindex=0;lindex<Nlayer;lindex++) {
      tmp_layer = find_average_layer(&(layer_wet[lindex]), 
				     &(layer_dry[lindex]), 
				     soil_con->depth[lindex], tmp_mu);
      moist[lindex] = tmp_layer.moist;
    }

    distribute_node_moisture_properties(energy->moist, energy->ice,
					energy->kappa_node, energy->Cs_node,
					soil_con->dz_node, energy->T,
					soil_con->max_moist_node,
					soil_con->expt_node,
					soil_con->bubble_node, 
					moist, soil_con->depth, 
					soil_con->soil_density,
					soil_con->bulk_density,
					soil_con->quartz, Nnodes, 
					Nlayer, soil_con->FS_ACTIVE,
					frozen_soil_flag);
    
  }

}



