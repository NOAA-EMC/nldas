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
#include "vicNl.h"


/***************************************************************************

  ARNO/ARNO Model of Evaporation

  Routine to compute evaporation based on the assumption that
  evaporation is at the potential for the area which is saturated,
  and at some percentage of the potential for the area which is partial
  saturated.

  Evaporation from bare soil calculated only from uppermost layer.

  Evaporation is in mm/(time step)  --> usually 1 day or 1 hour

  modifications:
  04-Jun-04 Changed logic of evap limit check to avoid creating spurious
          condensation.  Previously, when liquid moisture < residual
          moisture, (liquid moisture - residual moisture) would be
          negative.  Any non-negative evap would be greater than this,
          resulting in evap getting set to (liquid moisture - residual
          moisture), which would be negative (i.e. condensation).
          This artificially created condensation in whatever amount
          necessary to bring liquid moisture up to residual, causing
          1) large latent heat flux, 2) incorrect surface temperatures,
          3) occasional inability for calc_surf_energy_bal to converge
          in root_brent, and 4) spuriously high runoff and baseflow.
            Now there is an added condition that liquid moisture > residual
          moisture for evap to be capped at (liquid moisture - residual
          moisture).                                                  TJB

****************************************************************************/

float arno_evap(layer_data_struct *layer_wet,
		 layer_data_struct *layer_dry,
		 float             rad,
		 float             air_temp,
		 float             vpd,
		 float             net_short,
		 float             D1,
		 float             max_moist,
		 float             elevation,
		 float             b_infilt,
		 float             Tsurf,
		 float             displacement,
		 float             roughness,
		 float             ref_height,
		 float             ra,
		 float             dt,
		 float             mu,
		 float             moist_resid)
{

  int    num_term;
  int    i;
  int    Ndist;
  int    dist;
  float tmp,beta_asp,dummy;
  float ratio,as;
  float Epot;		/* potential bare soil evaporation */
  float moist;
  float evap;
  float max_infil;
  float Evap;
  float tmpsum;
  layer_data_struct *layer;

  //  if(options.DIST_PRCP) Ndist = 2;
  //  else Ndist = 1;
  Ndist = 1; 

  Evap = 0;

  for(dist=0;dist<Ndist;dist++) {

    if(dist>0) {
      mu = (1. - mu);
      layer = layer_dry;
    }
    else {
      layer = layer_wet;
    }
    
    /* moist = liquid soil moisture */
    moist = layer[0].moist - layer[0].ice;
    if(moist>max_moist) moist=max_moist;

    /* Calculate the potential bare soil evaporation (mm/time step) */
#if TIMESTEPSECS
    /* convert from mm/day to mm */
    Epot = vic_penman(rad, vpd * 1000., ra, 0.0, 0.0, 1.0, 1.0,
                  air_temp, net_short, elevation, 0.0) * dt / (24.0 * SECPHOUR);
#else
    Epot = vic_penman(rad, vpd * 1000., ra, 0.0, 0.0, 1.0, 1.0, 
		  air_temp, net_short, elevation, 0.0) * dt / 24.0;
#endif

    /**********************************************************************/
    /*  Compute temporary infiltration rate based on given soil_moist.    */
    /**********************************************************************/
    max_infil = (1.0+b_infilt)*max_moist;
    if(b_infilt == -1.0)
      tmp = max_infil;
    else {
      ratio = 1.0 - (moist) / (max_moist);
      /*****if(ratio < SMALL && ratio > -SMALL) ratio = 0.;*****/
      if(ratio > 1.0) {
	printf("\n  ERROR: SOIL RATIO GREATER THAN 1.0\n");
	printf("moisture %f   max_moisture %f -> ratio = %f\n",
	       moist,max_moist,ratio);
	exit(0);
      }
      else {
	if(ratio < 0.0) {
	  printf("\n  ERROR: SOIL RATIO LESS THAN 0.0\n");
	  printf("moisture %f   max_moisture %f -> ratio = %e\n",
		 moist,max_moist,ratio);
	  exit(0);
	}
	else
	  ratio = pow(ratio,(1.0 / (b_infilt + 1.0)));
      }
      tmp = max_infil*(1.0 - ratio);
    }

    /************************************************************************/
    /* Evaporate at potential rate, i.e., Eq.(10) in Liang's derivation.    */
    /************************************************************************/

    if(tmp >= max_infil)
      evap = Epot;     
    else {                 
    
      /********************************************************************/
      /*  Compute As. 'As' is % area saturated, '1-As' is % area          */
      /*  that is unsaturated.                                            */
      /********************************************************************/
    
      ratio = tmp/max_infil; 
      ratio = 1.0 - ratio;
      
      if(ratio > 1.0) {
	printf("\n ARNO ERROR: EVAP RATIO GREATER THAN 1.0");
	exit(0);
      }
      else {
	if(ratio < 0.0) {
	  printf("\n ARNO ERROR: EVAP RATIO LESS THAN 0.0");
	  exit(0);
	}
	else {
	  if(ratio != 0.0)
	    ratio = pow(ratio,b_infilt);
	}
      }
    
      as = 1 - ratio;

      /********************************************************************/
      /*  Compute the beta function in the ARNO evaporation model using   */
      /*  the first 30 terms in the power expansion expression.           */
      /********************************************************************/ 
	
      ratio = pow(ratio,(1.0/b_infilt));
    
      dummy = 1.0;
      for(num_term=1;num_term<=30;num_term++) {
	tmpsum = ratio;
	for(i=1;i<num_term;i++) tmpsum *= ratio;
	dummy += b_infilt * tmpsum /
	  (b_infilt + num_term);
      }
    
      beta_asp = as+(1.0-as)*(1.0-ratio)*dummy;
      evap = Epot*beta_asp;
    }
	
    /***********************************************************************/
    /*  Evaporation cannot exceed available soil moisture.                 */
    /*  Evaporation second soil layer = 0.0                                */
    /***********************************************************************/


    /* only consider positive evaporation; we won't put limits on condensation */
    if (evap > 0.0) {
      if (moist > moist_resid * D1 * 1000.) {
        /* there is liquid moisture available; cap evap at available liquid moisture */
        if (evap > moist -  moist_resid * D1 * 1000.) {
          evap = moist -  moist_resid * D1 * 1000.;
        }
      }
      else {
        /* no moisture available; cap evap at 0 */
        evap = 0.0;
      }
    }

    layer[0].evap = evap;
#if TIMESTEPSECS
        /* convert from mm to m/timestep */
        Evap += evap / 1000. / dt * mu;
#else
    Evap += evap / 1000. / dt / 3600. * mu;
#endif

  }
  
  return(Evap);

}
