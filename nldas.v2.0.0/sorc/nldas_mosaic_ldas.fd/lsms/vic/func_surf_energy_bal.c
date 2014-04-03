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
#include <stdarg.h>
#include "vicNl.h"

float func_surf_energy_bal(float Ts, va_list ap)
/**********************************************************************
	func_surf_energy_bal	Keith Cherkauer		January 3, 1996

  This subroutine computes the surface energy balance for bare soil
  and vegetation uncovered by snow.  It computes outgoing longwave,
  sensible heat flux, ground heat flux, and storage of heat in the thin
  upper layer, based on the given surface temperature.

  The Energy Balance Equation used comes from Xu Liang's Paper 
  "Insights of the Ground Heat Flux in Land Surface Parameterization
  Schemes."

  Modifications:
  04-14-98 modified to compute evapotranspiration within this routine
           in the hopes of reducing the number of iteration 
	  needed to find a solution surface temperature.       KAC
  07-13-98 modified to include elevation bands for vegetation 
           and snow                                             KAC
  01-20-00 modified to work with the updated radiation estimation
           routines, as well as the simplified frozen soil moisture
           storage                                              KAC

**********************************************************************/
{
  /** Thermal Properties **/
  float             T2;       	/** average soil temperature (C) **/
  float             Ts_old;	/** last temperature (C) **/
  float             T1_old;	/** last layer 1 soil temperature (C) **/
  float             Tair;     	/** Air Temperature **/
  float             ra;       	/** aerodynamic reisistance (s/m) **/
  float             atmos_density;	/** atmospheric density (kg/m^3) **/
  float             shortwave;
  float             longwave;
  float             albedo;
  float             emissivity;
  float             kappa1;	/** thermal conductivity of 1st layer */
  float             kappa2;	/** thermal conductivity of 2nd layer */
  float             Cs1;      	/** volumetric heat capacity of 1st layer **/
  float             Cs2;      	/** volumetric heat capacity of 2nd layer **/
  float             D1;       	/** thickness of 1st layer (m) **/
  float             D2;       	/** thickness of 2nd layer (m) **/
  float             dp;       	/** depth to constant temperature (m) */
  float             delta_t;	/** Time Step in Seconds **/
  float             Le;       	/** Latent heat of vapoization **/
  float             Ls;       	/** Latent heat of sublimation **/
  float             Vapor;    	/** Total vapor flux from snow in m/s **/
  float             moist;    	/** layer moisture content in m/m **/
  float             ice0;     	/** layer ice content in m/m **/
  float             max_moist;	/** layer maximum moisture content in m/m **/
  float             bubble;	/** bubbling pressure in cm **/
  float             expt;
  float             snow_depth;
  float             snow_density;
  float             Tsnow_surf;
  float             snow_cover_fraction;
  float             surf_atten;
  float             wind;
  float             displacement;
  float             roughness;
  float             ref_height;
  float             elevation;
  float             b_infilt;
  float             max_infil;
  float             dt;
  float             vpd;
  float             snow_energy;
  float             mu;
  float            *rainfall;
  float            *Wdew;
  float            *grnd_flux;
  float            *T1;
  float            *latent_heat;
  float            *sensible_heat;
  float            *deltaH;
  float            *snow_flux;
  float            *store_error;
  float             TMean;
  float             rad;
  float            *depth;
  float            *Wcr;
  float            *Wpwp;
  float            *resid_moist;
  float            *T_node;
  float            *Tnew_node;
  float            *dz_node;
  float            *kappa_node;
  float            *Cs_node;
  float            *moist_node;
  float            *bubble_node;
  float            *expt_node;
  float            *max_moist_node;
  float            *ice_node;
  float            *alpha;
  float            *beta;
  float            *gamma;
#if QUICK_FS
  float           **ufwc_table_layer;
  float          ***ufwc_table_node;
#endif
  float             *root;
  layer_data_struct *layer_wet;
  layer_data_struct *layer_dry;
  veg_var_struct    *veg_var_wet;
  veg_var_struct    *veg_var_dry;
  int                VEG;
  int                veg_class;
  int                month;
  int                Nnodes;
  int                Nlayer;
  char              *FIRST_SOLN;
  int                SNOWING;
  int                FS_ACTIVE;
  int                frozen_soil_flag;
  int                grnd_flux_flag;
  int                quick_flux_flag;

  float             error;
  float             ice;
  float             Evap;		/** Total evap in m/s **/
  float             kappa_snow;

  T2                  = (float) va_arg(ap, double);
  Ts_old              = (float) va_arg(ap, double);
  T1_old              = (float) va_arg(ap, double);
  Tair                = (float) va_arg(ap, double);
  ra                  = (float) va_arg(ap, double);
  atmos_density       = (float) va_arg(ap, double);
  shortwave           = (float) va_arg(ap, double);
  longwave            = (float) va_arg(ap, double);
  albedo              = (float) va_arg(ap, double);
  emissivity          = (float) va_arg(ap, double);
  kappa1              = (float) va_arg(ap, double);
  kappa2              = (float) va_arg(ap, double);
  Cs1                 = (float) va_arg(ap, double);
  Cs2                 = (float) va_arg(ap, double);
  D1                  = (float) va_arg(ap, double);
  D2                  = (float) va_arg(ap, double);
  dp                  = (float) va_arg(ap, double);
  delta_t             = (float) va_arg(ap, double);
  Le                  = (float) va_arg(ap, double);
  Ls                  = (float) va_arg(ap, double);
  Vapor               = (float) va_arg(ap, double);
  moist               = (float) va_arg(ap, double);
  ice0                = (float) va_arg(ap, double);
  max_moist           = (float) va_arg(ap, double);
  bubble              = (float) va_arg(ap, double);
  expt                = (float) va_arg(ap, double);
  snow_depth          = (float) va_arg(ap, double);
  snow_density        = (float) va_arg(ap, double);
  Tsnow_surf          = (float) va_arg(ap, double);
  snow_cover_fraction = (float) va_arg(ap, double);
  surf_atten          = (float) va_arg(ap, double);
  wind                = (float) va_arg(ap, double);
  displacement        = (float) va_arg(ap, double);
  roughness           = (float) va_arg(ap, double);
  ref_height          = (float) va_arg(ap, double);
  elevation           = (float) va_arg(ap, double);
  b_infilt            = (float) va_arg(ap, double);
  max_infil           = (float) va_arg(ap, double);
  dt                  = (float) va_arg(ap, double);
  vpd                 = (float) va_arg(ap, double);
  snow_energy         = (float) va_arg(ap, double);
  mu                  = (float) va_arg(ap, double);
  rainfall            = (float *) va_arg(ap, double *);
  Wdew                = (float *) va_arg(ap, double *);
  grnd_flux           = (float *) va_arg(ap, double *);
  T1                  = (float *) va_arg(ap, double *);
  latent_heat         = (float *) va_arg(ap, double *);
  sensible_heat       = (float *) va_arg(ap, double *);
  deltaH              = (float *) va_arg(ap, double *);
  snow_flux           = (float *) va_arg(ap, double *);
  store_error         = (float *) va_arg(ap, double *);
  depth               = (float *) va_arg(ap, double *);
  Wcr                 = (float *) va_arg(ap, double *);
  Wpwp                = (float *) va_arg(ap, double *);
  resid_moist         = (float *) va_arg(ap, double *);
  T_node              = (float *) va_arg(ap, double *);
  Tnew_node           = (float *) va_arg(ap, double *);
  dz_node             = (float *) va_arg(ap, double *);
  kappa_node          = (float *) va_arg(ap, double *);
  Cs_node             = (float *) va_arg(ap, double *);
  moist_node          = (float *) va_arg(ap, double *);
  bubble_node         = (float *) va_arg(ap, double *);
  expt_node           = (float *) va_arg(ap, double *);
  max_moist_node      = (float *) va_arg(ap, double *);
  ice_node            = (float *) va_arg(ap, double *);
  alpha               = (float *) va_arg(ap, double *);
  beta                = (float *) va_arg(ap, double *);
  gamma               = (float *) va_arg(ap, double *);
  root                = (float  *) va_arg(ap, double  *);
  layer_wet           = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  layer_dry           = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var_wet         = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  veg_var_dry         = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  VEG                 = (int) va_arg(ap, int);
  veg_class           = (int) va_arg(ap, int);
  month               = (int) va_arg(ap, int);
  Nnodes              = (int) va_arg(ap, int);
  Nlayer              = (int) va_arg(ap, int);
  FIRST_SOLN          = (char *)va_arg(ap, char *);
  SNOWING             = (int) va_arg(ap, int);
  FS_ACTIVE           = (int) va_arg(ap, int);
  frozen_soil_flag    = (int) va_arg(ap, int);
  grnd_flux_flag    = (int) va_arg(ap, int);
  quick_flux_flag    = (int) va_arg(ap, int);
  
  TMean = Ts; 
    if(grnd_flux_flag==1) {
  
    /**********************************************
      Compute Surface Temperature at Half Time Step
    **********************************************/
    if(snow_cover_fraction > 0) {

      /****************************************
        Compute energy flux through snow pack
      ****************************************/

      kappa_snow = 2.9302e-6 * (snow_density) * (snow_density);
 
      *snow_flux = kappa_snow * (TMean - Tsnow_surf) / snow_depth;

    }

    if(quick_flux_flag==1) {
      /**************************************************************
        Use Liang et al. 1999 Equations to Calculate Ground Heat 
	Flux
      **************************************************************/
      *T1 = estimate_T1(TMean, T1_old, T2, D1, D2, kappa1, kappa2, Cs1, 
			Cs2, dp, delta_t);
    
    }
    else {

    /*************************************************************
      Use Finite Difference Method to Explicitly Solve Ground Heat
      Flux at Soil Thermal Nodes (Cherkauer and Lettenmaier, 1999)
    *************************************************************/
      T_node[0] = TMean;
      solve_T_profile(Tnew_node, T_node, dz_node, kappa_node, Cs_node, 
		      moist_node, delta_t, max_moist_node, bubble_node, 
		      expt_node, ice_node, alpha, beta, gamma, Nnodes, 
		      FIRST_SOLN, FALSE, FS_ACTIVE,frozen_soil_flag);
      *T1 = Tnew_node[1];
    }

    /*****************************************************
      Compute the Ground Heat Flux from the Top Soil Layer
    *****************************************************/
/*     *grnd_flux = surf_atten * (kappa1/D1*((*T1) - (TMean))); */
/*     *grnd_flux = (kappa1/D1*((*T1) - (TMean))); */
    *grnd_flux = (snow_cover_fraction + (1. - snow_cover_fraction) 
		  * surf_atten) * (kappa1 / D1 * ((*T1) - TMean));

    /******************************************************
      Compute the Current Ice Content of the Top Soil Layer
    ******************************************************/
    if((FS_ACTIVE && (frozen_soil_flag==1)) && (TMean+ *T1)/2.<0.) {
      ice = moist - maximum_unfrozen_water((TMean+ *T1)/2.,
					   max_moist,bubble,expt);
      if(ice<0.) ice=0.;
    }
    else ice=0.;
//<Justin's frozen soil fix>
    //ice = 0.0;
//</Justin's frozen soil fix>
 
    *deltaH = Cs1 * ((Ts_old + T1_old)/2. - (TMean + *T1)/2.) * D1 / delta_t;
    /* Only adjust for ice if both FS_ACTIVE and FROZEN_SOIL are true */
    if(FS_ACTIVE && frozen_soil_flag)
      *deltaH -= ice_density*Lf*(ice0-ice)*D1/delta_t;

    /** Compute net surface radiation for evaporation estimates **/
    rad = (1.0 - albedo) * shortwave + longwave 
      - STEFAN_B * (TMean+KELVIN) * (TMean+KELVIN) * (TMean+KELVIN) 
      * (TMean+KELVIN) + *grnd_flux + *deltaH;

  } /* End computation for ground heat flux */
  else { /* ground heat flux not estimated */

#if TIMESTEPSECS
    if(dt < 24. * 3600.0)
#else
    if(dt < 24)
#endif

      /** Compute net surface radiation for evaporation estimates **/
      rad = (1.0 - albedo) * shortwave + longwave 
	- STEFAN_B * (TMean+KELVIN) * (TMean+KELVIN) * (TMean+KELVIN) 
	* (TMean+KELVIN);

    else

      /** Daily water balance model provides average shortwave and 
	  net longwave radiation **/
      rad = (1.0 - albedo) * shortwave + longwave;

  }

  /*************************************************
    Compute Evapotranspiration if not snow covered
  *************************************************/
    if(VEG && !SNOWING) {
    Evap = canopy_evap(layer_wet, layer_dry, veg_var_wet, veg_var_dry, TRUE, 
		       veg_class, month, mu, Wdew, dt, rad, vpd, 
		       (1.0 - albedo) * shortwave, Tair, ra, displacement, 
		       roughness, ref_height, elevation, rainfall, depth, 
		       Wcr, Wpwp, root, Nlayer);
  }
  else if(!SNOWING){
    Evap = arno_evap(layer_wet, layer_dry, rad, Tair, vpd, 
		     (1.0 - albedo) * shortwave, D1, 
		     max_moist * depth[0] * 1000., elevation, b_infilt,  
		     Tair, displacement, roughness, ref_height, ra, dt, mu,
		     resid_moist[0]);
  }
  else Evap = 0.;
  
  /**********************************************************************
    Compute the Latent Heat Flux from the Surface and Covering Vegetation
  **********************************************************************/
  *latent_heat  = -RHO_W*Le*Evap;
  *latent_heat += -atmos_density*Ls*Vapor;
  //printf("latent heat ..%i %f\n", SNOWING, Evap);

  if(grnd_flux_flag==1) {
  
    /************************************************
      Compute the Sensible Heat Flux from the Surface
    ************************************************/
    *sensible_heat = atmos_density*Cp*(Tair - (TMean))/ra;
    //    printf("sensible heat ..%f %f %f %f %f \n",atmos_density,Cp,Tair,TMean,ra);
    /*************************************
      Compute Surface Energy Balance Error
    *************************************/
    error = (1. - snow_cover_fraction) 
      * ((1.-albedo)*shortwave 
	 + emissivity*(longwave-STEFAN_B*(TMean + KELVIN)*(TMean + KELVIN)
		       *(TMean + KELVIN)*(TMean + KELVIN))
	 + *sensible_heat + *latent_heat + snow_energy) 
      - snow_cover_fraction * *snow_flux + *grnd_flux + *deltaH;
    
    *store_error = error;
  }
  else error = MISSING;
  
  return error;

}




