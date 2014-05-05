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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

float calc_surf_energy_bal(int                iveg,
			    int                nlayer,
			    int                Nveg,
			    int                dt,
			    int                Nnodes,
			    int                veg_class,
			    int                band,
			    int                hour,
			    float             ice0,
			    float             moist,
			    float             dp,
			    float             surf_atten,
			    float             T0,
			    float             shortwave,
			    float             longwave,
			    float             air_temp,
			    float             Le,
			    float             Ls,
			    float             mu,
			    float             displacement,
			    float             roughness,
			    float             ref_height,
			    float             snow_energy,
			    float             vapor_flux,
			    float             bare_albedo,
			    float            *aero_resist,
			    float            *wind,
			    float            *rainfall,
			    float            *ppt,
			    float             *root,
			    atmos_data_struct *atmos,
			    veg_var_struct    *veg_var_wet,
			    veg_var_struct    *veg_var_dry,
			    energy_bal_struct *energy,
			    snow_data_struct  *snow,
			    layer_data_struct *layer_wet,
			    layer_data_struct *layer_dry,
			    soil_con_struct   *soil_con,
			    dmy_struct        *dmy,
			    int               full_energy_flag,
			    int               quick_flux_flag,
			    int               grnd_flux_flag,
			    int               frozen_soil_flag)
/**************************************************************
  calc_surf_energy_bal.c  Greg O'Donnell and Keith Cherkauer  Sept 9 1997
  
  This function calculates the surface temperature, in the
  case of no snow cover.  Evaporation is computed using the
  previous ground heat flux, and then used to comput latent heat
  in the energy balance routine.  Surface temperature is found
  using the Root Brent method (Numerical Recipies).
  

  modifications:
    02-29-00  Included variables needed to compute energy flux
              through the snow pack.  The ground surface energy
              balance will now be a mixture of snow covered
	      and bare ground, controlled by the snow cover 
	      fraction set in solve_snow.c                 KAC

***************************************************************/
{
  extern veg_lib_struct *veg_lib;

  char     FIRST_SOLN[1];
  float   T2;
  float   Ts_old;
  float   T1_old;
  float   Tair;
  float   ra;
  float   atmos_density;
  float   albedo;
  float   emissivity;
  float   kappa1;
  float   kappa2;
  float   Cs1;
  float   Cs2;
  float   D1;
  float   D2;
  float   delta_t;
  float   Vapor;
  float   max_moist;
  float   bubble;
  float   expt;
  float   T1;
  float   snow_depth;
  float   snow_density;
  float   Tsnow_surf;
  float   snow_cover_fraction;
  /* float   snow_flux; */

  int      VEG;
  float   Tsurf;
  float   U;
  float   error;
  float   Wdew[2];
  float  *T_node;
  float  Tnew_node[MAX_NODES];
  float  *dz_node;
  float  *kappa_node;
  float  *Cs_node;
  float  *moist_node;
  float  *bubble_node;
  float  *expt_node;
  float  *max_moist_node;
  float  *ice_node;
  float  *alpha;
  float  *beta;
  float  *gamma;
  layer_data_struct layer[MAX_LAYERS];

  float   T_lower, T_upper;

  /**************************************************
    Correct Aerodynamic Resistance for Stability
  **************************************************/
  ra = aero_resist[0];
  U = wind[0];
  if (U > 0.0)
    ra /= StabilityCorrection(ref_height, displacement, T0,
          air_temp, U, roughness);
  else
    ra = HUGE_RESIST;
  
  /**************************************************
    Compute Evaporation and Transpiration 
  **************************************************/
  if(iveg!=Nveg) {
    if(veg_lib[veg_class].LAI[dmy->month-1] > 0.0) VEG = TRUE;
    else VEG = FALSE;
  }
  else VEG = FALSE;
#if TIMESTEPSECS
  Vapor = vapor_flux / (float)dt;
#else
  Vapor = vapor_flux / (float)dt / 3600.;
#endif

  /**************************************************
    Set All Variables For Use
  **************************************************/

  T2                  = energy->T[Nnodes-1];
  Ts_old              = energy->T[0];
  T1_old              = energy->T[1];
  Tair                = air_temp;
  atmos_density       = atmos->density;
  albedo              = bare_albedo;
  emissivity          = 1.;
  kappa1              = energy->kappa[0];
  kappa2              = energy->kappa[1];
  Cs1                 = energy->Cs[0];
  Cs2                 = energy->Cs[1];
  D1                  = soil_con->depth[0];
/*   D2                  = soil_con->depth[1]; */
  D2                  = soil_con->depth[0];
  delta_t             = (float)dt;
  max_moist           = soil_con->max_moist[0]/(soil_con->depth[0]*1000.);
  bubble              = soil_con->bubble[0];
  expt                = soil_con->expt[0];
  snow_depth          = snow->depth;
  snow_density        = snow->density;
  Tsnow_surf          = snow->surf_temp;
  snow_cover_fraction = snow->coverage;
  Wdew[WET]           = veg_var_wet->Wdew;
  //if(options.DIST_PRCP) Wdew[DRY] = veg_var_dry->Wdew;
  FIRST_SOLN[0] = TRUE;

  /*************************************************************
    Prepare soil node variables for finite difference solution
  *************************************************************/

 if(!quick_flux_flag) {

    bubble_node    = soil_con->bubble_node; 
    expt_node      = soil_con->expt_node; 
    max_moist_node = soil_con->max_moist_node;  
    alpha          = soil_con->alpha; 
    beta           = soil_con->beta; 
    gamma          = soil_con->gamma; 
    moist_node     = energy->moist;
    kappa_node     = energy->kappa_node;
    Cs_node        = energy->Cs_node;
    T_node         = energy->T;
    dz_node        = soil_con->dz_node;
    ice_node       = energy->ice;

 }
 else {

   bubble_node    = NULL; 
   expt_node      = NULL; 
   max_moist_node = NULL;  
   alpha          = NULL; 
   beta           = NULL; 
   gamma          = NULL; 
   moist_node     = NULL;
   kappa_node     = NULL;
   Cs_node        = NULL;
   T_node         = NULL;
   dz_node        = NULL;
   ice_node       = NULL;
   
 }

  /**************************************************
    Find Surface Temperature Using Root Brent Method
  **************************************************/
    if(full_energy_flag) {
      if(snow->swq > 0) {
	T_lower = 0.;
	T_upper = energy->T[0]-SURF_DT;
      }
      else {
	T_lower = 0.5*(energy->T[0]+air_temp)-SURF_DT;
	T_upper = 0.5*(energy->T[0]+air_temp)+SURF_DT;
      }
      
      Tsurf = root_brent(T_upper, T_lower, 
			 func_surf_energy_bal, T2, Ts_old, T1_old, Tair, 
			 ra, atmos_density, shortwave, longwave, albedo, 
			 emissivity, kappa1, kappa2, Cs1, Cs2, D1, D2, dp, 
			 delta_t, Le, Ls, Vapor, moist, ice0, max_moist, 
			 bubble, expt, snow_depth, snow_density, Tsnow_surf, 
			 snow_cover_fraction, surf_atten, U, displacement, 
			 roughness, ref_height, (float)soil_con->elevation, 
			 soil_con->b_infilt, soil_con->max_infil, 
			 (float)dt, atmos->vpd, snow_energy, mu, 
			 rainfall, Wdew, &energy->grnd_flux, &T1, 
			 &energy->latent, &energy->sensible, &energy->deltaH, 
			 &energy->snow_flux, &energy->error, soil_con->depth, 
			 soil_con->Wcr, soil_con->Wpwp, 
			 soil_con->resid_moist, T_node, Tnew_node, dz_node, 
			 kappa_node, Cs_node, moist_node, bubble_node, 
			 expt_node, max_moist_node, ice_node, alpha, beta, 
			 gamma, root, layer_wet, layer_dry, veg_var_wet, 
			 veg_var_dry, VEG, veg_class, 
			 dmy->month, Nnodes, nlayer,FIRST_SOLN, snow->snow, 
			 soil_con->FS_ACTIVE,frozen_soil_flag,
			 grnd_flux_flag,quick_flux_flag);
      //      printf("Tsurf from root brent..%f\n",Tsurf);
      if ( Tsurf == MISSING )
      {
         lis_log_msgC("ERR: calc_surf_energy_bal.c -- Tsurf equals MISSING value");
         exit(0);
      }

      if ( Tsurf <= MISSING+1 ) {  
	error = error_calc_surf_energy_bal(Tsurf, T2, Ts_old, T1_old, Tair, 
					   ra, atmos_density, shortwave, 
					   longwave, albedo, emissivity, kappa1, 
					   kappa2, Cs1, Cs2, D1, D2, dp, 
					   delta_t, Le, Ls, Vapor, moist, ice0, 
					   max_moist, bubble, expt, snow_depth, 
					   snow_density, Tsnow_surf, 
					   snow_cover_fraction, surf_atten, 
					   U, displacement, roughness, 
					   ref_height, 
					   (float)soil_con->elevation, 
					   soil_con->b_infilt, 
					   soil_con->max_infil, (float)dt, 
					   atmos->vpd, snow_energy, mu, 
					   rainfall, Wdew, &energy->grnd_flux, 
					   &T1, &energy->latent, 
					   &energy->sensible, &energy->deltaH, 
					   &energy->snow_flux,
					   &energy->error, soil_con->depth, 
					   soil_con->Wcr, soil_con->Wpwp, 
					   soil_con->resid_moist, T_node, 
					   Tnew_node, dz_node, kappa_node, 
					   Cs_node, moist_node, bubble_node, 
					   expt_node, max_moist_node, ice_node, 
					   alpha, beta, gamma, root, layer_wet, 
					   layer_dry, veg_var_wet, veg_var_dry, 
					   VEG, veg_class, 
					   dmy->month, iveg, Nnodes, FIRST_SOLN, 
					   snow->snow, soil_con->FS_ACTIVE);
      }
    }
    else {

    /** Frozen soil model run with no surface energy balance **/
    Tsurf = Tair;

    }

  /**************************************************
    Recalculate Energy Balance Terms for Final Surface Temperature
  **************************************************/

  error = solve_surf_energy_bal(Tsurf, T2, Ts_old, T1_old, Tair, ra, 
				atmos_density, shortwave, longwave, albedo, 
				emissivity, kappa1, kappa2, Cs1, Cs2, D1, D2, 
				dp, delta_t, Le, Ls, Vapor, moist, ice0, 
				max_moist, bubble, expt, snow_depth, 
				snow_density, Tsnow_surf, 
				snow_cover_fraction, surf_atten, U, 
				displacement, roughness, ref_height, 
				(float)soil_con->elevation, 
				soil_con->b_infilt, soil_con->max_infil, 
				(float)dt, atmos->vpd, snow_energy, mu, 
				rainfall, Wdew, &energy->grnd_flux, 
				&T1, &energy->latent, &energy->sensible, 
				&energy->deltaH, &energy->snow_flux, 
				&energy->error, soil_con->depth, 
				soil_con->Wcr, soil_con->Wpwp, 
				soil_con->resid_moist, T_node, 
				Tnew_node, dz_node, kappa_node, Cs_node, 
				moist_node, bubble_node, expt_node, 
				max_moist_node, ice_node, alpha, beta, gamma, 
				root, layer_wet, layer_dry, veg_var_wet, 
				veg_var_dry, VEG, 
				veg_class, dmy->month, Nnodes, nlayer, FIRST_SOLN, 
				snow->snow, soil_con->FS_ACTIVE,frozen_soil_flag,grnd_flux_flag,quick_flux_flag);
  
  energy->error = error;

  /***************************************************
    Recalculate Soil Moisture and Thermal Properties
  ***************************************************/
  if(grnd_flux_flag) {
    if(quick_flux_flag) {
      energy->T[0] = Tsurf;
      energy->T[1] = T1; 
    }
    else {
 
      finish_frozen_soil_calcs(energy, layer_wet, layer_dry, layer, soil_con, 
			       Nnodes, iveg, mu, Tnew_node, kappa_node, 
			       Cs_node, moist_node, nlayer, frozen_soil_flag);
      
    }
    
  }
  else {
  
    energy->T[0] = Tsurf;

  }

  /** Store precipitation that reaches the surface */
  if(!snow->snow) {
    if(iveg!=Nveg) {
      if(veg_lib[veg_class].LAI[dmy->month-1] <= 0.0) { 
	veg_var_wet->throughfall = rainfall[WET];
	//	if(options.DIST_PRCP) veg_var_dry->throughfall = rainfall[DRY];
	ppt[WET] = veg_var_wet->throughfall;
	//	if(options.DIST_PRCP) ppt[DRY] = veg_var_dry->throughfall;
      }
      else {
	ppt[WET] = veg_var_wet->throughfall;
	//	if(options.DIST_PRCP) ppt[DRY] = veg_var_dry->throughfall;
      }
    }
    else {
      ppt[WET] = rainfall[WET];
      //      if(options.DIST_PRCP) ppt[DRY] = rainfall[DRY];
    }
  }
   
  /** Store net longwave radiation **/
  if(hour < 24){
    energy->longwave = longwave 
      - STEFAN_B * (Tsurf+KELVIN) * (Tsurf+KELVIN) 
		    * (Tsurf+KELVIN) * (Tsurf+KELVIN);
  }
  else energy->longwave = longwave;

  /** Store surface albedo **/
  energy->albedo = bare_albedo;

  /** Store net shortwave radiation **/
  energy->shortwave = (1. - energy->albedo) * shortwave;

  /** Return soil surface temperature **/
  return (Tsurf);
    
}

float solve_surf_energy_bal(float Tsurf, ...) {

  va_list ap;

  float error;

  va_start(ap, Tsurf);
  error = func_surf_energy_bal(Tsurf, ap);
  va_end(ap);

  return error;

}

float error_calc_surf_energy_bal(float Tsurf, ...) {

  va_list ap;

  float error;

  va_start(ap, Tsurf);
  /* error = error_print_surf_energy_bal(Tsurf, ap); */
  va_end(ap);

  error = 0.0;
  return error;

}
#if 0
float error_print_surf_energy_bal(float Ts, va_list ap) {

  /* extern option_struct options; */

  float             T2; 
  float             Ts_old; 
  float             T1_old; 
  float             Tair; 
  float             ra; 
  float             atmos_density; 
  float             shortwave;
  float             longwave;
  float             albedo;
  float             emissivity;
  float             kappa1; 
  float             kappa2; 
  float             Cs1; 
  float             Cs2; 
  float             D1; 
  float             D2; 
  float             dp; 
  float             delta_t;
  float             Le; 
  float             Ls; 
  float             Vapor; 
  float             moist; 
  float             ice0; 
  float             max_moist; 
  float             bubble; 
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
  int                iveg;
  int                Nnodes;
  char              *FIRST_SOLN;
  int                SNOWING;
  int                FS_ACTIVE;


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
#if QUICK_FS
  ufwc_table_layer    = (float **) va_arg(ap, float **);
  ufwc_table_node     = (float ***) va_arg(ap, float ***);
#endif
  root                = (float  *) va_arg(ap, double  *);
  layer_wet           = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  layer_dry           = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var_wet         = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  veg_var_dry         = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  VEG                 = (int) va_arg(ap, int);
  veg_class           = (int) va_arg(ap, int);
  month               = (int) va_arg(ap, int);
  iveg                = (int) va_arg(ap, int);
  Nnodes              = (int) va_arg(ap, int);
  FIRST_SOLN          = (char *) va_arg(ap, char *);
  SNOWING             = (int) va_arg(ap, int);
  FS_ACTIVE           = (int) va_arg(ap, int);

//  rec                 = (int) va_arg(ap, int);

  sprintf(lismsg,"DBG: error_print -- T2 = %f [ %d ]",T2,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Ts_old = %f [ %d ]",Ts_old,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- T1_old = %f [ %d ]",T1_old,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Tair = %f [ %d ]",Tair,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- ra = %f [ %d ]",ra,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- atmos_density = %f [ %d ]",atmos_density,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- shortwave = %f [ %d ]",shortwave,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- longwave = %f [ %d ]",longwave,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- albedo = %f [ %d ]",albedo,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- emissivity = %f [ %d ]",emissivity,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- kappa1 = %f [ %d ]",kappa1,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- kappa2 = %f [ %d ]",kappa2,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Cs1 = %f [ %d ]",Cs1,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Cs2 = %f [ %d ]",Cs2,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- D1 = %f [ %d ]",D1,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- D2 = %f [ %d ]",D2,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- dp = %f [ %d ]",dp,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- delta_t = %f [ %d ]",delta_t,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Le = %f [ %d ]",Le,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Ls = %f [ %d ]",Ls,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Vapor = %f [ %d ]",Vapor,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- moist = %f [ %d ]",moist,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- ice0 = %f [ %d ]",ice0,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- max_moist = %f [ %d ]",max_moist,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- bubble = %f [ %d ]",bubble,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- expt = %f [ %d ]",expt,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- surf_atten = %f [ %d ]",surf_atten,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- wind = %f [ %d ]",wind,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- displacement = %f [ %d ]",displacement,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- roughness = %f [ %d ]",roughness,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- ref_height = %f [ %d ]",ref_height,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- elevation = %f [ %d ]",elevation,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- b_infilt = %f [ %d ]",b_infilt,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- max_infil = %f [ %d ]",max_infil,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- dt = %f [ %d ]",dt,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- vpd = %f [ %d ]",vpd,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- snow_energy = %f [ %d ]",snow_energy,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- mu = %f [ %d ]",mu,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- rainfall = %f [ %d ]",rainfall[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Wdew = %f [ %d ]",Wdew[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- grnd_flux = %f [ %d ]",grnd_flux[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- T1 = %f [ %d ]",T1[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- latent_heat = %f [ %d ]",latent_heat[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- sensible_heat = %f [ %d ]",sensible_heat[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- deltaH = %f [ %d ]",deltaH[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- store_error = %f [ %d ]",store_error[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- depth = %f [ %d ]",depth[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Wcr = %f [ %d ]",Wcr[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Wpwp = %f [ %d ]",Wpwp[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- Residual Moisture = %f [ %d ]",resid_moist[0],rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- VEG = %i [ %d ]",VEG,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- veg_class = %i [ %d ]",veg_class,rec);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- month = %i [ %d ]",month,rec);
  lis_log_msgC(lismsg);
#if 0
  write_layer(layer_wet,iveg,options.Nlayer,depth);
  if(options.DIST_PRCP) 
    write_layer(layer_dry,iveg,options.Nlayer,depth);
  write_vegvar(&(veg_var_wet[0]),iveg);
  if(options.DIST_PRCP) 
    write_vegvar(&(veg_var_dry[0]),iveg);

  if(!options.QUICK_FLUX) {
    fprintf(stderr,"Node\tT\tTnew\tdz\tkappa\tCs\tmoist\tbubble\texpt\tmax_moist\tice\n");
    for(i=0;i<Nnodes;i++) 
      fprintf(stderr,"%i\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
	      i, T_node[i], Tnew_node[i], dz_node[i], kappa_node[i], 
	      Cs_node[i], moist_node[i], bubble_node[i], expt_node[i], 
	      max_moist_node[i], ice_node[i]);
  }

  vicerror("Finished writing calc_surf_energy_bal variables.\nTry increasing SURF_DT to get model to complete cell.\nThen check output for instabilities.");
#endif

  return(0.0);
    
}
#endif
