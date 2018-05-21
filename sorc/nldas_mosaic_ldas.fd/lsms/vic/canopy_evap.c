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


float canopy_evap(layer_data_struct *layer_wet,
                   layer_data_struct *layer_dry,
                   veg_var_struct    *veg_var_wet, 
                   veg_var_struct    *veg_var_dry, 
		   char               CALC_EVAP,
                   int                veg_class, 
                   int                month, 
                   float             mu,
		   float            *Wdew,
                   float             dt,
                   float             rad,
		   float             vpd,
		   float             net_short,
		   float             air_temp,
                   float             ra,
                   float             displacement,
                   float             roughness,
                   float             ref_height,
		   float             elevation,
                   float            *prec,
		   float            *depth,
		   float            *Wcr,
		   float            *Wpwp,
		   float             *root,
		   int               Nlayer)
/**********************************************************************
	canopy_evap.c	Dag Lohmann		September 1995

  This routine computes the evaporation, traspiration and throughfall
  of the vegetation types for multi-layered model.

  The value of x, the fraction of precipitation that exceeds the 
  canopy storage capacity, is returned by the subroutine.

  UNITS:	moist (mm)
		evap (mm)
		prec (mm)
		melt (mm)

  VARIABLE TYPE        NAME          UNITS DESCRIPTION
  atmos_data_struct    atmos         N/A   atmospheric forcing data structure
  layer_data_struct   *layer         N/A   soil layer variable structure
  veg_var_struct      *veg_var       N/A   vegetation variable structure
  soil_con_struct      soil_con      N/A   soil parameter structure
  char                 CALC_EVAP     N/A   TRUE = calculate evapotranspiration
  int                  veg_class     N/A   vegetation class index number
  int                  month         N/A   current month
  global_param_struct  global        N/A   global parameter structure
  float               mu            fract wet (or dry) fraction of grid cell
  float               ra            s/m   aerodynamic resistance
  float               prec          mm    precipitation
  float               displacement  m     displacement height of surface cover
  float               roughness     m     roughness height of surface cover
  float               ref_height    m     measurement reference height

  Modifications:
  9/1/97	Greg O'Donnell
  4-12-98  Code cleaned and final version prepared, KAC
  06-25-98 modified for new distributed precipitation data structure KAC
  01-19-00 modified to function with new simplified soil moisture 
           scheme                                                  KAC

**********************************************************************/
{

  /** declare global variables **/
  extern veg_lib_struct *veg_lib; 

  /** declare local variables **/
  int                Ndist;
  int                dist;
  int                i;
  float             ppt;		/* effective precipitation */
  float             f;		/* fraction of time step used to fill canopy */
  float             throughfall;
  float             Evap;
  float             tmp_Evap;
  float             canopyevap;
  float             tmp_Wdew;
  float             layerevap[MAX_LAYERS];
  layer_data_struct *tmp_layer;
  veg_var_struct    *tmp_veg_var;

  /********************************************************************** 
     CANOPY EVAPORATION

     Calculation of evaporation from the canopy, including the
     possibility of potential evaporation exhausting ppt+canopy storage
     2.16 + 2.17
     Index [0] refers to current time step, index [1] to next one
     If f < 1.0 than veg_var->canopyevap = veg_var->Wdew + ppt and
                     Wdew = 0.0

     DEFINITIONS:
     Wdmax - max monthly dew holding capacity
     Wdew - dew trapped on vegetation

     Modified 
     04-14-98 to work within calc_surf_energy_balance.c  KAC
     07-24-98 fixed problem that caused hourly precipitation
              to evaporate from the canopy during the same
	      time step that it falls (OK for daily time step, 
	      but causes oscilations in surface temperature
	      for hourly time step)                      KAC, Dag
	      

  **********************************************************************/ 

  //  if(options.DIST_PRCP) Ndist = 2;
  //else 
  Ndist = 1;

  Evap = 0;
  for(dist=0;dist<Ndist;dist++) {

    /* Initialize variables */
    for(i=0;i<Nlayer;i++) layerevap[i] = 0;
    canopyevap = 0;
    throughfall = 0;

    /* Set parameters for distributed precipitation */
    if(dist==0) {
      tmp_layer   = layer_wet;
      tmp_veg_var = veg_var_wet;
      ppt         = prec[WET];
      tmp_Wdew    = Wdew[WET];
    }
    else {
      tmp_layer   = layer_dry;
      tmp_veg_var = veg_var_dry;
      ppt         = prec[DRY];
      mu          = (1. - mu);
      tmp_Wdew    = Wdew[DRY];
    }      

    if(mu > 0) {

      /****************************************************
        Compute Evaporation from Canopy Intercepted Water
      ****************************************************/

      /** Due to month changes ..... Wdmax based on LAI **/
      tmp_veg_var->Wdew = tmp_Wdew;
      if (tmp_Wdew > veg_lib[veg_class].Wdmax[month-1]) {
	throughfall = tmp_Wdew - veg_lib[veg_class].Wdmax[month-1];
	tmp_Wdew    = veg_lib[veg_class].Wdmax[month-1];
      }
      
      canopyevap = pow((tmp_Wdew / veg_lib[veg_class].Wdmax[month-1]),
		       (2.0/3.0))* vic_penman(rad, vpd * 1000., ra, (float) 0.0, 
					  veg_lib[veg_class].rarc,
					  veg_lib[veg_class].LAI[month-1], 
					  (float) 1.0, air_temp, 
					  net_short, elevation, 
#if TIMESTEPSECS
                                          veg_lib[veg_class].RGL) * dt / (24.0 * 3600.0);
#else					    
					  veg_lib[veg_class].RGL) * dt / 24.0;
#endif

#if TIMESTEPSECS
      if (canopyevap > 0.0 && dt==(24.0 *3600.0))
#else
      if (canopyevap > 0.0 && dt==24)
#endif
	/** If daily time step, evap can include current precipitation **/
	f = min(1.0,((tmp_Wdew + ppt) / canopyevap));
      else if (canopyevap > 0.0)
	/** If sub-daily time step, evap can not exceed current storage **/
	f = min(1.0,((tmp_Wdew) / canopyevap));
      else
	f = 1.0;
      canopyevap *= f;
    
      tmp_Wdew += ppt - canopyevap;
      if (tmp_Wdew < 0.0) 
	tmp_Wdew = 0.0;
      if (tmp_Wdew <= veg_lib[veg_class].Wdmax[month-1]) 
	throughfall += 0.0;
      else {
	throughfall += tmp_Wdew - veg_lib[veg_class].Wdmax[month-1];
	tmp_Wdew = veg_lib[veg_class].Wdmax[month-1];
      }

      /*******************************************
        Compute Evapotranspiration from Vegetation
      *******************************************/
      if(CALC_EVAP)
	transpiration(tmp_layer, veg_class, month, rad,
		      vpd, net_short, air_temp, ra,
		      ppt, f, dt, tmp_veg_var->Wdew, elevation,
		      depth, Wcr, Wpwp, &tmp_Wdew,
		      &canopyevap, layerevap, root,Nlayer);

    }

    tmp_veg_var->canopyevap = canopyevap;
    tmp_veg_var->throughfall = throughfall;
    tmp_veg_var->Wdew = tmp_Wdew;
    tmp_Evap = canopyevap;
    for(i=0;i<Nlayer;i++) {
      tmp_layer[i].evap  = layerevap[i];
      tmp_Evap          += layerevap[i];

    }
    
    /* Evap += tmp_Evap * mu / (1000. * dt * 3600.); */
    Evap += tmp_Evap * mu / ( 1000. * dt );


  }
  return (Evap);

}

/**********************************************************************
	EVAPOTRANSPIRATION ROUTINE
**********************************************************************/

void transpiration(layer_data_struct *layer,
		   int veg_class, 
		   int month, 
		   float rad,
		   float vpd,
		   float net_short,
		   float air_temp,
		   float ra,
		   float ppt,
		   float f,
		   float dt,
		   float Wdew,
		   float elevation,
		   float *depth,
		   float *Wcr,
		   float *Wpwp,
		   float *new_Wdew,
		   float *canopyevap,
		   float *layerevap,
		   float  *root,
		   int Nlayer)
/**********************************************************************
  Computes evapotranspiration for unfrozen soils
  Allows for multiple layers.
**********************************************************************/
{
  extern veg_lib_struct *veg_lib;

  int i;
  float gsm_inv;               	/* soil moisture stress factor */
  float moist1, moist2;                /* tmp holding of moisture */
  float evap;                          /* tmp holding for evap total */
  float Wcr1;                          /* tmp holding of critical water for upper layers */
  float root_sum;                      /* proportion of roots in moist>Wcr zones */
  float spare_evap;                    /* evap for 2nd distribution */
  float avail_moist[MAX_LAYERS];         /* moisture available for trans */

  /********************************************************************** 
     EVAPOTRANSPIRATION

     Calculation of the evapotranspirations
     2.18

     First part: Soil moistures and root fractions of both layers
     influence each other

     Re-written to allow for multi-layers.
  **********************************************************************/
 
  /**************************************************
    Compute moisture content in combined upper layers
    **************************************************/
  moist1 = 0.0;
  Wcr1 = 0.0;  
  //for(i=0; i<Nlayer-1; i++)
    //printf("root %d %f\n",i,root[i]);

  for(i=0;i<Nlayer-1;i++){
    if(root[i] > 0.) {
      avail_moist[i] = layer[i].moist - layer[i].ice;

      moist1+=avail_moist[i];
      Wcr1 += Wcr[i];
    }
    else avail_moist[i]=0.;
  }

  /*****************************************
    Compute moisture content in lowest layer
    *****************************************/
  i=Nlayer-1;
  moist2 = layer[i].moist - layer[i].ice;

  avail_moist[i]=moist2;

  /******************************************************************
    CASE 1: Moisture in both layers exceeds Wcr, or Moisture in
    layer with more than half of the roots exceeds Wcr.

    Potential evapotranspiration not hindered by soil dryness.  If
    layer with less than half the roots is dryer than Wcr, extra
    evaporation is taken from the wetter layer.  Otherwise layers
    contribute to evapotransipration based on root fraction.
  ******************************************************************/

  if( (moist1>=Wcr1 && moist2>=Wcr[Nlayer-1] && Wcr1>0.) ||
      (moist1>=Wcr1 && (1-root[Nlayer-1])>= 0.5) ||
      (moist2>=Wcr[Nlayer-1] &&
      root[Nlayer-1]>=0.5) ){
    gsm_inv=1.0;
    //printf("vic pen %d %d %f\n",veg_class, month, veg_lib[veg_class].LAI[month-1]);
    evap = vic_penman(rad, vpd * 1000., ra, veg_lib[veg_class].rmin,
		  veg_lib[veg_class].rarc, veg_lib[veg_class].LAI[month-1], 
		  gsm_inv, air_temp, net_short, elevation, 
#if TIMESTEPSECS
     		  veg_lib[veg_class].RGL) * dt / (24.0*3600.0) *
#else
		  veg_lib[veg_class].RGL) * dt / 24.0 *
#endif
      (1.0-f*pow((Wdew/veg_lib[veg_class].Wdmax[month-1]),
		 (2.0/3.0)));
    /** divide up evap based on root distribution **/
    /** Note the indexing of the roots **/
    root_sum=1.0;
    spare_evap=0.0;
    for(i=0;i<Nlayer;i++){
      if(avail_moist[i]>=Wcr[i]){
        layerevap[i]=evap*(float)root[i];
	//printf("Layer1 %d %f %f %f \n",i, evap, root[i], layerevap[i]);
      }
      else {
          
        if (avail_moist[i] >= Wpwp[i]) 
          gsm_inv = (avail_moist[i] - Wpwp[i]) /
                    (Wcr[i] - Wpwp[i]);
        else 
          gsm_inv=0.0;
	    
        layerevap[i]  = evap*gsm_inv*(float)root[i];
	//printf("Layer2 %f %f %f \n",evap,root[i],layerevap[i]);
        root_sum     -= root[i];
        spare_evap    = evap*(float)root[i]*(1.0-gsm_inv);
      }
    }

    /** Assign excess evaporation to wetter layer **/
    if(spare_evap>0.0){
      for(i=0;i<Nlayer;i++){
        if(avail_moist[i] >= Wcr[i]){
          layerevap[i] += (float)root[i]*spare_evap/root_sum;
        }
      }
    }
  }

  /*********************************************************************
    CASE 2: Independent evapotranspirations

    Evapotranspiration is restricted by low soil moisture. Evaporation
    is computed independantly from each soil layer.
  *********************************************************************/

  else {

    for(i=0;i<Nlayer;i++){
      /** Set evaporation restriction factor **/
      if(avail_moist[i] >= Wcr[i])
	gsm_inv=1.0;
      else if(avail_moist[i] >= Wpwp[i])
	gsm_inv=(avail_moist[i] - Wpwp[i]) /
	  (Wcr[i] - Wpwp[i]);
      else 
	gsm_inv=0.0;

      if(gsm_inv > 0.0){
	/** Compute potential evapotranspiration **/
	//printf("bef Layer1 %d %f %f \n", i, root[i], layerevap[i]);
        layerevap[i] = vic_penman(rad, vpd * 1000., ra, veg_lib[veg_class].rmin,
			      veg_lib[veg_class].rarc, 
			      veg_lib[veg_class].LAI[month-1], gsm_inv, 
			      air_temp, net_short, elevation, 
#if TIMESTEPSECS
                              veg_lib[veg_class].RGL) * dt / (24.0 * 3600.0)
#else
			      veg_lib[veg_class].RGL) * dt / 24.0
#endif 
	  * (float)root[i] * (1.0-f*pow((Wdew/
					  veg_lib[veg_class].Wdmax[month-1]),
					 (2.0/3.0)));
	//printf("Layer1 %f %f \n", root[i], layerevap[i]);
      }
      else layerevap[i] = 0.0;

    }
  }
    
  /****************************************************************
    Check that evapotransipration does not cause soil moisture to 
    fall below wilting point.
  ****************************************************************/
  for(i=0;i<Nlayer;i++){
    //printf("Layer1 %f %f %f \n",layerevap[i],layer[i].moist, Wpwp[i]);
    if(layerevap[i] > layer[i].moist - Wpwp[i]) {
      layerevap[i] = layer[i].moist - Wpwp[i];
    }
    if ( layerevap[i] < 0.0 ) {
      layerevap[i] = 0.0;
    }
  }

}













