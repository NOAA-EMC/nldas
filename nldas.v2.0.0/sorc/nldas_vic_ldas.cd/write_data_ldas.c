#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

#include "vicNl_ldas.h"

static char vcid[] = "$Id: write_data.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void write_data_ldas(out_data_struct *out_data,
/** Start of changes - NGMOD **/
		     /* ldas_outfiles_struct *outfiles, */
		     ldas_outfiles_struct *outfiles,
/** End of changes - NGMOD **/
		     dmy_struct      *dmy,
		     int              dt)
/**********************************************************************
	write_data	Dag Lohmann		Janurary 1996

  This subroutine writes all energy and moisture balance parameters to
  output files.

  OUTPUT:
	evaporation and vapor fluxes in mm/time step
	layer moisture in mm/time step
	runoff in mm/time step
	baseflow in mm/time step
	freezing and thawing depths in cm
	snow depth in cm
	snow water equivlence in mm
	all energy fluxes are in W/m^2

  Modifications:
  5/20/96	Program was modified to account for a variable
		number of soil layers.  It was also modified to
		write out frozen soils data per time step.	KAC
  1/15/97	Program modified to output daily sums, or values
		independant of selected time step.  This aids in
		comparisons between model versions.		KAC
  3/98          Routine modified to output fluxes in PILPS2c 
                ASCII column format                             Dag
  4/30/98       Routine modified to add binary output options for
                improved file speed, and less disk usage for large
		model basins                                    KAC
  7/19/99       modified to output a single binary file containing
                the data selected for the LDAS project         KAC
  8/3/99        modified again to reduce the storage space needed
                for the LDAS output files.  
  1/4/2000      modified to allow both standard and LDAS formatted
                output using a compiler flag                    KAC

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  char               *tmp_cptr;
  short int          *tmp_siptr;
  unsigned short int *tmp_usiptr;
  int                 band, j;
  int                 fidx;
  int                *tmp_iptr;
  float              *tmp_fptr;

  /***************************************************************


  ***************************************************************/

  tmp_iptr = (int *)  calloc(1,sizeof(int));
  tmp_fptr = (float *)calloc(1,sizeof(float));

  /**************************************
    Ouput Snow Variables
  **************************************/
#if 0
  if(options.BINARY_OUTPUT) {
    /***** Write Binary Snow Output File *****/
    tmp_fptr[0] = (float)out_data->swq[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SWE);
    tmp_fptr[0] = (float)out_data->snow_depth[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SnowDepth);
  }
  else{
    /***** Write ASCII full energy snow output file *****/
    fprintf(outfiles->SWE,"%.4f\n",out_data->swq[0]);
    fprintf(outfiles->SnowDepth,"%.4f\n",out_data->snow_depth[0]);
  }
#endif


  /************************************
    Output Standard Energy and Moisture Flux Variables
  ************************************/

  if(options.BINARY_OUTPUT) {
  
    /***** Write Binary Fluxes File *****/
#if 1    
    /* energy balance components */
    tmp_fptr[0] = (float)out_data->net_short;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SWnet);
    tmp_fptr[0] = (float)out_data->net_long;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->LWnet);
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      tmp_fptr[0] = (float)out_data->latent;
      fwrite(tmp_fptr,1,sizeof(float),outfiles->Qle);
      tmp_fptr[0] = (float)out_data->sensible;
      fwrite(tmp_fptr,1,sizeof(float),outfiles->Qh);
      tmp_fptr[0] = (float)out_data->grnd_flux;
      fwrite(tmp_fptr,1,sizeof(float),outfiles->Qg);
      tmp_fptr[0] = (float)out_data->snow_fluxes[0];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->Qf);
    }
    tmp_fptr[0] = (float)out_data->shortwave;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SWdown);
    tmp_fptr[0] = (float)out_data->in_long;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->LWdown);
  
    /* water balance components */
    tmp_fptr[0] = (float)out_data->snowfall;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Snowf);
    tmp_fptr[0] = (float)out_data->rainfall;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Rainf);
    tmp_fptr[0] = (float)out_data->evap;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Evap);
    tmp_fptr[0] = (float)out_data->runoff;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Qs);
    tmp_fptr[0] = (float)out_data->baseflow;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Qsb);
    tmp_fptr[0] = (float)out_data->snow_melt[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Qsm);

    /* surface state variables */
    tmp_fptr[0] = (float)out_data->snow_pack_temp[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SnowT);
    tmp_fptr[0] = (float)out_data->surf_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->AvgSurfT);
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      tmp_fptr[0] = (float)out_data->rad_temp;
      fwrite(tmp_fptr,1,sizeof(float),outfiles->RadT);
    }
    tmp_fptr[0] = (float)out_data->albedo;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Albedo);
    tmp_fptr[0] = (float)out_data->swq[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SWE);
    tmp_fptr[0] = (float)out_data->Wdew;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->CanopInt);
    
    /* subsurface state variables */
    for(j=0;j<options.Nlayer;j++) {    
      tmp_fptr[0] = (float)out_data->soil_temp[j];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->SoilTemp[j]);
    }
    tmp_fptr[0] = (float)out_data->moist_total;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SoilMoistTotal);
    tmp_fptr[0] = (float)out_data->moist_root;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SoilMoistRoot);
    tmp_fptr[0] = (float)out_data->moist_1m;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SoilMoist1m);
    for(j=0;j<options.Nlayer;j++) {    
      tmp_fptr[0] = (float)out_data->moist[j];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->SoilMoist[j]);
      tmp_fptr[0] = (float)out_data->liquid[j];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->LSoilMoist[j]);
    }
    tmp_fptr[0] = (float)out_data->wetness_total;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SoilWetTotal);
    tmp_fptr[0] = (float)out_data->wetness_root;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SoilWetRoot);
    
    /* evaporation components */
    tmp_fptr[0] = (float)out_data->evap_canop;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->ECanop);
    tmp_fptr[0] = (float)out_data->evap_veg;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->TVeg);
    tmp_fptr[0] = (float)out_data->evap_bare;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->ESoil);
    tmp_fptr[0] = (float)out_data->sub_snow + (float)out_data->sub_canop; /* sub = pack + canopy */
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SubSnow);    
    tmp_fptr[0] = (float)out_data->aero_conduct;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->ACond);
    tmp_fptr[0] = (float)out_data->lai;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->LAI);
    
    /* cold season processes */
    tmp_fptr[0] = (float)out_data->snow_depth[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SnowDepth);
    tmp_fptr[0] = (float)out_data->coverage[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SnowFrac);
    tmp_fptr[0] = (float)out_data->snow_albedo[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SAlbedo);
    
    /* non-LDAS, but useful */
    tmp_fptr[0] = (float)out_data->prec;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Precip);
    tmp_fptr[0] = (float)out_data->r_net;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->RNet);
    tmp_fptr[0] = (float)out_data->air_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Air_Temp);
    tmp_fptr[0] = (float)out_data->wind;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Wind);
    tmp_fptr[0] = (float)out_data->rel_humid;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->RHumid);

#else
    /***** Write Binary Fluxes File *****/
    tmp_fptr[0] = (float)out_data->prec;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Precip);
    tmp_fptr[0] = (float)out_data->evap;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Evap);
    tmp_fptr[0] = (float)out_data->runoff;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Qs);
    tmp_fptr[0] = (float)out_data->baseflow;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Qsb);
    for(j=0;j<options.Nlayer;j++) {    
      tmp_fptr[0] = (float)out_data->moist[j];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->SoilMoist[j]);
    }
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      tmp_fptr[0] = (float)out_data->rad_temp;
      fwrite(tmp_fptr,1,sizeof(float),outfiles->RadT);
    }
    tmp_fptr[0] = (float)out_data->net_short;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->SWnet);
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      tmp_fptr[0] = (float)out_data->latent;
      fwrite(tmp_fptr,1,sizeof(float),outfiles->Qle);
    }
    tmp_fptr[0] = (float)out_data->evap_canop;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->ECanop);
    tmp_fptr[0] = (float)out_data->evap_veg;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->TVeg);
    tmp_fptr[0] = (float)out_data->evap_bare;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->ESoil);
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      tmp_fptr[0] = (float)out_data->sensible;
      fwrite(tmp_fptr,1,sizeof(float),outfiles->Qh);
      tmp_fptr[0] = (float)out_data->grnd_flux;
      fwrite(tmp_fptr,1,sizeof(float),outfiles->Qg);
    }
    tmp_fptr[0] = (float)out_data->aero_resist;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->ACond);
    tmp_fptr[0] = (float)out_data->surf_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->AvgSurfT);
    tmp_fptr[0] = (float)out_data->albedo;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Albedo);

    tmp_fptr[0] = (float)out_data->r_net;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->RNet);
    tmp_fptr[0] = (float)out_data->in_long;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->InLong);
    tmp_fptr[0] = (float)out_data->net_long;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->NetLong);
    tmp_fptr[0] = (float)out_data->air_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Air_Temp);
    tmp_fptr[0] = (float)out_data->wind;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->Wind);
    tmp_fptr[0] = (float)out_data->rel_humid;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->RHumid);
#endif

  }
  else {
    /***** Write ASCII energy balance fluxes file *****/

#if 1
    /* energy balance components */
    fprintf(outfiles->SWnet,"%.4f\n",out_data->net_short);
    fprintf(outfiles->LWnet,"%.4f\n",out_data->net_long);
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      fprintf(outfiles->Qle,"%.4f\n",out_data->latent);
      fprintf(outfiles->Qh,"%.4f\n",out_data->sensible);
      fprintf(outfiles->Qg,"%.4f\n",out_data->grnd_flux);
      fprintf(outfiles->Qf,"%.4f\n",out_data->snow_fluxes[0]);
    }
    fprintf(outfiles->SWdown,"%.4f\n",out_data->shortwave);
    fprintf(outfiles->LWdown,"%.4f\n",out_data->in_long);

    /* water balance components */
    fprintf(outfiles->Snowf,"%.4f\n",out_data->snowfall);
    fprintf(outfiles->Rainf,"%.4f\n",out_data->rainfall);
    fprintf(outfiles->Evap,"%.4f\n",out_data->evap);
    fprintf(outfiles->Qs,"%.4f\n",out_data->runoff);
    fprintf(outfiles->Qsb,"%.4f\n",out_data->baseflow);
    fprintf(outfiles->Qsm,"%.4f\n",out_data->snow_melt[0]);

    /* surface state variables */
    fprintf(outfiles->SnowT,"%.4f\n",out_data->snow_pack_temp[0]);
    fprintf(outfiles->AvgSurfT,"%.4f\n",out_data->surf_temp);
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      fprintf(outfiles->RadT,"%.4f\n",out_data->rad_temp);
    }
    fprintf(outfiles->Albedo,"%.4f\n",out_data->albedo);
    fprintf(outfiles->SWE,"%.4f\n",out_data->swq[0]);
    fprintf(outfiles->CanopInt,"%.4f\n",out_data->Wdew);

    /* subsurface state variables */
    for(j=0;j<options.Nlayer;j++) {    
      fprintf(outfiles->SoilTemp[j],"%.4f\n",out_data->soil_temp[j]);
    }
    fprintf(outfiles->SoilMoistTotal,"%.4f\n",out_data->moist_total);
    fprintf(outfiles->SoilMoistRoot,"%.4f\n",out_data->moist_root);
    fprintf(outfiles->SoilMoist1m,"%.4f\n",out_data->moist_1m);
    for(j=0;j<options.Nlayer;j++) {    
      fprintf(outfiles->SoilMoist[j],"%.4f\n",out_data->moist[j]);
      fprintf(outfiles->LSoilMoist[j],"%.4f\n",out_data->liquid[j]);
    }
    fprintf(outfiles->SoilWetTotal,"%.4f\n",out_data->wetness_total);
    fprintf(outfiles->SoilWetRoot,"%.4f\n",out_data->wetness_root);
    
    /* evaporation components */
    fprintf(outfiles->ECanop,"%.4f\n",out_data->evap_canop);
    fprintf(outfiles->TVeg,"%.4f\n",out_data->evap_veg);
    fprintf(outfiles->ESoil,"%.4f\n",out_data->evap_bare);
    fprintf(outfiles->SubSnow,"%.4f\n",out_data->sub_snow + out_data->sub_canop); /* sub = pack + canopy */
    fprintf(outfiles->ACond,"%.4f\n",out_data->aero_conduct);
    fprintf(outfiles->LAI,"%.4f\n",out_data->lai);
    
    /* cold season processes */
    fprintf(outfiles->SnowDepth,"%.4f\n",out_data->snow_depth[0]);
    fprintf(outfiles->SnowFrac,"%.4f\n",out_data->coverage[0]);
    fprintf(outfiles->SAlbedo,"%.4f\n",out_data->snow_albedo[0]);

    /* non-LDAS, but useful */
    fprintf(outfiles->Precip,"%.4f\n",out_data->prec);
    fprintf(outfiles->RNet,"%.4f\n",out_data->r_net);
    fprintf(outfiles->Air_Temp,"%.4f\n",out_data->air_temp);
    fprintf(outfiles->Wind,"%.4f\n",out_data->wind);
    fprintf(outfiles->RHumid,"%.4f\n",out_data->rel_humid);
    
#else   
    fprintf(outfiles->Precip,"%.4f\n",out_data->prec);
    fprintf(outfiles->Evap,"%.4f\n",out_data->evap);
    fprintf(outfiles->Qs,"%.4f\n",out_data->runoff);
    fprintf(outfiles->Qsb,"%.4f\n",out_data->baseflow);
    for(j=0;j<options.Nlayer;j++) {    
      fprintf(outfiles->SoilMoist[j],"%.4f\n",out_data->moist[j]);
    }
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      fprintf(outfiles->RadT,"%.4f\n",out_data->rad_temp);
    }
    fprintf(outfiles->SWnet,"%.4f\n",out_data->net_short);
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      fprintf(outfiles->Qle,"%.4f\n",out_data->latent);
    }
    fprintf(outfiles->ECanop,"%.4f\n",out_data->evap_canop);
    fprintf(outfiles->TVeg,"%.4f\n",out_data->evap_veg);
    fprintf(outfiles->ESoil,"%.4f\n",out_data->evap_bare);
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      fprintf(outfiles->Qh,"%.4f\n",out_data->sensible);
      fprintf(outfiles->Qg,"%.4f\n",out_data->grnd_flux);
    }
    fprintf(outfiles->ACond,"%.4f\n",out_data->aero_resist);
    fprintf(outfiles->AvgSurfT,"%.4f\n",out_data->surf_temp);
    fprintf(outfiles->Albedo,"%.4f\n",out_data->albedo);

    fprintf(outfiles->RNet,"%.4f\n",out_data->r_net);
    fprintf(outfiles->InLong,"%.4f\n",out_data->in_long);
    fprintf(outfiles->NetLong,"%.4f\n",out_data->net_long);
    fprintf(outfiles->Air_Temp,"%.4f\n",out_data->air_temp);
    fprintf(outfiles->Wind,"%.4f\n",out_data->wind);
    fprintf(outfiles->RHumid,"%.4f\n",out_data->rel_humid);

#endif
 
  }
  

  free((char *)tmp_iptr);
  free((char *)tmp_fptr);


}

void calc_water_balance_error(int    rec,
			      double inflow,
			      double outflow,
			      double storage) {
/***************************************************************
  calc_water_balance_error  Keith Cherkauer        April 1998

  This subroutine computes the overall model water balance, and 
  warns the model user if large errors are found.
***************************************************************/

  static double last_storage;
  static double cum_error;
  static double max_error;
  static int    error_cnt;
  static int    Nrecs;

  double error;

  if(rec<0) {
    last_storage = storage;
    cum_error    = 0.;
    max_error    = 0.;
    error_cnt    = 0;
    Nrecs        = -rec;
  }
  else {
    error = inflow - outflow - (storage - last_storage);
    cum_error += error;
    if(fabs(error)>fabs(max_error) && fabs(error)>1e-5) {
      max_error = error;
      fprintf(stderr,"Maximum Moist Error:\t%i\t%.5f\t%.5f\n",
	      rec,error,cum_error);
    }
    if(rec==Nrecs-1) {
      if(cum_error>1e-5)
	fprintf(stderr,"Total Cumulative Water Error for Grid Cell = %.4f\n",
		cum_error);
    }
    last_storage = storage;
  }
}

void calc_energy_balance_error(int    rec,
			       double net_rad,
			       double latent,
			       double sensible,
			       double grnd_flux,
			       double snow_fluxes) {
/***************************************************************
  calc_energy_balance_error   Keith Cherkauer     April 1998

  This subroutine computes the overall model energy balance, and
  reports the maximum time step error above a thresehold to the
  user.  The total cumulative error for the grid cell is also 
  computed and reported at the end of the model run.
***************************************************************/

  static double cum_error;
  static double max_error;
  static int    Nrecs;

  double error;

  if(rec<0) {
    cum_error = 0;
    Nrecs     = -rec;
    max_error = 0;
  }
  else {
    error = net_rad - latent - sensible - grnd_flux + snow_fluxes;
    cum_error += error;
    if(fabs(error)>fabs(max_error) && fabs(error)>0.001) {
      max_error = error;
      if ( rec > 0 ) 
	fprintf(stderr,"Maximum Energy Error:\t%i\t%.4f\t%.4f\n",
		rec,error,cum_error/(double)rec);
      else 
	fprintf(stderr,"Maximum Energy Error:\t%i\t%.4f\t%.4f\n",
		rec,error,cum_error);
    }
    if(rec==Nrecs-1) {
      if(cum_error>0.001)
	fprintf(stderr,"Total Cumulative Energy Error for Grid Cell = %.4f\n",
		cum_error/(double)rec);
    }
  }
}

