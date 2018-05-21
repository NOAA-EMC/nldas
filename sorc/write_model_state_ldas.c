#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

#include "vicNl_ldas.h"

static char vcid[] = "$Id: write_model_state.c,v 4.2 2000/05/16 21:57:54 vicadmin Exp vicadmin $";

#if SAVE_STATE

void write_model_state_ldas(dist_prcp_struct    *prcp,
			    global_param_struct *gp,
			    int                  Nveg,
			    int                  cellnum,
			    ldas_outfiles_struct *outfiles,
			    soil_con_struct      *soil_con,
			    veg_con_struct       *veg_con,
			    int                  DRY_TIME
/** Start of changes -- JS **/
                       ,int                 still_storm
/** End of changes -- JS **/
			    ) 
/*********************************************************************
  write_model_state      Keith Cherkauer           April 14, 2000

  This subroutine saves the model state at hour 0 of the date 
  defined in the global control file using STATEDAY, STATEMONTH,
  and STATEYEAR.  The saved files can then be used to initialize 
  the model to the same state as when the files were created.

  Soil moisture, soil thermal, and snowpack variables  are stored 
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

  LDAS MODIFICATIONS:
  Write out state for wet and dry fractions.
                                                       ngmod

  Modifications:
	Sep 2002 JS: added pack_water and surf_water to state file

*********************************************************************/
{
  extern option_struct options;

  double tmpval;
  double Nsum;
  int    veg;
  int    band;
  int    lidx;
  int    nidx;
  int    dist;
  int    Ndist;
  int    Nbands;

  int    active_bands;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;
  Nbands = options.SNOW_BAND;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;

  /* find number of active snow bands */
  active_bands=0;
  for ( band = 0; band < Nbands; band++ )
      if(soil_con->AreaFract[band] > 0)
	active_bands++;

  /* find if bare soil present */
  if(veg_con[0].Cv_sum == 1.)Nveg--;

  /* write cell information */
/** Start of changes -- JS **/
/*  fprintf(outfiles->statefile,"%i %i %i %d", 
	  cellnum, Nveg+1, active_bands, DRY_TIME );*/
  fprintf(outfiles->statefile,"%i %i %i %d %d", 
	  cellnum, Nveg+1, active_bands, DRY_TIME, still_storm );
/** End of changes -- JS **/

  /* Write soil thermal node depths */
  Nsum = 0;
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    fprintf(outfiles->statefile," %.15g", Nsum);
    if ( nidx < options.Nnode - 1 )
      Nsum += (soil_con->dz_node[nidx] + soil_con->dz_node[nidx+1]) / 2.;
  }

  fprintf(outfiles->statefile,"\n");

    /* Output for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {

    /* Output for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {
      
      /* check snow band has an area > 0 */
      if(soil_con->AreaFract[band] > 0){

	/* Write cell identification information and mu */
	/* mu may vary under different veg - ONLY if snow still present */
	fprintf(outfiles->statefile,"%i %i %.15g", veg, band, prcp->mu[veg]);
      
	/* Write average total soil moisture */
	for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	  fprintf(outfiles->statefile," %.15g ", 
		  cell[WET][veg][band].layer[lidx].moist);
	  fprintf(outfiles->statefile," %.15g ",
		  cell[DRY][veg][band].layer[lidx].moist);
/** Start of changes -- JS **/
	  fprintf(outfiles->statefile," %.15g ", 
		  cell[WET][veg][band].layer[lidx].ice);
	  fprintf(outfiles->statefile," %.15g ",
		  cell[DRY][veg][band].layer[lidx].ice);
/** End of changes -- JS **/
	}
	
	/* Write average dew storage */
	/* check that this is not bare soil */
	if ( veg_con[0].Cv_sum == 1 || veg < Nveg ) {
	  fprintf(outfiles->statefile," %.15g ", 
		  veg_var[WET][veg][band].Wdew);
	  fprintf(outfiles->statefile," %.15g ", 
		  veg_var[DRY][veg][band].Wdew);
	}
	
	/* Write snow data */
/** Start of changes -- JS **/
/*
	fprintf(outfiles->statefile," %i %.15g %.15g %.15g %.15g %.15g", 
		snow[veg][band].last_snow, snow[veg][band].swq, 
		snow[veg][band].surf_temp, snow[veg][band].pack_temp, 
		snow[veg][band].density, snow[veg][band].snow_canopy);
*/
        /* add pack_water and surf_water to state file */
 	fprintf(outfiles->statefile," %i %.15g %.15g %.15g %.15g %.15g %.15g %.15g", 
		snow[veg][band].last_snow, snow[veg][band].swq, 
		snow[veg][band].surf_temp, snow[veg][band].pack_temp, 
 		snow[veg][band].density, snow[veg][band].snow_canopy,
 		snow[veg][band].pack_water, snow[veg][band].surf_water);
/** End of changes -- JS **/
	
	/* Write soil thermal node temperatures */
	for ( nidx = 0; nidx < options.Nnode; nidx++ ) 
	  fprintf(outfiles->statefile," %.15g", energy[veg][band].T[nidx]);
	
	fprintf(outfiles->statefile,"\n");

      } /* end of band>0 loop */

    } /* end of band loop */

  } /* end of veg loop */

  /* Force file to be written */
  fflush(outfiles->statefile);

}

#endif
