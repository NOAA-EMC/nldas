#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

#include "vicNl_ldas.h"

static char vcid[] = "$Id: read_initial_model_state.c,v 4.2 2000/05/16 21:57:54 vicadmin Exp vicadmin $";

void read_initial_model_state_ldas(FILE                *statefile,
				   dist_prcp_struct    *prcp,
				   global_param_struct *gp,
				   int                  Nveg,
				   int                  Nbands,
				   int                  cellnum,
				   soil_con_struct     *soil_con,
				   veg_con_struct      *veg_con,
				   int                 *dry_time
/** Start of changes -- JS **/
                          ,int                *still_storm
/** End of changes -- JS **/
				   ) 
/*********************************************************************
  read_initial_model_state   Keith Cherkauer         April 14, 2000

  This subroutine initializes the model state at hour 0 of the date 
  defined in the given state file.  

  Soil moisture, soil thermal, and snowpack variables  are stored 
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

  LDAS MODIFICATIONS
  Modified so mu is saved, and the storages are not merged.

  Modifications:
	Sep 2002 JS: added pack_water and surf_water to state file

*********************************************************************/
{
  extern option_struct options;

  char   tmpstr[MAXSTRING];
  char   ErrStr[MAXSTRING];
  double tmpval;
  double Nsum;
  double depth_node[MAX_NODES];
  int    veg, iveg;
  int    band, iband;
  int    lidx;
  int    nidx;
  int    dist;
  int    tmp_cellnum;
  int    tmp_Nveg;
  int    tmp_Nband;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;

  /* skip this stuff */
  /* rewind(statefile); */
  /* skip header */
  /*
  fgets(tmpstr, MAXSTRING, statefile);
  fgets(tmpstr, MAXSTRING, statefile);
  */

  /* read cell information */
/** Start of changes -- JS **/
/*  fscanf(statefile, "%i %i %i %d", 
	 &tmp_cellnum, &tmp_Nveg, &tmp_Nband, dry_time);*/
  fscanf(statefile, "%i %i %i %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband, dry_time, still_storm);
/** End of changes -- JS **/

    /*fprintf(stderr, ">>>>%i %i %i %d %d\n", 
	 tmp_cellnum, tmp_Nveg, tmp_Nband, *dry_time, *still_storm);
	 fflush(NULL); */
  
  /* read the veg file to get 'current cell */
  while ( tmp_cellnum != cellnum && !feof(statefile) ) {
    /* skip over the rest of the string */
    fgets(tmpstr, MAXSTRING, statefile);
    /*fprintf(stderr, "%s\n", tmpstr);*/
    for ( veg = 0; veg < tmp_Nveg; veg++ ) 
      for ( band = 0; band < tmp_Nband; band++ ) {
	fgets(tmpstr,MAXSTRING,statefile);
	/*fprintf(stderr, "%d %s\n", veg, tmpstr);*/
      }
/** Start of changes -- JS **/
/* need to read in dry_time as well */
/*    fscanf(statefile, "%i %i %i", &tmp_cellnum, &tmp_Nveg, 
	   &tmp_Nband);*/
    fscanf(statefile, "%i %i %i %d %d", &tmp_cellnum, &tmp_Nveg, 
	   &tmp_Nband, dry_time, still_storm);
    /*fprintf(stderr, ">>>>%i %i %i %d %d\n", 
	 tmp_cellnum, tmp_Nveg, tmp_Nband, *dry_time, *still_storm);*/
	 fflush(NULL); 
/** End of changes -- JS **/
  }


  /* could not find cell */
  if ( feof(statefile) ) {
    sprintf(ErrStr, "Requested grid cell (%i) is not in the model state file.", 
	    cellnum);
    nrerror(ErrStr);
  }

  /* Read soil thermal node depths */
  for ( nidx = 0; nidx < options.Nnode; nidx++ ){
    fscanf(statefile," %lf", &depth_node[nidx]);
  }

  if ( options.Nnode > 1 )
    compute_dz(soil_con->dz_node, depth_node, options.Nnode, soil_con->dp);
  else soil_con->dz_node[0] = 0;
    
  /* Input for all vegetation types */
  for ( veg = 0; veg < tmp_Nveg; veg++ ) {
    
    /* Input for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {
      
      /* Check if band computational */
      if(soil_con->AreaFract[band] > 0){

	/* read in mu */
	fscanf(statefile, "%i %i %lf", &veg, &iband, &prcp->mu[veg]);

	/* Read WET and DRY soil moisture */
	for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	  if ( fscanf(statefile," %lf", 
		      &cell[WET][veg][band].layer[lidx].moist) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
	  if ( fscanf(statefile," %lf", 
		      &cell[DRY][veg][band].layer[lidx].moist) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
/** Start of changes -- JS **/
	  if ( fscanf(statefile," %lf", 
		      &cell[WET][veg][band].layer[lidx].ice) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
	  if ( fscanf(statefile," %lf", 
		      &cell[DRY][veg][band].layer[lidx].ice) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
/** End of changes -- JS **/ 
	}
	
	/* Read average dew storage */
	if ( veg < tmp_Nveg-1 || veg_con[0].Cv_sum == 1 ) {
	  if ( fscanf(statefile," %lf", &veg_var[WET][veg][band].Wdew) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
	  if ( fscanf(statefile," %lf", &veg_var[DRY][veg][band].Wdew) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
	}
	
	/* Read snow data */
/** Start of changes -- JS **/
/*
	if ( fscanf(statefile," %i %lf %lf %lf %lf %lf", 
		    &snow[veg][band].last_snow, &snow[veg][band].swq, 
		    &snow[veg][band].surf_temp, &snow[veg][band].pack_temp, 
		    &snow[veg][band].density, &snow[veg][band].snow_canopy) 
*/
        /* add pack_water and surf_water to state file */
 	if ( fscanf(statefile," %i %lf %lf %lf %lf %lf %lf %lf", 
		    &snow[veg][band].last_snow, &snow[veg][band].swq, 
		    &snow[veg][band].surf_temp, &snow[veg][band].pack_temp, 
  		    &snow[veg][band].density, &snow[veg][band].snow_canopy,
  		    &snow[veg][band].pack_water, &snow[veg][band].surf_water) 
/** End of change -- JS **/
	     == EOF ) 
	  nrerror("End of model state file found unexpectedly");
	if(snow[veg][band].density > 0.) 
	  snow[veg][band].depth = 1000. * snow[veg][band].swq 
	    / snow[veg][band].density;
	
	/* Read soil thermal node temperatures */
	for ( nidx = 0; nidx < options.Nnode; nidx++ ) 
	  if ( fscanf(statefile," %lf", &energy[veg][band].T[nidx]) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
	
      }
    }
  }
  /** Start of changes -- JS **/
  /* make sure we skip over the rest of the string */
  fgets(tmpstr, MAXSTRING, statefile);
  /** Start of changes -- JS **/

}
