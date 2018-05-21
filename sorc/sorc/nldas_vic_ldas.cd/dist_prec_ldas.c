#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vicNl.h"
/** Start of changes - NGMOD **/
#include "vicNl_ldas.h"
/** End of changes - NGMOD **/

static char vcid[] = "$Id: dist_prec.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void dist_prec_ldas(atmos_data_struct   *atmos,
		    dist_prcp_struct    *prcp,
		    soil_con_struct     *soil_con,
		    veg_con_struct      *veg_con,
		    dmy_struct          *dmy,
		    global_param_struct *global_param,
		    /** Start of changes - NGMOD **/
		    /* outfiles_struct     *outfiles, */
		    ldas_outfiles_struct *outfiles,
 		    /** End of changes - NGMOD **/
		    int                  rec,
		    int                  cellnum,
		    char                 NEWCELL,
		    char                 LASTREC,
		    /** Start of changes - NGMOD **/
		    /* add variable */
		    int                  statehour,
		    int                  dry_time,
		    int                 still_storm) {
 		    /** End of changes - NGMOD **/
/**********************************************************************
  dist_prec		Keith Cherkauer		October 9, 1997

  This subroutine calls the solution routines for a single grid cell
  for one time step.  It also controls the distribution of precipitation
  and reassembles grid cell data for output.

  The fractional coverage of precipitation over an area or grid cell, 
  mu, is estimated using the equation from Fan et. al. (1996).  The 
  coefficient, 0.6, was selected for the Arkansas - Red River Basin and
  was found using precipitation records on a 100km x 100km area.  It
  may not be applicable to all regions, please check the reference

  References:

  Modifications:
  11-30-98 Added counter to assure that a storm has been stopped
           for at least one day, before allowing the model to 
	   average soil moisture when a new precipitation event
	   arrives.                                             KAC

**********************************************************************/

  extern option_struct   options;
  extern veg_lib_struct *veg_lib; 
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  static char STILL_STORM;
  static int  DRY_TIME;

  char    ANY_SNOW;
  int     veg, i;
  int     month;
  double  Wdmax;
  double  NEW_MU;

  if(options.DIST_PRCP) {

    /*******************************************
      Controls Distributed Precipitation Model
    *******************************************/
	NEW_MU = 1.0 - exp(-options.PREC_EXPT*atmos->prec[NR]);
    for(veg=0; veg<=veg_con[0].vegetat_type_num; veg++) {
/*printf("TEMP rec = %d, veg = %d\n", rec, veg);*/
      ANY_SNOW = FALSE;
      for(i=0; i<options.SNOW_BAND; i++)
        /* Check for snow on ground or falling */
	if(prcp->snow[veg][i].swq > 0 
	   || prcp->snow[veg][i].snow_canopy > 0.) 
	  ANY_SNOW = TRUE;
      if(ANY_SNOW || atmos->snowflag[NR]) {
        /* If snow present, mu must be set to 1. */
	NEW_MU = 1.;
	if(rec == 0) {
          /* Set model variables if first time step */
	  /** Start of changes - NGMOD **/
	  if(options.INIT_STATE==FALSE){
	    DRY_TIME = 0;
	    prcp->mu[veg]=NEW_MU;
	  }
	  else{ /* INIT_STATE == TRUE */
	    DRY_TIME = dry_time;
	  }
	  if(atmos->prec[NR] > 0) 
	    STILL_STORM=TRUE;
	  else 
	    STILL_STORM=FALSE;
	} 
	/** End of changes - NGMOD **/
	ANY_SNOW = TRUE;
      }
      else {
      
	/** Start of changes - NGMOD **/
	/* change condition */
	
/*printf("TEMP Start, still_storm = %1d, rain = %7f, dry_time = %2d\n", STILL_STORM, atmos->prec[NR], DRY_TIME);*/

	if(rec==0 && options.INIT_STATE==TRUE) { /* INIT_STATE == TRUE */
/** Start of changes -- JS **/
	  /* DRY_TIME = dry_time; */
	  if(veg==0) {
	    DRY_TIME = dry_time;
	    if(still_storm)
	      STILL_STORM = TRUE;
	    else
	      STILL_STORM = FALSE;
	  }

/** End of changes -- JS **/
	}
	
	if(rec==0 && options.INIT_STATE==FALSE) {
	    if(atmos->prec[NR] == 0) {
	      /* If first time step has no rain, than set mu to 1. */
	      prcp->mu[veg] = 1.;
	      NEW_MU=1.;
	      STILL_STORM = TRUE;
	      DRY_TIME = 24;
	    }
	    else {
	      /* If first time step has rain, then set mu based on intensity */
	      prcp->mu[veg]=NEW_MU;
	      STILL_STORM=TRUE;
	      DRY_TIME = 0;
	    }
	}	
	
	/** End of changes - NGMOD **/
	else if(atmos->prec[NR] == 0 && DRY_TIME >= 24./(float)global_param->dt) {
          /* Check if storm has ended */
/*printf("TEMP 	no precip for 24hours\n");*/
	  NEW_MU=prcp->mu[veg];
	  STILL_STORM=FALSE;
          DRY_TIME = 0;
	}
	
        else if(atmos->prec[NR] == 0) {
	  /* May be pause in storm, keep track of pause length */
/*printf("TEMP 	pause in storm\n");*/
	  NEW_MU=prcp->mu[veg];
/** Start of changes -- JS **/
	  /*DRY_TIME += global_param->dt;*/
	  if(veg == 0) DRY_TIME += global_param->dt;
/** End of changes -- JS **/
	}
        else {
/*printf("TEMP 	we are in a storm\n");*/
        }
      }
      
/** Start of changes -- JS **/
      if(atmos->prec[NR] > STORM_THRES || ANY_SNOW) {
        DRY_TIME = 0;
      }
/** End of changes -- JS **/
      
      if(!STILL_STORM && (atmos->prec[NR] > STORM_THRES || ANY_SNOW)) {
	/** Average soil moisture before a new storm **/
/*printf("TEMP 	average moisture b4 storm start\n");*/
	initialize_new_storm(prcp->cell,prcp->veg_var,
			     veg,veg_con[0].vegetat_type_num,rec,
			     prcp->mu[veg],NEW_MU);
/** Start of changes -- JS **/
	/*STILL_STORM=TRUE;*/
	if(veg==veg_con[0].vegetat_type_num) STILL_STORM=TRUE;
/** End of changes -- JS **/
	prcp->mu[veg] = NEW_MU;
      }
      else if(NEW_MU != prcp->mu[veg] && STILL_STORM) {
	/** Redistribute soil moisture during the storm if mu changes **/
/*printf("TEMP 	redistr moisture during storm\n");*/
	if ( dmy[rec].day == 1 && dmy[rec].hour == 0 ) {
	  month = dmy[rec].month - 2;
	  if ( month < 0 ) month = 11;
	}
	else month = dmy[rec].month - 1;
	if (veg < veg_con[0].vegetat_type_num) 
	  Wdmax = veg_lib[veg_con[veg].veg_class].Wdmax[month];
	else 
	  Wdmax = 0;
	redistribute_during_storm(prcp->cell, prcp->veg_var, veg, 
				  veg_con[0].vegetat_type_num, rec, Wdmax, 
				  prcp->mu[veg], NEW_MU, soil_con->max_moist);
	prcp->mu[veg] = NEW_MU;
      }
/*printf("TEMP End,   still_storm = %1d, rain = %7f, dry_time = %2d\n", STILL_STORM, atmos->prec[NR], DRY_TIME);*/
    }

    /** Solve model time step **/
    full_energy(rec, atmos, soil_con, veg_con, prcp, dmy, global_param, 
		cellnum, NEWCELL);

  }

  else {

    /**************************************************
      Controls Grid Cell Averaged Precipitation Model
    **************************************************/

    full_energy(rec, atmos, soil_con, veg_con, prcp, dmy, global_param, 
		cellnum, NEWCELL);

  }

  /**************************************************
    Write cell average values for current time step
  **************************************************/

/** Start of changes - NGMOD **/
  /*function name change*/
    put_data_ldas(prcp, atmos, veg_con, outfiles, soil_con->depth, 
	   soil_con->dz_node, soil_con->dp, soil_con->AreaFract, 
	   soil_con->porosity,
	   &dmy[rec], rec, global_param->dt, options.Nnode, 
	   global_param->skipyear, 
/** Start of changes -- JS **/
	   global_param->MAX_SNOW_TEMP,
	   global_param->MIN_RAIN_TEMP); 
/** End of changes - JS **/
/** End of changes - NGMOD **/

    /* this has been moved to the end of the step */
#if SAVE_STATE

    /************************************
    Save model state at assigned date
    ************************************/
    /** Start of changes - NGMOD **/
    /*  fprintf(stderr,"%d %d\n", dmy[rec].hour, statehour); */
    if ( dmy[rec].year == global_param->stateyear
	 && dmy[rec].month == global_param->statemonth 
	 && dmy[rec].day == global_param->stateday 
	 && dmy[rec].hour == statehour ){
      write_model_state_ldas(prcp, global_param, veg_con[0].vegetat_type_num, 
			     soil_con->gridcel, outfiles, soil_con, veg_con, 
			   DRY_TIME
/** Start of changes -- JS **/
                       , (STILL_STORM==TRUE?1:0)
/** End of changes -- JS **/
			   );
    }
    /** End of changes - NGMOD **/
    
#endif

}
