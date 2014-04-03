#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vicNl.h"
#include "vicNl_ldas.h"

static char vcid[] = "$";

#define MET_ATTR_USE		0
#define MET_ATTR_DESC		1
#define MET_ATTR_FNAME		2
#define MET_ATTR_GRIB		3
#define MET_ATTR_MINBOUND	4
#define MET_ATTR_MAXBOUND	5
#define MET_ATTR_MAX		6

enum MET_VARS {
	MET_VAR_TMP=0,
	MET_VAR_SPFH,
	MET_VAR_PRES,
	MET_VAR_UGRD,
	MET_VAR_VGRD,
	MET_VAR_DSWRF1,
	MET_VAR_DLWRF,
	MET_VAR_APCP1,
	MET_VAR_ACPCP,
	MET_VAR_CAPE,
	MET_VAR_DSWRF,
	MET_VAR_BRTMP,
	MET_VAR_HTSGW,
	MET_VAR_WVDIR,
	MET_VAR_APCP,
	MET_VAR_APCP2,
        MET_VAR_MAX        
};

/* start of changes -- JS */
char met_attr[MET_VAR_MAX][MET_ATTR_MAX][BUFSIZ+1] = { 
	{"1", "Edas temp", 		"TMP_",		 "11",	"213.0",	"333.0"},	/* edas 2m air temp [k] */
	{"1", "Edas spf humidity", 	"SPFH_",	 "51",	"0.0000001",	"0.03"},	/* edas specific humidity [kg/kg] */
	{"1", "Edas pressure", 		"PRES_",	  "1",	"5000.0",	"110000.0"},	/* edas surface pressure [mb] */
	{"1", "Edas wind u", 		"UGRD_",	 "33",	"-75.0",	"75.0"},	/* edas wind u [m/s] */
	{"1", "Edas wind v", 		"VGRD_",	 "34",	"-75.0",	"75.0"},	/* edas wind v [m/s] */
	{"1", "Edas shortwave", 	"DSWRF1_",	"204",	"0.0",		"1372.0"},	/* edas downward shortwave [w/m2] */
/*	{"1", "Edas longwave", 		"DLWRF_",	"205",	"75.0",		"750.0"},*/	/* edas downward longwave [w/m2] */
	{"1", "Edas longwave", 		"DLWRF_",	"205",	"50.0",		"750.0"},	/* edas downward longwave [w/m2] */
	{"1", "Edas total precip", 	"APCP1_",	 "61",	"0.0",		"140.0"},	/* edas total precip [kg/m2] */
	{"1", "Edas conv precip", 	"ACPCP_",	 "63",	"0.0",		"140.0"},	/* convective precip [kg/m2] */
	{"0", "Cape", 			"CAPE_",	"157",	"MISSING",	"MISSING"},	/* convective available potential energy */
	{"1", "Goes shortwave", 	"DSWRF_",	"204",	"0.0",		"1372.0"},	/* goes (pinker) downward shortwave [w/m2] */
	{"0", "Goes bright temp", 	"BRTMP_",	"118",	"MISSING",	"MISSING"},	/* goes skin temp [K] */
	{"0", "Diffuse radiation", 	"HTSGW_",	"100",	"MISSING",	"MISSING"},	/* diffuse radiation */
	{"0", "Par", 			"WVDIR_",	"101",	"MISSING",	"MISSING"},	/* par */
	{"1", "Combo total precip", 	"APCP_",	 "61",	"0.0",		"140.0"},	/* total precip */
	{"0", "StageIV total precip", 	"APCP2_",	 "61",	"0.0",		"140.0"}	/* total precip */
  };

void qc_check(int, double);

/* end of changes -- JS */


void initialize_ldas_atmos(atmos_data_struct    *atmos,
			   ldas_metfiles_struct *ldas_metfiles,
			   double               *Tfactor,
			   int                  index,
			   int                  ncells_met,
			   int			cell_pweight)
/**********************************************************************
  initialize_ldas_atmos	Greg O'Donnell		June, 2000

  This routine is loosely based on initialize_atmoc.c.

  The mtclim and other parameterizartions have been removed - this
  decision has been made due to the anticipated difficulties in 
  running such code over short time periods.

  At a later stage error checking etc will be implemented into this 
  routine.

  Modifications:
  	092000	The cell location (USA, Canada, Mexico) passed as an
  		input argument and used to determine which precip
  		data to use: station data (cell located in USA)
  		or EDAS model prediction (cell not in USA or station
  		data are missing).
  								JS
  	060201	Added QC checks on precip					
								JS
  	030501	Precip weighting map replaces location (USA..) map.
		Precip values are now blended between observed station
		data and EDAS predicted precip.
  								JS
	25 Jun 01
		Precip blending and other changes removed and are now 
		done in pre-processor. QC checks done on each forcing 
		variable with new function qc_check().

**********************************************************************/
{
  extern option_struct       options;
  extern param_set_struct    param_set;
  extern global_param_struct global_param;
  extern int                 NR, NF;

  int     rec;
  int     band;

  long    seek;

  float   tmp_met, tmp_met1;  /* temporary storage for met data */
  float   u_wind;
  float   v_wind;

  double  min_Tfactor;
  double  svp_tair;

#if VERBOSE 
  fprintf(stderr,"\nRead meteorological forcing file\n");
#endif

   /************************************************* 
     Create sub-daily precipitation if not provided 
   *************************************************/ 

  /* place here for now - should really be elsewhere */
  if(param_set.FORCE_DT[param_set.TYPE[PREC].SUPPLIED-1] != 1) {
    fprintf(stderr,"initialize_ldas_atmos requires input at the ");
    fprintf(stderr,"same temporal resolution as the model run.\n");
    exit(EXIT_FAILURE);
  }
  
  for(rec = 0; rec < global_param.nrecs; rec++) {

    seek = (index+rec*ncells_met)*sizeof(float);

    /* only check fseek status on precip - if this is OK, assume all are */
    if(fseek(ldas_metfiles->apcp,seek,SEEK_SET)){
      fprintf(stderr,"fseek failure in intialize_atmos_data: %ld.\n", seek);
      exit(EXIT_FAILURE);
    }
    fseek(ldas_metfiles->dlwrf,seek,SEEK_SET);
    fseek(ldas_metfiles->dswrf,seek,SEEK_SET);
    fseek(ldas_metfiles->pres, seek,SEEK_SET);
    fseek(ldas_metfiles->tmp,  seek,SEEK_SET);
    fseek(ldas_metfiles->ugrd, seek,SEEK_SET);
    fseek(ldas_metfiles->vgrd, seek,SEEK_SET);
    fseek(ldas_metfiles->spfh, seek,SEEK_SET);
#if 0
/** Start of changes -- JS **/
    /* seek the EDAS model predicted precip */
    fseek(ldas_metfiles->apcp1, seek,SEEK_SET);
/** End of changes -- JS **/
#endif

    /* PROCESS MET DATA */
    /* Note, the input files are in float and the atmos structure is double */

    /* precipitation - no unit conversion required */
    fread(&tmp_met,sizeof(float),1,ldas_metfiles->apcp);
    qc_check(MET_VAR_APCP, tmp_met);
    atmos[rec].prec[0] = (double)tmp_met;

    /* downward longwave - no unit conversion required */
    fread(&tmp_met,sizeof(float),1,ldas_metfiles->dlwrf);
    qc_check(MET_VAR_DLWRF, tmp_met);
    atmos[rec].longwave[0] = (double)tmp_met;

    /* downward shortwave - no unit conversion required */
    fread(&tmp_met,sizeof(float),1,ldas_metfiles->dswrf);
    qc_check(MET_VAR_DSWRF, tmp_met);
    atmos[rec].shortwave[0] = (double)tmp_met;

    /* pressure - convert from Pa to KPa */
    fread(&tmp_met,sizeof(float),1,ldas_metfiles->pres);
    qc_check(MET_VAR_PRES, tmp_met);
    atmos[rec].pressure[0] = (double)tmp_met/1000.0;

    /* temperature - convert from K to C */
    fread(&tmp_met,sizeof(float),1,ldas_metfiles->tmp);
    qc_check(MET_VAR_TMP, tmp_met);
    atmos[rec].air_temp[0] = (double)tmp_met-KELVIN;

    /* wind - no unit conversion required */
    /* 2-wind components */
    fread(&v_wind,sizeof(float),1,ldas_metfiles->vgrd);
    qc_check(MET_VAR_VGRD, v_wind);
    fread(&u_wind,sizeof(float),1,ldas_metfiles->ugrd);
    qc_check(MET_VAR_UGRD, u_wind);
    atmos[rec].wind[0] = sqrt((double)(u_wind*u_wind + v_wind*v_wind));

    /* specific humidity related calcs */
    fread(&tmp_met,sizeof(float),1,ldas_metfiles->spfh);
    qc_check(MET_VAR_SPFH, tmp_met);
    /* vp */
    atmos[rec].vp[0]  = tmp_met * atmos[rec].pressure[0] / 0.622;  
    /* vpd */
    atmos[rec].vpd[0] = svp(atmos[rec].air_temp[0]) - atmos[rec].vp[0];
    /* density */
    atmos[rec].density[0] = 3.486*atmos[rec].pressure[0]/(275.0+atmos[rec].air_temp[0]); 

/* printf("%f %f %f %f %f %f %f\n", atmos[rec].prec[0] , atmos[rec].longwave[0], atmos[rec].shortwave[0] , atmos[rec].pressure[0] , atmos[rec].air_temp[0] , atmos[rec].wind[0] , atmos[rec].vpd[0] );
*/

  }

  /* determine if it is a snow step */
  /* due to Tfactor overhead do outside rec loop */

  /* determine the largest temp offset */
   min_Tfactor = Tfactor[0];
   for (band = 1; band < options.SNOW_BAND; band++) {
     if (Tfactor[band] < min_Tfactor)
       min_Tfactor = Tfactor[band];
   }

   for (rec = 0; rec < global_param.nrecs; rec++) {
     if ((atmos[rec].air_temp[0] + min_Tfactor) < global_param.MAX_SNOW_TEMP
	 &&  atmos[rec].prec[0] > 0) 
       atmos[rec].snowflag[0] = TRUE;
     else
       atmos[rec].snowflag[0] = FALSE;
   }   
 
}

/****
* qc_check: 
****/
void qc_check(int var, double value)
{

  double minbound;
  double maxbound;
  
  minbound = atof(met_attr[var][MET_ATTR_MINBOUND]);
  maxbound = atof(met_attr[var][MET_ATTR_MAXBOUND]);
  if(value < minbound || value > maxbound) {
    
    fprintf(stderr, "ERROR: %s %f out of range (%f - %f) in intialize_atmos_data (qc_check)\n", met_attr[var][MET_ATTR_DESC], value, minbound, maxbound);
    exit(EXIT_FAILURE);

  }

}

