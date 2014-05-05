#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vicNl.h"
#include "vicNl_ldas.h"
#include "global.h"
#include "preproc.h"
#include "gribfuncs.h"


/* AUTHOR: Greg O'Donnell [tempgd@hydro.washington.edu]
   DATE:   May 9 2000

   This program is a preprocessor to convert NCEP LDAS files
   into binary forcing files suitable for the revised Princeton VIC code 
   
   Modifications:
   09/14/00	initialise the error message to a null string. Otherwise the 
   		string may contain undefined data and the test of the string
   		contents later in the code can give incorrect results.
   									JS

   06/06/01	re-written to handle variable numbers of records in a grib file.
									JS

   07/07/01	added code to blend precip and solar radiation datasets
									JS

   10/05/01	added code to swap between GOES+EDAS solar radiation and EDAS
		only, based on the date. This is because of bad GOES calibration 
		during the period 01 July 1999 to 23 February 2000 inclusive.
									JS

   10/10/02	- complete re-write of main code.
		- removed code which reads in grib data and put this in a separate 
		function (read_grib). this function searches through a grib file 
		for a given grib parameter id and grib proc id and returns the 
		grib data.
		- now the code just loops through the variable list, extracting 
		data for each variable required instead of reading in every grib
		record and relying on the variables to be in the same order in 
		the grib file
		- added an extra field in the variables table which indicates that
		the data are instantaneous and should be converted to time step 
		averaged values using the function process_ins2ave
		- added the function process_ins2ave
		- replaced code for creating and copying float arrays with
		functions (arr_new, arr_copy)
									JS

   11/25/02	- all EDAS variables are instantaneous except for precipitation
		and so need to be converted by changing the flag in the variable
		attributes table
		- GOES and EDAS radiation need to be combined first and then 
		converted from instantaneous to time averaged. Therefore, store
		the current and old timestep DSWRF and DSWRF1, combine then at the
		end as per usual and then convert them to time averaged values.

*/


/* local funcs */
float *arr_new(int size);
float *arr_copy(float *arr_in, int size);


/* off we go */
int main( int argc, char *argv[])
{

  /* general structs */
  extern global_param_struct global_param;
  filenames_struct           filenames; 

  /* LDAS structs */
  ldas_option_struct	ldas_options;
  ldas_index_struct		ldas_index;

  /* file pointers */
  FILE *gp;                          /* global input file pointer */
  FILE *op[MET_VAR_MAX];             /* output file pointers */
  FILE *qcp;                         /* qc log file pointer */

  dmy_struct          *dmy;          /* day-month-year timeseries struct */

  int i, j, var=0;            /* counters */
  int match;

  Stats stats;

  char fname[BUFSIZ+1];      /* filename for NCEP grib files */
  char fname2[BUFSIZ+1];      /* filename for NCEP grib files */
  char prefix[14+1];

  /* naming convention for output files - names taken from NCEP grib convention */
  /* eg,  wgrib ../2000043013.lsmforce_anal | awk 'BEGIN{FS=":"}{print $4}' */
  /* note - duplicate names have a "1" appended" */ 
  /* fields:
	1: 0 = don't dump to file, 1 = dump to file
	2: variable description
	3: dump file name prefix
	4: grib id number - used to identify the variable in the grib file
	5: variable type, 0 = average/accumulation over the timestep, 1 = instantaneous
	6: min allowable value for QA
	7: max allowable value for QA
  */
  char met_attr[MET_VAR_MAX][MET_ATTR_MAX][BUFSIZ+1] = { 
	{"1", "Edas temp", 			"TMP_",		"84",	"11",	"0",	"213.0",		"333.0"},	/* edas 2m air temp [k] */
	{"1", "Edas spf humidity", 	"SPFH_",	"84",	"51",	"0",	"0.0000001",	"0.03"},	/* edas specific humidity [kg/kg] */
	{"1", "Edas pressure", 		"PRES_",	"84",	"1",	"0",	"5000.0",		"110000.0"},/* edas surface pressure [mb] */
	{"1", "Edas wind u", 		"UGRD_",	"84",	"33",	"0",	"-75.0",		"75.0"},	/* edas wind u [m/s] */
	{"1", "Edas wind v", 		"VGRD_",	"84",	"34",	"0",	"-75.0",		"75.0"},	/* edas wind v [m/s] */
	{"1", "Edas shortwave", 	"DSWRF1_",	"84",	"204",	"0",	"0.0",			"1372.0"},	/* edas downward shortwave [w/m2] */
	{"1", "Edas longwave", 		"DLWRF_",	"84",	"205",	"0",	"75.0",			"750.0"},	/* edas downward longwave [w/m2] */
	{"1", "Edas total precip", 	"APCP1_",	"84",	"61",	"0",	"0.0",			"140.0"},	/* edas total precip [kg/m2] */
	{"0", "Edas conv precip", 	"ACPCP_",	"84",	"63",	"0",	"0.0",			"140.0"},	/* edas convective precip [kg/m2] */
	{"0", "Edas Cape", 			"CAPE_",	"84",	"157",	"0",	"MISSING",		"MISSING"},	/* edas conv. avail. pot. energy */
	{"1", "Goes shortwave", 	"DSWRF_",	"154",	"204",	"0",	"0.0",			"1372.0"},	/* goes downward shortwave [w/m2] */
	{"1", "Goes bright temp", 	"BRTMP_",	"154",	"118",	"0",	"-99999",		"99999"},	/* goes skin temp [K] */
	{"0", "Diffuse radiation", 	"HTSGW_",	"154",	"100",	"0",	"MISSING",		"MISSING"},	/* diffuse radiation */
	{"0", "Par", 				"WVDIR_",	"154",	"101",	"0",	"MISSING",		"MISSING"},	/* par */
	{"1", "Combo total precip", "APCP_",	"155",	"61",	"0",	"0.0",			"140.0"},	/* total precip */
	{"0", "StageIV total precip", 	"APCP2_",	"156",	"61",	"0",	"0.0",			"140.0"}	/* total precip */
  };

  int found[MET_VAR_MAX];
  int dump=TRUE;

  /* variables associated with grib */
  int      nReturn=0;        /* return status from Grib_dec   */
  float *grib_data = (float *) NULL; /* pointer to block of decoded data */
  float *grib_data_old = (float *) NULL;
  float *grib_data_apcp = (float *) NULL;
  float *grib_data_apcp1 = (float *) NULL;
  float *grib_data_dswrf = (float *) NULL;
  float *grib_data_dswrf1 = (float *) NULL;
  float *grib_data_dswrf_old = (float *) NULL;
  float *grib_data_dswrf1_old = (float *) NULL;
  int grib_id, grib_proc;
  
  /* check input argumnets */
  if(argc!=2)
    usage(0);

  fprintf(stderr, "VIC-LDAS pre-processor\n");

  /* open global file using standard VIC open_file */
  initialize_global();
  gp = open_file(argv[1], "r"); 

  /* process "standard" global file options */
  global_param = get_global_param(&filenames, gp);   

  /* process "ldas" global file options */
  get_ldas_global_param(gp, &ldas_options);

  /* make dmy struct - required for file stamps */
  dmy = make_dmy_ldas(&global_param, ldas_options.endhour);

  /* process the masks */
  process_masks(&ldas_options, &ldas_index);

  /* are we writing to file or not? */
  dump = (!strcmp("FALSE", ldas_options.forcing))?FALSE:TRUE;
#if VERBOSE
  if(!dump) fprintf(stderr, "FORCING1 = FALSE so data will NOT be written to file\n");
#endif

  /* create the output files */
  if(dump) {
    for(i=0; i<MET_VAR_MAX; i++){
      if( atoi(met_attr[i][MET_ATTR_USE]) ){
#if INS2AVE
        op[i] = open_file(make_fname(ldas_options.forcing, met_attr[i][MET_ATTR_FNAME], "", &dmy[1]), "w");
#else
        op[i] = open_file(make_fname(ldas_options.forcing, met_attr[i][MET_ATTR_FNAME], "", &dmy[0]), "w");
#endif
     }
    }
  }
  
  /* open the QC log file */
  qcp = open_file(ldas_options.qclog, "w");  

  /* loop over all time periods */
#if INS2AVE
  for(i=1; i<global_param.nrecs; i++){
#else
  for(i=0; i<global_param.nrecs; i++){
#endif

    /* record the number of records */
    stats.nsteps++;

    /* create the forcing filenames */
#ifdef NCEP
    strcpy(fname, make_fname(ldas_options.ncep_met, "", ldas_options.ncep_ext, &dmy[i]));
#else
    sprintf(prefix, "%4.4d/%4.4d%2.2d%2.2d/", dmy[i].year, dmy[i].year, dmy[i].month, dmy[i].day);
    strcpy(fname, make_fname(ldas_options.ncep_met, prefix, ldas_options.ncep_ext, &dmy[i]));
#endif
    fprintf(stderr,"Processing file \"%s\"\n",fname);

    /* initialise */
    for(var=0; var<MET_VAR_MAX; var++)
      found[var] = FALSE;

    /* loop over our database */
    for(var=0; var<MET_VAR_MAX; var++) {

      grib_id = atoi(met_attr[var][MET_ATTR_GRIBID]);
      grib_proc = atoi(met_attr[var][MET_ATTR_GRIBPROC]);

      if(atoi(met_attr[var][MET_ATTR_USE])) {

        stats.ngribs++;

        /* get the data */
        nReturn = read_grib(fname, grib_id, grib_proc, &grib_data, (var==0)?1:0);

        /* the data could not be found */
        if(!nReturn) {

          match=FALSE;
          found[var] = FALSE;
          stats.nmiss++;
          stats.miss[var]+=ldas_index.ncell_comp;

        }

        /* the data were read ok */
        else {

          /* we have a match */
          match = TRUE;
          found[var] = TRUE;
          stats.nmatch++;
          stats.match[var]+=ldas_index.ncell_comp;

#if VERBOSE
          fprintf(stderr, "match: var %d, id = %d proc = %d, desc = \"%s\"\n", var, grib_id, grib_proc, met_attr[var][MET_ATTR_DESC]);
#endif

          /* if it's EDAS or Combo precip then save it */
          if(var == MET_VAR_APCP1) {
            grib_data_apcp1 = arr_copy(grib_data, ldas_index.nr*ldas_index.nc);
          }
          if(var == MET_VAR_APCP) {
            grib_data_apcp = arr_copy(grib_data, ldas_index.nr*ldas_index.nc);
          }

          /* if it's EDAS dswrf then save it */
          if(var == MET_VAR_DSWRF1) {
            grib_data_dswrf1 = arr_copy(grib_data, ldas_index.nr*ldas_index.nc);
          }
           
          /* if it's GOES dswrf then save it */
          if(var == MET_VAR_DSWRF) {
            grib_data_dswrf = arr_copy(grib_data, ldas_index.nr*ldas_index.nc);
          }   

#if INS2AVE
          /* convert from instantaneous to time averaged */
          if(atoi(met_attr[var][MET_ATTR_TYPE]) == 1) {
#if VERBOSE
            fprintf(stderr, "converting %s for record %d using records %d and %d\n", met_attr[var][MET_ATTR_DESC], i, i-1, i);
#endif

#ifdef NCEP
            strcpy(fname2, make_fname(ldas_options.ncep_met, "", ldas_options.ncep_ext, &dmy[i-1]));
#else
            sprintf(prefix, "%4.4d/%4.4d%2.2d%2.2d/", dmy[i-1].year, dmy[i-1].year, dmy[i-1].month, dmy[i-1].day);
            strcpy(fname2, make_fname(ldas_options.ncep_met, prefix, ldas_options.ncep_ext, &dmy[i-1]));
#endif

            /* get the old data from the previous time step */
            if(!read_grib(fname2, grib_id, grib_proc, &grib_data_old, 0)) {
              fprintf(stderr, "Error reading file %s, var %d, id = %d proc = %d, desc = \"%s\"\n", fname, var, grib_id, grib_proc, met_attr[var][MET_ATTR_DESC]);
              exit(1);
            }

            /* store DSWRF and DSWRF1 from previous timestep */
            if(var == MET_VAR_DSWRF1) {
              grib_data_dswrf1_old = arr_copy(grib_data_old, ldas_index.nr*ldas_index.nc);
            }
            if(var == MET_VAR_DSWRF) {
              grib_data_dswrf_old = arr_copy(grib_data_old, ldas_index.nr*ldas_index.nc);
            }   

            /* convert from instantaneous to ave, but not DSWRF as this needs to be combined with DSWRF1 first */
            if(var != MET_VAR_DSWRF)
              process_inst2ave(ldas_index, grib_data, grib_data_old, var, &stats);

          }
#endif

          /* write the data - but not for Combo precip or GOES solar because we may want to blend them with Edas later */
#if 1
          if(var != MET_VAR_APCP && var != MET_VAR_DSWRF) {
#else
          if(var != MET_VAR_APCP) {
#endif
            if(dump) {
#if VERBOSE
              fprintf(stderr, "writing data to %s\n", met_attr[var][MET_ATTR_FNAME]);
#endif
              write_grib(grib_data, ldas_index.cell_num, ldas_index.cell_index, ldas_index.ncell_comp, ldas_options.ncep_out, op[var], var, i, qcp, atof(met_attr[var][MET_ATTR_MINBOUND]), atof(met_attr[var][MET_ATTR_MAXBOUND]), met_attr[var][MET_ATTR_FNAME], &stats);
              stats.written[var]+=ldas_index.ncell_comp;
              stats.nwritten++;
            }
          }

        }

      }

    }

    /* free up array */
    if(grib_data != NULL) { free(grib_data); grib_data = NULL; }  
    if(grib_data_old != NULL) { free(grib_data_old); grib_data_old = NULL; }  

    /* if BRTMP not available for this file then dump field of all missing data */
    var = MET_VAR_BRTMP;
    if(atoi(met_attr[var][MET_ATTR_USE])) {
      if(!found[var]) {
        grib_data = arr_new(ldas_index.nr*ldas_index.nc);
        for(j=0; j<ldas_index.nr*ldas_index.nc; j++)
          grib_data[i] = MISSING;
        found[var] = TRUE;

        if(dump) {
#if VERBOSE
          fprintf(stderr, "writing data to %s\n", met_attr[var][MET_ATTR_FNAME]);
#endif
          write_grib(grib_data, ldas_index.cell_num, ldas_index.cell_index, ldas_index.ncell_comp, ldas_options.ncep_out, op[var], var, i, qcp, atof(met_attr[var][MET_ATTR_MINBOUND]), atof(met_attr[var][MET_ATTR_MAXBOUND]), met_attr[var][MET_ATTR_FNAME], &stats);
          stats.written[var]+=ldas_index.ncell_comp;
          stats.nwritten++;
        }
      }
    }

    if(grib_data != NULL) { free(grib_data); grib_data = NULL; }

    /* process Combo precip */
    var = MET_VAR_APCP;
    if(atoi(met_attr[var][MET_ATTR_USE])) {

      grib_data = arr_new(ldas_index.nr*ldas_index.nc);

      /* use Edas precip if Combo not found */
      if(!found[var] && found[MET_VAR_APCP1]) {
#if VERBOSE
        fprintf(stderr, "*** substituting EDAS precip for Combo precip\n");
#endif
        grib_data = arr_copy(grib_data_apcp1, ldas_index.nr*ldas_index.nc);
        stats.nsubs++;
        stats.subs[var]+=ldas_index.ncell_comp;
        stats.match[var]+=ldas_index.ncell_comp;
        found[var] = TRUE;
      }
      /* if Combo found then blend with Edas */
      else if(found[var] && found[MET_VAR_APCP1]) {
#if VERBOSE
        fprintf(stderr, "*** blending EDAS and Combo precip\n");
#endif
        /* weight the precip at the USA/Canada and USA/Mexico borders */
        process_precip(ldas_index, grib_data_apcp, grib_data_apcp1, var, &stats);
        grib_data = arr_copy(grib_data_apcp, ldas_index.nr*ldas_index.nc);
        stats.nsubs++;
        stats.nblended++;
      }

      if(dump) {
#if VERBOSE
        fprintf(stderr, "writing data to %s\n", met_attr[var][MET_ATTR_FNAME]);
#endif
        write_grib(grib_data, ldas_index.cell_num, ldas_index.cell_index, ldas_index.ncell_comp, ldas_options.ncep_out, op[var], var, i, qcp, atof(met_attr[var][MET_ATTR_MINBOUND]), atof(met_attr[var][MET_ATTR_MAXBOUND]), met_attr[var][MET_ATTR_FNAME], &stats);
        stats.written[var]+=ldas_index.ncell_comp;
        stats.nwritten++;
      }

    }

    /* free up array */
    if (grib_data != NULL) { free(grib_data); grib_data = NULL; }  

#if 1
    /* use EDAS shortwave radiation for missing Goes radiation */
    var = MET_VAR_DSWRF;
    if(atoi(met_attr[var][MET_ATTR_USE]) ) {

      grib_data = arr_new(ldas_index.nr*ldas_index.nc);

      if(!found[var] && found[MET_VAR_DSWRF1]) {
#if VERBOSE
      fprintf(stderr, "*** substituting EDAS dswrf for Goes dswrf\n");
#endif
        grib_data = arr_copy(grib_data_dswrf1, ldas_index.nr*ldas_index.nc);
        found[var] = TRUE;
        stats.nsubs++;
        stats.match[var]+=ldas_index.ncell_comp;
        stats.subs[var]+=ldas_index.ncell_comp;
      }
#ifdef REALTIME
      /* use EDAS between 01 July 1999 to 23 February 2000 inclusive */
      else if((dmy[i].year==1999 && dmy[i].month>=7) || (dmy[i].year==2000 && dmy[i].month==1) ||
               (dmy[i].year==2000 && dmy[i].month==2 && dmy[i].day<=23)) {
printf("SWITCHING!!! %d/%d/%d\n", dmy[i].year, dmy[i].month, dmy[i].day);
        grib_data = arr_copy(grib_data_dswrf1, ldas_index.nr*ldas_index.nc);
        found[var] = TRUE;
        stats.nsubs++;
        stats.match[var]+=ldas_index.ncell_comp;
        stats.subs[var]+=ldas_index.ncell_comp;
      }
#endif
      /* use EDAS if any GOES cells are missing */
      else {
#if VERBOSE
        fprintf(stderr, "*** blending GOES and EDAS dswrf\n");
#endif

        /* infill missing Goes cells with Edas */
        process_dswrf(ldas_index, grib_data_dswrf, grib_data_dswrf1, var, &stats);
        grib_data = arr_copy(grib_data_dswrf, ldas_index.nr*ldas_index.nc);
        stats.nsubs++;
        stats.nblended++;

      }

#if INS2AVE
      grib_data_old = arr_new(ldas_index.nr*ldas_index.nc);

      if(!found[var] && found[MET_VAR_DSWRF1]) {
#if VERBOSE
      fprintf(stderr, "*** substituting EDAS dswrf for Goes dswrf\n");
#endif
        grib_data_old = arr_copy(grib_data_dswrf1_old, ldas_index.nr*ldas_index.nc);
        found[var] = TRUE;
        stats.nsubs++;
        stats.match[var]+=ldas_index.ncell_comp;
        stats.subs[var]+=ldas_index.ncell_comp;
      }
#ifdef REALTIME
      /* use EDAS between 01 July 1999 to 23 February 2000 inclusive */
      else if((dmy[i-1].year==1999 && dmy[i-1].month>=7) || (dmy[i-1].year==2000 && dmy[i-1].month==1) ||
               (dmy[i-1].year==2000 && dmy[i-1].month==2 && dmy[i-1].day<=23)) {
printf("SWITCHING!!! %d/%d/%d\n", dmy[i-1].year, dmy[i-1].month, dmy[i-1].day);
        grib_data_old = arr_copy(grib_data_dswrf1_old, ldas_index.nr*ldas_index.nc);
        found[var] = TRUE;
        stats.nsubs++;
        stats.match[var]+=ldas_index.ncell_comp;
        stats.subs[var]+=ldas_index.ncell_comp;
      }
#endif
      /* use EDAS if any GOES cells are missing */
      else {
#if VERBOSE
        fprintf(stderr, "*** blending GOES and EDAS dswrf\n");
#endif

        /* infill missing Goes cells with Edas */
        process_dswrf(ldas_index, grib_data_dswrf_old, grib_data_dswrf1_old, var, &stats);
        grib_data_old = arr_copy(grib_data_dswrf_old, ldas_index.nr*ldas_index.nc);
        stats.nsubs++;
        stats.nblended++;

      }

      /* convert from instantaneous to time averaged */
      if(atoi(met_attr[var][MET_ATTR_TYPE]) == 1) {
#if VERBOSE
        fprintf(stderr, "converting %s for record %d using records %d and %d\n", met_attr[var][MET_ATTR_DESC], i, i-1, i);
#endif

        /* convert from instantaneous to ave, but not DSWRF as this needs to be combined with DSWRF1 first */
        process_inst2ave(ldas_index, grib_data, grib_data_old, var, &stats);

      }
#endif

      if(dump) {
#if VERBOSE
        fprintf(stderr, "writing data to %s\n", met_attr[var][MET_ATTR_FNAME]);
#endif 
        write_grib(grib_data, ldas_index.cell_num, ldas_index.cell_index, ldas_index.ncell_comp, ldas_options.ncep_out, op[var], var, i, qcp, atof(met_attr[var][MET_ATTR_MINBOUND]), atof(met_attr[var][MET_ATTR_MAXBOUND]), met_attr[var][MET_ATTR_FNAME], &stats);
        stats.written[var]+=ldas_index.ncell_comp;
        stats.nwritten++;
      }

    }
#endif

    /* free some memory */
    if (grib_data != (float *)NULL) { free(grib_data); grib_data = (float *)NULL; }  
    if (grib_data_old != (float *)NULL) { free(grib_data_old); grib_data_old = (float *)NULL; }  
    if (grib_data_apcp != (float *)NULL) { free(grib_data_apcp); grib_data_apcp = NULL; } 
    if (grib_data_apcp1 != (float *)NULL) { free(grib_data_apcp1); grib_data_apcp1 = NULL; } 
    if (grib_data_dswrf != (float *)NULL) { free(grib_data_dswrf); grib_data_dswrf = NULL; } 
    if (grib_data_dswrf1 != (float *)NULL) { free(grib_data_dswrf1); grib_data_dswrf1 = NULL; } 
    if (grib_data_dswrf_old != (float *)NULL) { free(grib_data_dswrf_old); grib_data_dswrf_old = NULL; } 
    if (grib_data_dswrf1_old != (float *)NULL) { free(grib_data_dswrf1_old); grib_data_dswrf1_old = NULL; } 

    /* check that we have output the data that we asked for */
    match = FALSE;
    for(var=0; var<MET_VAR_MAX; var++) {

      /* print out the current variable in the database */
      if(atoi(met_attr[var][MET_ATTR_USE]) && !found[var]) {
        fprintf(stderr, "record NOT found: ");
        fprintf(stderr, "var %d, use = %s, id = %s, desc = \"%s\"\n", var, met_attr[var][MET_ATTR_USE], met_attr[var][MET_ATTR_GRIBID], met_attr[var][MET_ATTR_DESC]);
        match++;
      }

    }

    if(match) {
      fprintf(stderr, "Error: %d requested variables NOT found\n", match);
      return(EXIT_FAILURE);
    }

  } /* end of record loop */

  /* print stats */
  fprintf(stderr, "\n>>> statistics\n");
  fprintf(stderr, "# time steps            = %d\n", stats.nsteps);
  fprintf(stderr, "# grib records read     = %d\n", stats.ngribs);
  fprintf(stderr, "# records matched       = %d\n", stats.nmatch);
  fprintf(stderr, "# records not found     = %d\n", stats.nmiss);
  fprintf(stderr, "# records substituted   = %d\n", stats.nsubs);
  fprintf(stderr, "# records blended       = %d\n", stats.nblended);
  fprintf(stderr, "# records ins2ave       = %d\n", stats.nins2ave);
  fprintf(stderr, "# records written       = %d\n", stats.nwritten);
  fprintf(stderr, "# qc failures           = %ld\n", stats.nqcfails);
  fprintf(stderr, "\n");
  fprintf(stderr, "Variable                   Read  Missed    Subs Blended Ins2Ave Written QCfails\n");
  fprintf(stderr, "-------------------------------------------------------------------------------\n");
  for(var=0; var<MET_VAR_MAX; var++)
    fprintf(stderr, "(%20s) %8d%8d%8d%8d%8d%8d%8ld\n", met_attr[var][MET_ATTR_DESC], stats.match[var], stats.miss[var], stats.subs[var], stats.blended[var], stats.ins2ave[var], stats.written[var], stats.qcfails[var]);

  return(0);

}


/****
*
****/
float *arr_new(int size)
{

  float *arr_out;

  if( (arr_out = (float *)malloc(sizeof(float)*size)) == NULL) {
    fprintf(stderr, "Error: out of memory in arr_new\n"); 
    exit(EXIT_FAILURE);
  }

  return (arr_out);

} /* arr_new */


/****
*
****/
float *arr_copy(float *arr_in, int size)
{

  float *arr_out;

  arr_out = arr_new(size);
  memcpy(arr_out, arr_in, sizeof(float)*size);

  return (arr_out);

} /* arr_copy */


/****
*
****/
void usage(char *msg)
{

  fprintf(stderr, "Usage:\t preproc <control file>\n");
  exit(EXIT_FAILURE);

} /* usage */




