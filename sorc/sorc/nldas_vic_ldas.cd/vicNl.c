#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>

/** Start of changes - NGMOD **/
#include <vicNl_ldas.h>
/** End of changes - NGMOD **/

static char vcid[] = "$Id: vicNl.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

/** Main Program **/

int main(int argc, char *argv[])
/**********************************************************************
	vicNl.c		Dag Lohmann		January 1996

  This program controls file I/O and variable initialization as well as
  being the primary driver for the model.

  For details about variables, input files and subroutines check:
	http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC_home.html

  UNITS: unless otherwise marked:
         all water balance components are in mm
	 all energy balance components are in mks
	 depths, and lengths are in m

  modifications:
  1997-98 Model was updated from simple 2 layer water balance to 
          an extension of the full energy and water balance 3 layer
	  model.                                                  KAC

  May, 2000
          This is a revised code for the LDAS project. 
	                                                         GMOD
**********************************************************************/
{

  extern veg_lib_struct *veg_lib;
  extern option_struct options;

#if LINK_DEBUG
  extern debug_struct debug;
#endif
  extern Error_struct Error;
  extern global_param_struct global_param;

  /** Variable Declarations **/

  char                     NEWCELL;
  char                     LASTREC;
  char                     MODEL_DONE;
  int                      rec, i, j;
  int                      veg;
  int                      band;
  int                      Ndist;
  int                      Nveg_type;
  int                      cellnum;
  int                      index;
  int                      RUN_MODEL;
  int                      Ncells;
  int                      cell_cnt;
  double                   storage;
  double                   veg_fract;
  double                   band_fract;
  dmy_struct              *dmy;
  atmos_data_struct       *atmos;
  veg_con_struct          *veg_con;
  soil_con_struct          soil_con;
  dist_prcp_struct         prcp; /* stores information about distributed 
				    precipitation */

  filenames_struct         filenames;
  infiles_struct           infiles;
  /** Start of changes - NGMOD **/
  /* obsolete structures */
  /* outfiles_struct          outfiles; */
  /* filenames_struct         builtnames; */
  /** End of changes - NGMOD **/

  /** Start of changes - NGMOD **/
  /** Additional variables required for VIC LDAS **/
  ldas_metfiles_struct ldas_metfiles;
  ldas_outfiles_struct ldas_outfiles;
  ldas_option_struct   ldas_options;
  ldas_index_struct    ldas_index;
  int                  dry_time;
/** Start of changes -- JS **/
  int                 still_storm;
/** End of changes -- JS **/
  /** End of changes - NGMOD **/

#if VERBOSE
  fprintf(stderr,"Running Model Version: %s\n",vcid);
#endif

  /** Read Model Options **/
  initialize_global();
  filenames = cmd_proc(argc, argv);

  /** Read Global Control File **/
  infiles.globalparam = open_file(filenames.global,"r");
  global_param = get_global_param(&filenames, infiles.globalparam);
  /** Start of changes - NGMOD **/
  /** get the real time LDAS parameters seperately and process */
   get_ldas_global_param(infiles.globalparam, &ldas_options);
   process_masks(&ldas_options, &ldas_index);
  /** End of changes - NGMOD **/

  /** Check and Open Files **/
  check_files(&infiles, &filenames);

  /** Check and Open Debugging Files **/
#if LINK_DEBUG
  open_debug();
#endif

  /** Read Vegetation Library File **/
  veg_lib = read_veglib(infiles.veglib,&Nveg_type);

  /** Initialize Parameters **/
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  cellnum = -1;

  /** Make Date Data Structure **/
  /** Start of changes - NGMOD **/
  /*  slight code change to allow part days to be specified */
  dmy      = make_dmy_ldas(&global_param, ldas_options.endhour);
  /** End of changes - NGMOD   **/

  /** allocate memory for the atmos_data_struct **/
  alloc_atmos(global_param.nrecs, &atmos);

  /** Start of changes - NGMOD **/
  /** open the ldas metfiles and output files, note, this must be done after
      initializing dmy - date used in filenames */
  open_ldas_files( &ldas_metfiles, &ldas_outfiles, &filenames, &dmy[0] );
  /** End of changes - NGMOD   **/

#if SAVE_STATE
  /** open state file if model state is to be saved **/
  if ( global_param.stateyear > 0 ) 
    ldas_outfiles.statefile = open_state_file(&global_param, options.Nlayer, 
					 options.Nnode);
#endif


  /** Start of changes - NGMOD **/
  /* obsolete
  if ( options.INIT_STATE ) 
    infiles.statefile = check_state_file(filenames.init_state, dmy[0], 
					 &global_param, options.Nlayer, 
					 options.Nnode);
  */
  if ( options.INIT_STATE ) 
    infiles.statefile = open_file(filenames.init_state,"r");
  /** End of changes - NGMOD **/

  /************************************
    Run Model for all Active Grid Cells
    ************************************/
  MODEL_DONE = FALSE;
  cell_cnt=0;
  while(!MODEL_DONE) {
    if(!options.ARC_SOIL) {

      if((fscanf(infiles.soilparam, "%d", &flag))!=EOF) {
	if(flag) RUN_MODEL=TRUE;
	else     RUN_MODEL=FALSE;
      }
      else {
	MODEL_DONE = TRUE;
	RUN_MODEL = FALSE;
      }
      /** Start of changes - NGMOD **/
      /* ldas_index.cell_num[] contains a list of gridcel-s to run */
      /* Leave code above ie RUN_MODEL..., but clobber result here. */
      /* This avoids changes to soil file or read_soil_param */
      if(ldas_index.ncell_comp==cell_cnt){
	MODEL_DONE = TRUE;
	RUN_MODEL  = FALSE;
      }
      if(!MODEL_DONE) {
	soil_con = read_soilparam(infiles.soilparam);
	if(soil_con.gridcel == ldas_index.cell_num[cell_cnt]){ 
	  RUN_MODEL=TRUE;
	  cell_cnt++;
	}
	else if (soil_con.gridcel > ldas_index.cell_num[cell_cnt]){
	  fprintf(stderr,"Soilfile gridcel exceeds ldas_index.cell_num[] at:\t%d\n",
		  ldas_index.cell_num[cell_cnt]);
	  fprintf(stderr,"Probable miss-ordering of soilfile.\n");
	  exit(EXIT_FAILURE);
	  RUN_MODEL=FALSE;
	}
	else
	  RUN_MODEL=FALSE;
      }
      /** End of changes - NGMOD **/
    }
    else {
      soil_con = read_soilparam_arc(infiles.soilparam, 
				    filenames.soil_dir, &Ncells, 
				    &RUN_MODEL, cell_cnt);
      cell_cnt++;
      if(cell_cnt==Ncells) MODEL_DONE = TRUE;
    }
    if(RUN_MODEL) {

#if LINK_DEBUG
      if(debug.PRT_SOIL) write_soilparam(&soil_con); 
#endif

#if QUICK_FS
      /** Allocate Unfrozen Water Content Table **/
      if(options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  soil_con.ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
	}
	for(i=0;i<MAX_NODES;i++) {
	  soil_con.ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
	}
      }
#endif

      NEWCELL=TRUE;
      cellnum++;

      /** Read Grid Cell Vegetation Parameters **/
      veg_con = read_vegparam(infiles.vegparam, soil_con.gridcel,
                              Nveg_type);
      calc_root_fractions(veg_con, &soil_con);
#if LINK_DEBUG
      if(debug.PRT_VEGE) write_vegparam(veg_con); 
#endif

      /** Build Gridded Filenames, and Open **/
      /** Start of changes - NGMOD **/
      /** Comment: not required, met data no longer held on a cell basis **/
      /*
      builtnames = make_in_and_outfiles(&infiles, &filenames, &soil_con,
					&outfiles);  
      */
      /** End of changes - NGMOD **/
      
      /** Read Elevation Band Data if Used **/
      read_snowband(infiles.snowband,soil_con.gridcel,
		    (double)soil_con.elevation,
		    &soil_con.Tfactor,&soil_con.Pfactor,&soil_con.AreaFract);

      /** Make Precipitation Distribution Control Structure **/
      prcp     = make_dist_prcp(veg_con[0].vegetat_type_num, 
				&options.Nnode);

      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/

#if VERBOSE
      fprintf(stderr,"Initializing Forcing Data\n");
#endif

      /** Start of changes - NGMOD **/
      /* remove the following, new met data initilaisation routines included */
      /* initialize_atmos(atmos, dmy, infiles.forcing,
		       (double)soil_con.time_zone_lng, (double)soil_con.lng,
		       (double)soil_con.lat, soil_con.elevation,
		       soil_con.annual_prec, global_param.wind_h, 
		       soil_con.rough, soil_con.Tfactor); */

      /*      fread(&tmp,sizeof(float),1,ldas_metfiles.tmp); */
      /*fprintf(stderr,"%f\n",tmp); */
      /*exit(0); */

      initialize_ldas_atmos(atmos,
			    &ldas_metfiles,
			    soil_con.Tfactor,
/** Start of changes -- JS **/
/* read in met data for only the active cells as the preprocessor only 
   processes met data for the active cells and not for the whole mask */
                            /*ldas_index.cell_index[cell_cnt-1],*/
			    /*ldas_index.ncell_met,*/ 
                            cell_cnt-1,
			    ldas_index.ncell_comp, 
/** End of changes - JS **/
			    ldas_index.cell_pweight[cell_cnt-1]);

      /** End of changes - NGMOD **/
      
#if LINK_DEBUG
      if(debug.PRT_ATMOS) write_atmosdata(atmos, global_param.nrecs);
#endif

      /**************************************************
        Initialize Energy Balance and Snow Variables 
      **************************************************/

#if VERBOSE
      fprintf(stderr,"Model State Initialization\n");
#endif
      initialize_model_state_ldas(&prcp, dmy[0], atmos[0].air_temp[NR], 
			     &global_param, infiles, soil_con.gridcel, 
			     veg_con[0].vegetat_type_num, 
			     options.Nnode, Ndist, &soil_con, veg_con
			     /** Start of changes - NGMOD **/
			     /* add new variable to call */
			     , &dry_time
/** Start of changes -- JS **/
			     , &still_storm
/** End of changes -- JS **/
			     /** End of changes - NGMOD **/
			     );


#if VERBOSE
      fprintf(stderr,"Running Model\n");
#endif

      /** Update Error Handling Structure **/
      /** Start of changes - NGMOD **/
      /* Error.outfp = outfiles; */
      /** End of changes - NGMOD **/
      Error.infp = infiles;

      /***************************************************
	Intialize Moisture and Energy Balance Error Checks
	***************************************************/
      storage = 0.;
      for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ ) {
	if ( veg < veg_con[0].vegetat_type_num ) veg_fract = veg_con[veg].Cv;
	else veg_fract = ( 1.0 - veg_con[0].Cv_sum );
	for ( band = 0; band < options.SNOW_BAND; band++ ) {
	  band_fract = soil_con.AreaFract[band];

	  if ( veg_fract > SMALL && band_fract > SMALL ) {
	    for(index=0;index<options.Nlayer;index++)
	      /** Start of changes - NGMOD **/
	      /*  changes needed here for mass balance with INIT_STATE */

	      if(options.INIT_STATE)
		storage +=
		  ( prcp.mu[veg]*prcp.cell[WET][veg][band].layer[index].moist +
		   (1-prcp.mu[veg])*prcp.cell[DRY][veg][band].layer[index].moist)
		  * veg_fract * band_fract;
	    /** End of changes - NGMOD **/
	      else
		storage += prcp.cell[WET][veg][band].layer[index].moist 
		  * veg_fract * band_fract;

	    storage += prcp.snow[veg][band].swq * 1000. * veg_fract 
	      * band_fract;

	    if ( veg != veg_con[0].vegetat_type_num ) {
	      /** Start of changes - NGMOD **/
	      /*  changes needed here for mass balance with INIT_STATE */

	      if(options.INIT_STATE)
		storage += (prcp.mu[veg]*prcp.veg_var[WET][veg][band].Wdew+
			    (1-prcp.mu[veg])*prcp.veg_var[DRY][veg][band].Wdew)
		  * veg_fract * band_fract;
	      /** End of changes - NGMOD **/
	      else
		storage += prcp.veg_var[WET][veg][band].Wdew 
		  * veg_fract * band_fract;
	      storage += prcp.snow[veg][band].snow_canopy * 1000. 
		* veg_fract * band_fract;

	    }
	  }
	}
      }
      calc_water_balance_error(-global_param.nrecs,0.,0.,storage);
      calc_energy_balance_error(-global_param.nrecs,0.,0.,0.,0.,0.);

      /******************************************
	Run Model in Grid Cell for all Time Steps
	******************************************/

      for ( rec = 0 ; rec < global_param.nrecs; rec++ ) {
       if ( rec == global_param.nrecs - 1 ) LASTREC = TRUE;
        else LASTREC = FALSE;
        
	/** Start of changes - NGMOD **/
	/* dist_prec(&atmos[rec], &prcp, &soil_con, veg_con,
	   dmy, &global_param, &outfiles, rec, cellnum,
	   NEWCELL, LASTREC); */
	dist_prec_ldas(&atmos[rec], &prcp, &soil_con, veg_con,
		       dmy, &global_param, &ldas_outfiles, rec, cellnum,
		       NEWCELL, LASTREC, ldas_options.statehour, dry_time, still_storm);
	/** End of changes - NGMOD **/
        NEWCELL=FALSE;
	
      }	/* End Rec Loop */

      /* Start of changes - NGMOD **/
      /* close_files(&infiles,&outfiles,&builtnames); */
      /** End of changes - NGMOD **/

#if QUICK_FS
      if(options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_layer[i][j]);
	  free((char *)soil_con.ufwc_table_layer[i]);
	}
	for(i=0;i<MAX_NODES;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_node[i][j]);
	  free((char *)soil_con.ufwc_table_node[i]);
	}
      }
#endif
      free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
      free_vegcon(&veg_con);
      free((char *)soil_con.AreaFract);
      free((char *)soil_con.Tfactor);
      free((char *)soil_con.Pfactor);
      for(index=0;index<=options.Nlayer;index++) 
	free((char*)soil_con.layer_node_fract[index]);
      free((char*)soil_con.layer_node_fract);
    }	/* End Run Model Condition */
    
  } 	/* End Grid Loop */

  /** cleanup **/
  free_atmos(global_param.nrecs, &atmos);

  return EXIT_SUCCESS;
}	/* End Main Program */
