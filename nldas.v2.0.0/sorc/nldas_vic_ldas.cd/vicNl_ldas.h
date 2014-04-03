/*********************************************************************** 
   Author: Greg O'Donnell [tempgd@hydro.washington.edu], June 2000

   Function prototypes and structures for use in the LDAS code 
   
   Modifications:
   	09-20-00	Added EDAS precip (acpc1) to ldas_metfiles_struct.
   													JS
   	09-20-00	Added location (loc) array for grid cell to ldas_index_struct
   			and defined location macros for USA, Canada and Mexico.
   													JS
   	10-10-00	Re-wrote ldas_outfiles_struct to hold file pointers for
   			all LDAS outputs and original output.
   													
***********************************************************************/

/** Start of changes -- JS **/
/* macros for the ldas mask locations */
#define LOCATION_USA 1
#define LOCATION_CANADA 2
#define LOCATION_MEXICO 3
/** End of changes -- JS **/

/*
  Contains file pointers for the NCEP processed met files.
  For deatils of LDAS forcings see GEWEX News, November, 1999.
  Units refer to those in the NCEP file, NOT those required by VIC.
*/
typedef struct {

  FILE *apcp;      /* precipitation, mm's, gauge only, with radar timing */
  FILE *dlwrf;     /* downward longwave, W/m2, GOES satellite based */
  FILE *dswrf;     /* downward shortwave, W/m2, EDAS */
  FILE *pres;      /* surface pressure, Pa, EDAS */
  FILE *spfh;      /* 2-m specific humidity, Kg/Kg, EDAS */
  FILE *tmp;       /* 2-m temperature, K, EDAS */
  FILE *ugrd;      /* 10-m U wind component, m/s, EDAS */
  FILE *vgrd;      /* 10-m V wind component, m/s, EDAS */
/** Start of changes -- JS **/
  FILE *apcp1;      /* precipitation, mm's, EDAS */
/** End of changes -- JS **/

} ldas_metfiles_struct;


/* Pointer for the output filenames. Filenames are based on ALMA standard. Sign conventions are LDAS driven. */

typedef struct {

/** Start of changes -- JS **/
  /* 1) Energy Balance Components */
  FILE *SWnet;			/* NSWRS: Net surface shortwave radiation, W/m2, +ve upward, A-man */
  FILE *LWnet;			/* NLWRS: Net surface longwave radiation, W/m2, +ve upward, A-man */
  FILE *Qle;			/* LHTFL: Latent heat flux, W/m2, +ve upward, A-man */
  FILE *Qh;			/* SHTFL: Sensible heat flux, W/m2, +ve upward, A-man */
  FILE *Qg;			/* GFLUX: Ground heat flux, W/m2, +ve upward, A-man */
  FILE *Qf;			/* SNOHF: Snow phase change heat flux, W/m2, +ve solid to liquid, A-rec */
  FILE *SWdown;			/* DSWRF: Downward surface shortwave radiation, W/m2, +ve downward, L-man */
  FILE *LWdown;			/* DLWRF: Downward surface longwave radiation, W/m2, +ve downward, L-man */
  
  /* 2) Water Balance Components */
  FILE *Snowf;			/* ASNOW: Snowfall (frozen precipitation), kg/m2, +ve always, A-man */		
  FILE *Rainf;			/* ARAIN: Rainfall (unfrozen precipitation), kg/m2, +ve always, A-man */	
  FILE *Evap;			/* EVP  : Total evapotranspiration (all sources), kg/m2, +ve upward, A-man */
  FILE *Qs;			/* SSRUN: Surface runoff, kg/m2, +ve out of cell, A-man */
  FILE *Qsb;			/* BGRUN: Subsurface runoff, kg/m2, +ve out of cell, A-man */
  FILE *Qsm;			/* SNOM : Snowmelt, kg/m2, +ve solid to liquid, A-man */
  FILE *Qgwrec;			/* GWREC: Groundwater recharge, kg/m2, +ve out of cell, L-opt */
  FILE *Qrec;			/* QREC : Flood plain recharge, kg/m2, +ve out of cell, A-opt */
  
  /* 3) Surface State Variables */
  FILE *SnowT;			/* SNOWT: Snow temperature, depth averaged, K, +ve always, A-man */
  FILE *VegT;			/* VEGT : Vegetation canopy temperature, K, +ve always, A-man */
  FILE *BareT;			/* BARET: Bare soil surface temperature, K, +ve always, A-man */
  FILE *AvgSurfT;		/* AVSFT: Average surface temperature, K, +ve always, A-man */
  FILE *RadT;			/* RADT : Effective radiative surface temperature, K, +ve always, A-man */
  FILE *Albedo;			/* ALBDO: Surface albedo, all wavelengths, %, +ve always, A-man */
  FILE *SWE;			/* WEASD: Snow water equivalent, kg/m2, +ve always, A-man */
  FILE *SurfStor;		/* SSTOR: Surface water storage, kg/m2, +ve always, A-man */
  FILE *CanopInt;		/* CWAT : Plant canopy surface water storage, kg/m2, +ve always, L-man */

  /* 4) Subsurface State Variables */
  FILE *SoilTemp[MAX_LAYERS];	/* SOILT: Soil temperature, each layer, K, +ve always, L-man */  
  FILE *SoilMoistTotal;			/* SOILM: Total column soil moisture, kg/m2, +ve always, L-man */
  FILE *SoilMoistRoot;			/* SOILM: Root zone soil moisture, kg/m2, +ve always, A-man */
  FILE *SoilMoist1m;			/* SOILM: Top 1m layer soil moisture, kg/m2, +ve always, L-man */
  FILE *SoilMoist[MAX_LAYERS];	/* SOILM: Average layer soil moisture, kg/m2, +ve always, A-man */
  FILE *LSoilMoist[MAX_LAYERS];	/* LSOIL: Liquid soil moisture, each layer, kg/m2, +ve always, A-opt */
  FILE *SoilWetTotal;			/* MSTAV: Total soil column wetness, %, +ve always, A-man */
  FILE *SoilWetRoot;			/* MSTAV: Root zone wetness, %, +ve always, L-man */

  /* 5) Evaporation Components */
  FILE *ECanop;			/* EVCW : Intercepted evaporation, kg/m2, +ve upwards, L-man */
  FILE *TVeg;			/* TRANS: Vegetation transpiration, kg/m2, +ve upwards, A-man */
  FILE *ESoil;			/* EVBS : Bare soil evaporation, kg/m2, +ve upwards, A-man */
  FILE *EWater;			/* EWATR: Open water evaporation, kg/m2, +ve upwards, A-rec */
  FILE *SubSnow;		/* SBSNO: Snow evaporation (sublimation), kg/m2, +ve upwards, A-rec */
  FILE *EPot;			/* PEVPR: Potential evaporation, kg/m2, +ve upwards, L-rec */  
  FILE *ACond;			/* ACOND: Aerodynamic conductance, m/s, +ve always, A-man */
  FILE *CCond;			/* CCOND: Canopy conductance, m/s, +ve always, L-opt */
  FILE *VGreen;			/* VEG  : Vegetation greenness, %, +ve always, L-opt */
  FILE *LAI;			/* LAI  : Leaf area index, -, +ve always, L-opt */
  FILE *MRough;			/* SFCR : Surface momentum roughness length, m, +ve always, L-opt */
  FILE *HRough;			/* SFCRH: Surface heat roughness length, m, +ve always, L-opt */
  
  /* 6) Streamflow */
  /* none */
  
  /* 7) Cold Season Processes */
  FILE *SnowDepth;		/* SNOD : Snow depth, m, +ve always, A-rec */
  FILE *SnowFrac;		/* SNOC : Snow cover, %, +ve always, A-opt */
  FILE *SAlbedo;		/* SALBD: Snow albedo, %, +ve always, A-rec */
  
  /* Non-LDAS outputs, but are useful anyway */
  FILE *Precip;			/* Total precipitation, kg/m2, +ve always, NON-LDAS */
  FILE *RNet;			/* Net radiation, W/m2, +ve upwards???, NON-LDAS  */
  FILE *Air_Temp;		/* Air temperature, C, -, NON-LDAS*/
  FILE *Wind;			/* Wind speed, m/s, +ve always, NON-LDAS */
  FILE *RHumid;			/* Relative Humidity, %, +ve always, NON-LDAS */

  FILE *InLong;			/* Incoming longwave W/m2 */
  FILE *NetLong;		/* Net longwave W/m2 */

  /* others */
  FILE *statefile;    /* State file */
/** End of changes -- JS **/

} ldas_outfiles_struct;

/* LDAS global options */
typedef struct {
  char mask[BUFSIZ+1];           /* mask of domain with met cells identified */
  char pp_output[BUFSIZ+1];      /* postprocessed output */
  char processor_mask[BUFSIZ+1]; /* processor mask */
/** Start of changes -- JS **/
  char pweight_mask[BUFSIZ+1]; 		 /* precip weights */
/** End of changes -- JS **/
  char lats_table[BUFSIZ+1];     /* lats table */
  int  processor;                /* processor mask id to run */
  int  endhour;                  /* end hour of simulation */
  int statehour;                 /* hour to write out state file */

  /* preprocessor stuff */
  char ncep_met[BUFSIZ+1];	/* path of input files */
  char forcing [BUFSIZ+1];	/* path of output files */
/*  char mask    [BUFSIZ+1];*/	/* mask path/filename */
  char ncep_ext[BUFSIZ+1];	/* extension for NCEP filenames */
  char qclog     [BUFSIZ+1];	/* QC log filename */
  char ncep_out;                /* output format for NCEP met 1=binary */
                                /* 0=ascii */
} ldas_option_struct;

/* mask data information */
/* NOTE, ldas_mask.cells contains the cell # to  calculate
   ie from 1:N - whilst array C indexing is 0:N-1 */
typedef struct {
  float xll;      	/* longitude of bottom right corner - cell centre */
  float yll;      	/* latitude of bottom right corner - cell centre */
  float size;     	/* cell resolution */
  int   nc;       	/* number of columns */
  int   nr;       	/* number of rows */
  
  int ncell_comp; 	/* number of cells to be computed */
  int ncell_met;  	/* number of cells to be computed */
  int *cell_num;  	/* contains the cell #s as identified in the soil file [1:N] */
  int *cell_index;	/* elements in met data file [0:N-1] corresponding to cell_num */
/** Start of changes -- JS **/
  int *cell_loc;  /* contains the location (USA=1, Canada=2 or Mexico=3) for each cell num */
  int *cell_pweight;/* contains the precip blending weights for each cell num */
/** End of changes -- JS **/
} ldas_index_struct;

/* function prototypes */
/* calls open_ldas_metfiles and open_ldas_outfiles */
void open_ldas_files    (ldas_metfiles_struct *ldas_metfiles, 
			 ldas_outfiles_struct *ldas_outfiles,
			 filenames_struct     *fname,
			 dmy_struct           *dmy );
void open_ldas_metfiles( ldas_metfiles_struct *ldas_metfiles, 
			 filenames_struct     *fname,
			 dmy_struct           *dmy,
			 char                 *mode);
void open_ldas_outfiles( ldas_outfiles_struct *ldas_outfiles, 
			 filenames_struct     *fname,
			 dmy_struct           *dmy,
			 char                 *mode);


/* get ldas specific options from global file */
void get_ldas_global_param(FILE *gp, ldas_option_struct *ldas_options );
/* Create filenames, which include a datestamp */
char * make_fname(char *path, char *prefix, char *suffix, dmy_struct *dmy );
/* process the ldas mask files */
void process_masks(ldas_option_struct *options, ldas_index_struct *ldas_index);

/* read in and initialize met data */
void initialize_ldas_atmos(atmos_data_struct    *atmos,
			   ldas_metfiles_struct *ldas_metfiles,
			   double               *Tfactor,
			   int                  index,
			   int                  ncells_met,
			   int					cell_loc);
/* slight modifications to allow subdaily runs */
dmy_struct *make_dmy_ldas(global_param_struct *global, int endhour);
/* minor changes */
void dist_prec_ldas(atmos_data_struct   *atmos,
		    dist_prcp_struct    *prcp,
		    soil_con_struct     *soil_con,
		    veg_con_struct      *veg_con,
		    dmy_struct          *dmy,
		    global_param_struct *global_param,
		    ldas_outfiles_struct     *outfiles,
		    int                  rec,
		    int                  cellnum,
		    char                 NEWCELL,
		    char                 LASTREC,
		    int                  statehour,
		    int                  dry_time,
/** Start of changes -- JS **/
		    int                 still_storm);
/** End of changes -- JS **/
/* minor changes */
void put_data_ldas(dist_prcp_struct  *prcp,
		atmos_data_struct 	*atmos,
		veg_con_struct    	*veg_con,
		ldas_outfiles_struct   *outfiles,
		double            	*depth,
		double            	*dz,
		double             	dp,
		double            	*AreaFract,
		double            	*porosity,
		dmy_struct       	*dmy,
		int                	rec,
		int                	dt,
		int                	Nnodes,
		int					skipyear,
/** Start of changes -- JS **/
		double               MAX_SNOW_TEMP,
		double               MIN_RAIN_TEMP);
/** End of changes -- JS **/
  /* new write data routine */
void write_data_ldas(out_data_struct *out_data,
		     ldas_outfiles_struct *outfiles,
		     dmy_struct      *dmy,
		     int              dt);

void write_model_state_ldas(dist_prcp_struct *, global_param_struct *, int, 
			    int, ldas_outfiles_struct *, soil_con_struct *,
			    veg_con_struct *veg_con, int DRY_TIME
/** Start of changes -- JS **/
		        ,int still_storm);
/** End of changes -- JS **/

void read_initial_model_state_ldas(FILE *, dist_prcp_struct *, 
				global_param_struct *, int, int, int, 
				soil_con_struct *, veg_con_struct *veg_con, int *dry_time
/** Start of changes -- JS **/
		        ,int *still_storm);
/** End of changes -- JS **/

void initialize_model_state_ldas(dist_prcp_struct    *prcp,
				 dmy_struct           dmy,
				 double               surf_temp,
				 global_param_struct *global_param,
				 infiles_struct       infiles, 
				 int                  cellnum,
				 int                  Nveg,
				 int                  Nnodes,
				 int                  Ndist,
				 soil_con_struct     *soil_con,
				 veg_con_struct      *veg_con,
				 int                 *dry_time
/** Start of changes -- JS **/
		         ,int               *still_storm);
/** End of changes -- JS **/
