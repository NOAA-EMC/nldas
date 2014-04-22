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
/***** Version Information *****/
#define VERSION         "VIC Release 4.0.4"
#include "misc.h"
#if (defined SPMD)
#include "mpi.h"
#endif
#include "user_def.h"
#include "snow.h"

/***** Model Constants *****/
#define MAXSTRING    512
#define MINSTRING    20
#define HUGE_RESIST  1.e20      /* largest allowable float number */
#define SMALL        1.e-12     /* smallest allowable float number */
#define MISSING      ( -9999. ) /* missing value for multipliers in
                                   BINARY format */
#define LITTLE 1              /* little-endian flag */
#define BIG 2                 /* big-endian flag */
#define BROOKS 1              /* Brooks-Corey parameters for unsaturated flow */
#define ASCII 1     /* met file format flag */
#define BINARY 2    /* met file format flag */
#define MAX_SNOW_TEMP 1.5 
#define MIN_RAIN_TEMP -0.5 

/***** Forcing Variable Types *****/
#define N_FORCING_TYPES 13
#define AIR_TEMP  0  /* air temperature per time step (C) */
#define ALBEDO    1  /* surface albedo (fraction) */
#define DENSITY   2  /* atmospheric density (kg/m^3) */
#define LONGWAVE  3  /* incoming longwave radiation (W/m^2) */
#define PREC      4  /* precipitation (mm) */
#define PRESSURE  5  /* atmospheric pressure (kPa) */
#define SHORTWAVE 6  /* incoming shortwave (W/m^2) */
#define TMAX      7  /* maximum daily temperature (C) */
#define TMIN      8  /* minimum daily temperature (C) */
#define TSKC      9  /* cloud cover (fraction) */
#define VP        10 /* vapor pressure (kPa) */
#define WIND      11 /* wind speed (m/s) */
#define SKIP      12 /* place holder for unused data columns */

/***** Physical Constants *****/
#define BARE_SOIL_ALBEDO 0.2        /* albedo for bare soil */
#define RESID_MOIST      0.0        /* define residual moisture content 
                                       of soil column */
#define ice_density      917.       /* density of ice (kg/m^3) */
#define T_lapse          6.5        /* tempreature lapse rate of US Std 
                                       Atmos in C/km */
#define von_K        0.40      /* Von Karmin constant for evapotranspiration */
#define KELVIN       273.15    /* conversion factor C to K */
#define STEFAN_B     5.6696e-8 /* stefan-boltzmann const in unit W/m^2/K^4 */
#define Lf           3.337e5   /* Latent heat of freezing (J/kg) at 0C */
#define RHO_W        1000.0    /* Density of water (kg/m^3) at 0C */
#define Cp           1004.0    /* Specific heat at constant pressure of air 
                                  (J/deg/K) */
#define CH_ICE       2100.0e3  /* Volumetric heat capacity (J/(m3*C)) of ice */

#define SECPHOUR     3600      /* seconds per hour */
#define SNOW_DT       5.0      /* Used to bracket snow surface temperatures
                                  while computing the snow surface energy 
                                  balance (C) */
#define SURF_DT       1.0      /* Used to bracket soil surface temperatures 
                                   while computing energy balance (C) */
#define SOIL_DT       0.25      /* Used to bracket soil temperatures while
                                   solving the soil thermal flux (C) */
#define HOURSPERDAY   24        /* number of hours per day */
#define HOURSPERYEAR  24*365    /* number of hours per year */

/***** Physical Constraints *****/
#define MINSOILDEPTH 0.001      /* minimum layer depth with which model can
                                   work (m) */
#define STORM_THRES  0.001      /* thresehold at which a new storm is 
                                   declared */
/***** Define Boolean Values *****/
#ifndef FALSE
#define FALSE 0
#define TRUE !FALSE
#endif

#ifndef WET
#define WET 0
#define DRY 1
#endif

#ifndef SNOW
#define RAIN 0
#define SNOW 1
#endif

#define min(a,b) ( (a < b) ? a : b )
#define max(a,b) ( (a > b) ? a : b )


#define DAYS_PER_YEAR 365.
#define DtoR 0.017453293      /* degrees to radians */
#ifndef PI
#define PI 3.1415927
#endif
#define STEFAN 5.6696e-8      /* Stefan boltzmann constant */
#define SOLAR_CONSTANT 1400.0 /* Solar constant in W/m^2 */
#define SEC_PER_DAY 86400.    /* seconds per day */

/* define constants for saturated vapor pressure curve (kPa) */
#define A_SVP 0.61078
#define B_SVP 17.269
#define C_SVP 237.3

/* define constants for penman evaporation */
#define CP_PM 1013       /* specific heat of moist air J/kg/C 
                            (Handbook of Hydrology) */
#define PS_PM 101300     /* sea level air pressure in Pa */
#define LAPSE_PM -0.006  /* environmental lapse rate in C/m */

/***********************************************************
  This structure stores the soil parameters for a grid cell.
  ***********************************************************/
typedef struct{
  int rank; 
  int npes;
#if (defined SPMD)
  MPI_Datatype MPI_VEG_CON; 
  MPI_Datatype MPI_SOIL_CON; 
  MPI_Datatype MPI_VEG_LIB;
  MPI_Datatype MPI_OUT_STRUCT; 
  MPI_Datatype MPI_ATMOS_STRUCT; 
  MPI_Datatype MPI_MS_STRUCT; 
#endif
} par_struct; 

typedef struct{
  int *cdi_array;
  int *cdispls;
} ctile_spmd; 

typedef struct {
  float   Ds;                        /* fraction of maximum subsurface flow 
                                        rate */
  float   Dsmax;                     /* maximum subsurface flow rate 
                                        (mm/day) */
  float   Ksat[MAX_LAYERS];          /* saturated hydraulic  conductivity 
                                        (mm/day) */
  float   Wcr[MAX_LAYERS];           /* critical moisture level for soil 
                                        layer, evaporation is no longer 
                                        affected moisture stress in the 
                                        soil (mm) */
  float   Wpwp[MAX_LAYERS];          /* soil moisture content at permanent 
                                        wilting point (mm) */
  float   Ws;                        /* fraction of maximum soil moisture */
  float   alpha[MAX_NODES];          /* thermal solution constant */
  //  float   annual_prec;               /* annual average precipitation (mm) */
  float   avg_temp;                  /* average soil temperature (C) */
  float   b_infilt;                  /* infiltration parameter */
  float   beta[MAX_NODES];           /* thermal solution constant */
  float   bubble[MAX_LAYERS];        /* bubbling pressure, HBH 5.15 (cm) */
  float   bubble_node[MAX_NODES];    /* bubbling pressure (cm) */
  float   bulk_density[MAX_LAYERS];  /* soil bulk density (kg/m^3) */
  float   c;                         /* exponent */
  float   depth[MAX_LAYERS];         /* thickness of each soil moisture 
                                        layer (m) */
  float   dp;                        /* soil thermal damping depth (m) */
  float   dz_node[MAX_NODES];        /* thermal node thickness (m) */
  float   depth_node[MAX_NODES];     /* thermal node depths (m) */
  float   expt[MAX_LAYERS];          /* pore-size distribution per layer, 
                                        HBH 5.15 */
  float   expt_node[MAX_NODES];      /* pore-size distribution per node */
  float   gamma[MAX_NODES];          /* thermal solution constant */
  float   init_moist[MAX_LAYERS];    /* initial layer moisture level (mm) */
  float   max_infil;                 /* maximum infiltration rate */
  float   max_moist[MAX_LAYERS];     /* maximum moisture content (mm) per 
                                        layer */
  float   max_moist_node[MAX_NODES]; /* maximum moisture content (mm/mm) per 
                                        node */
  //float   phi_s[MAX_LAYERS];         /* soil moisture diffusion parameter 
  /**(mm/mm) */
  float   porosity[MAX_LAYERS];      /* porosity (fraction) */
  float   quartz[MAX_LAYERS];        /* quartz content of soil (fraction) */
  float   resid_moist[MAX_LAYERS];   /* residual moisture content of soil 
                                        layer */
  float   rough;                     /* soil surface roughness (m) */
  float   snow_rough;                /* snow surface roughness (m) */
  float   soil_density[MAX_LAYERS];  /* soil partical density (kg/m^3) */
  //  float  *AreaFract;                 /* Fraction of grid cell included in 
  //each elevation band */
  //  float  *Pfactor;                   /* Change in Precipitation due to 
  //                                        elevation (fract) */
  //float  *Tfactor;                   /* Change in temperature due to 
  //                                      elevation (C) */
  float    elevation;                 /* grid cell elevation (m) */
  //float    lat;                       /* grid cell central latitude */
  //float    lng;                       /* grid cell central longitude */
  //float    time_zone_lng;             /* central meridian of the time zone */
//<kluge>
//  float  **layer_node_fract;          /* fraction of all nodes within each 
//                                         layer */
  float layer_node_fract[MAX_LAYERS+1][MAX_NODES];
//</kluge>
  int      FS_ACTIVE;                 /* if TRUE frozen soil algorithm is 
                                         active in current grid cell */

  //  int      gridcel;                   /* grid cell number */
} soil_con_struct;



/***************************************************************************
   This structure stores the atmospheric forcing data for each model time 
   step for a single grid cell.  Each array stores the values for the 
   SNOW_STEPs during the current model step and the value for the entire model
   step.  The latter is referred to by array[NR].  Looping over the SNOW_STEPs
   is done by for (i = 0; i < NF; i++) 
***************************************************************************/
typedef struct {
  /**char   snowflag;  TRUE if there is snowfall in any of the snow 
                       bands during the timestep, FALSE otherwise*/
  float air_temp;  /* air temperature (C) */
  float density;   /* atmospheric density (kg/m^3) */
  float longwave;  /* incoming longwave radiation (W/m^2) (net incoming 
                      longwave for water balance model) */
  float out_prec;      /* Total precipitation for time step - accounts
                          for corrected precipitation totals */
  float prec;      /* average precipitation in grid cell (mm) */
  float pressure;  /* atmospheric pressure (kPa) */
  float shortwave; /* incoming shortwave radiation (W/m^2) */
  float vp;        /* atmospheric vapor pressure (kPa) */
  float vpd;       /* atmospheric vapor pressure deficit (kPa) */
  float wind;      /* wind speed (m/s) */
} atmos_data_struct;

/*******************************************************************
  This structure stores information about the vegetation coverage of
  the current grid cell.
  *******************************************************************/
typedef struct {
  //float  Cv;               /* fraction of vegetation coverage */ 
  //float  Cv_sum;           /* total fraction of vegetation coverage */
  float   zone_depth[2];       /* depth of root zone (hard coded for now..) */
  float   zone_fract[2];       /* fraction of roots within root zone */
  float   root[MAX_LAYERS]; /* percent of roots in each soil layer (fraction) */
  int     veg_class;        /* vegetation class reference number */
  //int     vegetat_type_num; /* number of vegetation types in the grid cell */
} veg_con_struct;

/***************************************************************
  This structure stores all soil variables for each layer in the
  soil column.
  ***************************************************************/
typedef struct {
  float Cs;                /* average volumetric heat capacity of the 
                              current layer (J/m^3/K) */
  float T;                 /* temperature of the unfrozen sublayer (C) */
  float evap;              /* evapotranspiration from soil layer (mm) */
  float ice;               /* ice content of the frozen sublayer (mm) */
  float kappa;             /* average thermal conductivity of the current 
                              layer (W/m/K) */
  float moist;             /* moisture content of the unfrozen sublayer 
                              (mm) */
  float phi;               /* moisture diffusion parameter */
} layer_data_struct;

/******************************************************************
  This structure stores soil variables for the complete soil column 
  for each grid cell.
  ******************************************************************/
typedef struct {
  float aero_resist[3];               /* aerodynamic resistane (s/m) 
                                         [0] = over bare vegetation or soil 
                                         [2] = over snow */
  float baseflow;                     /* baseflow from current cell (mm/TS) */
  float inflow;                       /* moisture that reaches the top of 
                                         the soil column (mm) */
  float runoff;                       /* runoff from current cell (mm/TS) */
  layer_data_struct layer[MAX_LAYERS]; /* structure containing soil variables 
                                          for each layer (see above) */
} cell_data_struct;

typedef struct{
  float swnet;
  float lwnet;
  float qle; 
  float qh; 
  float qg; 
  float qfz; 
  float snowf; 
  float rainf; 
  float evap; 
  float qs; 
  float qsb;
  float snowt;
  float avgsurft;
  float radt;
  float albedo;
  float soilt[MAX_LAYERS];
  float moist[MAX_LAYERS];
  float moist_prev; 
  float swe_prev; 
  float tveg; 
  float esoil; 
  float soilwet; 
  float rootmoist; 
  float swe;
  float qsm; 
  float acond; 
  int   count;
} out_data_struct; 

/***********************************************************************
  This structure stores energy balance components, and variables used to
  solve the thermal fluxes through the soil column.
  ***********************************************************************/
typedef struct {
  char    frozen;                /* TRUE = frozen soil present */
  float  Cs[2];                 /* heat capacity for top two layers 
                                   (J/m^3/K) */
  float  Cs_node[MAX_NODES];    /* heat capacity of the soil thermal nodes 
                                   (J/m^3/K) */
  float  T[MAX_NODES];          /* thermal node temperatures (C) */
  float  Trad[2];               /* surface temperature of energy balance 
                                   (C) */
  float  advection;             /* advective flux (Wm-2) */
  float  albedo;                /* surface albedo (fraction) */
  float  deltaCC;               /* change in snow heat storage (Wm-2) */
  float  deltaH;                /* change in soil heat storage (Wm-2) */
  float  error;                 /* energy balance error (W/m^2) */
  float  fdepth[MAX_FRONTS];    /* all simulated freezing front depths */
  float  grnd_flux;             /* ground heat flux (Wm-2) */
  float  ice[MAX_NODES];        /* thermal node ice content */
  float  kappa[2];              /* soil thermal conductivity for top two 
                                   layers (W/m/K) */
  float  kappa_node[MAX_NODES]; /* thermal conductivity of the soil thermal 
                                   nodes (W/m/K) */
  float  latent;                /* net latent heat flux (Wm-2) */
  float  longwave;              /* net longwave flux (Wm-2) */
  float  moist[MAX_NODES];      /* thermal node moisture content */
  float  refreeze_energy;       /* energy used to refreeze the snowpack 
                                   (Wm-2) */
  float  sensible;              /* net sensible heat flux (Wm-2) */
  float  shortwave;             /* incoming shortwave heat (Wm-2) */
  float  snow_flux;             /* thermal flux through the snow pack 
                                   (Wm-2) */
  float  tdepth[MAX_FRONTS];    /* all simulated thawing front depths */
  float  unfrozen;              /* frozen layer water content that is 
                                   unfrozen */
  int     Nfrost;                /* number of simulated freezing fronts */
  int     Nthaw;                 /* number of simulated thawing fronts */
  //  int     T1_index;              /* soil node at the bottom of the top layer */
} energy_bal_struct;

/***********************************************************************
  This structure stores vegetation variables for each vegetation type in 
  a grid cell.
  ***********************************************************************/
typedef struct {
  float canopyevap;  /* evaporation from canopy (mm/TS) */
  float throughfall; /* water that reaches the ground through 
                        the canopy (mm/TS) */
  float Wdew;        /* dew trapped on vegetation (mm) */
} veg_var_struct;

/************************************************************************
  This structure stores snow pack variables needed to run the snow model.
  ************************************************************************/
typedef struct {
  char MELTING;            /* flag indicating that snowpack melted previously */
  int    snow;             /* TRUE = snow, FALSE = no snow */
  float Qnet;              /* New energy at snowpack surface */
  float albedo;            /* snow surface albedo (fraction) */
  float canopy_vapor_flux; /* depth of water evaporation, sublimation, or 
                              condensation from intercepted snow (m) */
  float coldcontent;       /* cold content of snow pack */
  float coverage;          /* fraction of snow band that is covered with 
                              snow */
  float density;           /* snow density (kg/m^3) */
  float depth;             /* snow depth (m) */
  //float mass_error;        /* snow mass balance error */
  float melt; 

  float pack_temp;         /* depth averaged temperature of the snowpack 
                              (C) */
  float pack_water;        /* liquid water content of the snow pack (m) */
  float snow_canopy;       /* amount of snow on canopy (m) */
  float surf_temp;         /* depth averaged temperature of the snow pack 
                              surface layer (C) */
  float surf_water;        /* liquid water content of the surface layer (m) */
  float swq;               /* snow water equivalent of the entire pack (m) */
  float tmp_int_storage;   /* temporary canopy storage, used in snow_canopy */
  float vapor_flux;        /* depth of water evaporation, sublimation, or 
                              condensation from snow pack (m) */
  int    last_snow;         /* time steps since last snowfall */
} snow_data_struct;         /* an array of size Nrec */

/*****************************************************************
  This structure stores all variables needed to solve, or save 
  solututions for all versions of this model.  Vegetation and soil
  variables are created for both wet and dry fractions of the grid
  cell (for use with the distributed precipitation model).
*****************************************************************/
typedef struct {
  cell_data_struct  **cell[2];    /* Stores soil layer variables (wet and 
                                     dry) */
  //float             *mu;         /* fraction of grid cell that receives 
  //                                  precipitation */
  energy_bal_struct **energy;     /* Stores energy balance variables */
  snow_data_struct  **snow;       /* Stores snow variables */
  veg_var_struct    **veg_var[2]; /* Stores vegetation variables (wet and 
                                     dry) */
} dist_prcp_struct;
/*************************************************************************
  This structure stores information about the time and date of the current
  time step.
  *************************************************************************/
typedef struct {
  int day;                      /* current day */
  int day_in_year;              /* julian day in year */
  int hour;                     /* beginning of current hour */
  int month;                    /* current month */
  int year;                     /* current year */
  float dt; 
} dmy_struct;                   /* array of length nrec created */
/******************************************************************
  This structure stores parameters for individual vegetation types.
  ******************************************************************/
typedef struct {
  char  overstory;        /* TRUE = overstory present, important for snow 
                             accumulation in canopy */
  float rarc;             /* architectural resistance (s/m) */
  float rmin;             /* minimum stomatal resistance (s/m) */
  float LAI[12];          /* monthly leaf area index */
  float Wdmax[12];        /* maximum monthly dew holding capacity (mm) */
  float albedo[12];       /* vegetation albedo (added for full energy) 
                             (fraction) */
  float displacement[12]; /* vegetation displacement height (m) */
  float emissivity[12];   /* vegetation emissivity (fraction) */
  float rad_atten;        /* radiation attenuation due to canopy, 
                             default = 0.5 (N/A) */
  float roughness[12];    /* vegetation roughness length (m) */
  float trunk_ratio;      /* ratio of trunk height to tree height, 
                             default = 0.2 (fraction) */
  float wind_atten;       /* wind attenuation through canopy, 
                             default = 0.5 (N/A) */
  float wind_h;           /* height at which wind is measured (m) */
  float root_depth[2];    /* root depths(m) */
  float root_frac[2];     /* root fractions */
  float  RGL;              /* Value of solar radiation below which there 
                              will be no transpiration (ranges from 
                              ~30 W/m^2 for trees to ~100 W/m^2 for crops) */
  int    veg_class;        /* vegetation class reference number */
} veg_lib_struct;

typedef struct {
#define MS_NVEG ( 2 )
   /* Note, in write_model_state, dist and veg are hard-coded to 1 */
   char MELTING[MS_NVEG][MAX_BANDS];              /* [veg][band] */
   int last_snow[MS_NVEG][MAX_BANDS];             /* [veg][band] */

   float moist[1][MS_NVEG][MAX_BANDS][MAX_LAYERS];/* [dist][veg][band][layer] */
   float ice[1][MS_NVEG][MAX_BANDS][MAX_LAYERS];  /* [dist][veg][band][layer] */
   float Wdew[1][MS_NVEG][MAX_BANDS];             /* [dist][veg][band] */
   float swq[MS_NVEG][MAX_BANDS];                 /* [veg][band] */
   float surf_temp[MS_NVEG][MAX_BANDS];           /* [veg][band] */
   float pack_temp[MS_NVEG][MAX_BANDS];           /* [veg][band] */
   float density[MS_NVEG][MAX_BANDS];             /* [veg][band] */
   float snow_canopy[MS_NVEG][MAX_BANDS];         /* [veg][band] */
   float pack_water[MS_NVEG][MAX_BANDS];          /* [veg][band] */
   float surf_water[MS_NVEG][MAX_BANDS];          /* [veg][band] */
   float T[MS_NVEG][MAX_BANDS][MAX_NODES];        /* [veg][band][node] */

   float coverage[MS_NVEG][MAX_BANDS];            /* [veg][band] */
   float coldcontent[MS_NVEG][MAX_BANDS];         /* [veg][band] */

} model_state_struct;

