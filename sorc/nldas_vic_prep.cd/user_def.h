#ifndef USER_DEF_H
#define USER_DEF_H

/**********************************************************************
  This header file contains model parameters that can be modified by
  the user to control model performance.  When this file is modified 
  the model needs to be recompiled for the changes to take effect.
**********************************************************************/

/***** If TRUE include all model messages to stdout, and stderr *****/
#define VERBOSE FALSE

/***** If TRUE limit output data to runoff and baseflow for optimization *****/
#define OPTIMIZE FALSE

/***** If TRUE include all debugging code - debugging options still
       have to be activated to get extra output.  When set to FALSE
       all debugging if-then statements are removed from the compiled 
       code *****/
#define LINK_DEBUG TRUE

/***** If TRUE output will be in LDAS binary format, which is a single
       file with limited variables, most of which are truncated to
       conserve disk space *****/
#define LDAS_OUTPUT FALSE

/***** If TRUE VIC uses a system of linear equations defined in global.h
       to estimate the maximum unfrozen water content equation.  This 
       significantly reduces the run time with frozen soil, but may
       introduce new errors (STILL UNDER TESTING) *****/
#define QUICK_FS FALSE
#define QUICK_FS_TEMPS 7

/***** If TRUE VIC uses the linear interpolation of the logarithm of the
       matric potential from the two surrounding layers to estimate the 
       soil moisture drainage from each layer (Boone and Wetzel, 1996).
       This should improve the soil moisture drainage predicted by the
       low resolution solution computed by VIC. *****/
#define LOW_RES_MOIST FALSE

/***** If TRUE VIC code to save the model state is included in the 
       compiled code.  If STATEYEAR, STATEMONTH and STATEDAY are
       defined in the global control file the model state will
       be written to a file. *****/
#define SAVE_STATE TRUE

/***** Define maximum array sizes for model source code *****/
#define MAX_VEG      10          /* maximum number of vegetation types per 
				   cell */
#define MAX_LAYERS   3          /* maximum number of soil moisture layers */
#define MAX_NODES    18         /* maximum number of soil thermal nodes */
#define MAX_BANDS    10          /* maximum number of snow bands */
#define MAX_FRONTS   3          /* maximum number of freezing and thawing 
				   front depths to store */

/***** Number of iterations to use in solving the surface energy balance.
       The original VIC model uses only 1 iteration for speed.  Increasing
       the number of iterations improves precision, and is recommended
       for single point comparisons with frozen soils *****/
#define MAXIT_FE        25
 
/***** Coefficient multiplied by the LAI to determine the amount of
       water that can be storde in the canopy *****/
#define LAI_WATER_FACTOR 0.2

/***** Longwave correction factor, used to correct estimated incoming
       longwave radiation (use 1, unless measured longwave available for
       calibration) *****/
#define LWAVE_COR	1.

/***** Snow albedo curve parameters.  Defaults are from Bras p263.
       Should not be changed except for serious problems with snow melt *****/
#define NEW_SNOW_ALB		0.85
#define SNOW_ALB_ACCUM_A	0.94
#define SNOW_ALB_ACCUM_B	0.58
#define SNOW_ALB_THAW_A		0.82
#define SNOW_ALB_THAW_B		0.46

#endif

