#ifndef PREPROC_H
#define PREPROC_H

/*************************************************************/
/* header file:	ncep_forcing_pp.h                            */
/*************************************************************/

#include "user_def.h"
#include "vicNl.h"
#include "vicNl_ldas.h"

/*************************************************************/
/* #define statements                                        */
/*************************************************************/

/*************************************************************/
/* compilation flags                                         */
/*************************************************************/
#define VERBOSE		0
#define DEBUG		0
#define VERBOSE_QC 1
#define DEBUG_GRIB 0	/* for grib routines       */

#define FILL_VALUE      -9999.00       /* default for missing data pts in grib */

#define INS2AVE			1

#define MET_ATTR_USE		0
#define MET_ATTR_DESC		1
#define MET_ATTR_FNAME		2
#define MET_ATTR_GRIBPROC	3
#define MET_ATTR_GRIBID		4
#define MET_ATTR_TYPE		5
#define MET_ATTR_MINBOUND	6
#define MET_ATTR_MAXBOUND	7
#define MET_ATTR_MAX	8
/*#define MET_VARS 16  */           /* number of met types contained in NCEP
				  grib files */

#if NARR_RETRO
enum MET_VARS {
	MET_VAR_TMP=0,
	MET_VAR_SPFH,
	MET_VAR_PRES,
	MET_VAR_UGRD,
	MET_VAR_VGRD,
	MET_VAR_DLWRF,
	MET_VAR_RAINSNOW,
	MET_VAR_CAPE,
	MET_VAR_PEVAP,
	MET_VAR_APCP,
	MET_VAR_DSWRF,
        MET_VAR_MAX        
};
#else
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
#endif

/*************************************************************/
/* structure definitions                                     */
/*************************************************************/

/* struct for statistics */
typedef struct {
    int nsteps;
    int ngribs;
    int nmatch;
    int nmiss;
    int nsubs;
    int nblended;
    int nins2ave;
    int nwritten;
    int match[MET_VAR_MAX];
    int miss[MET_VAR_MAX];
    int subs[MET_VAR_MAX];
    int blended[MET_VAR_MAX];
    int ins2ave[MET_VAR_MAX];
    int written[MET_VAR_MAX];
    long nqcfails;
    long qcfails[MET_VAR_MAX];
} Stats;

/***********************************/
/* header for global options file  */
/***********************************/


/*************************************************************/
/* function prototypes                                       */
/*************************************************************/

/* awrite grib met */
int read_grib(char *fname, int id, int proc, float **grib_data, int reset);
void write_grib( float *grib_data, int *cell_num, int *cell_index, int ncell_comp, 
		 int format, FILE *fp, int var, int rec, FILE *qcp, float min, float max, char *fname, Stats *stats );
/* qc data checks */
int qc_data(FILE *qcp, float *grib_data, int *cell_num, int ncell_comp, float min_bound, float max_bound);

/* process precip */
void process_precip(ldas_index_struct, float *, float *, int var, Stats *stats);
void process_dswrf(ldas_index_struct, float *, float *, int var, Stats *stats);
void process_inst2ave(ldas_index_struct ldas_index, float *grib_data, float *grib_data_old, int var, Stats *stats);

#endif

