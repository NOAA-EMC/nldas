/****************************************************************************
*   File:  gsv5d.h
*   V5D related structures
****************************************************************************/
#include <stdio.h>



#ifdef V5D_ON
#include "v5d.h"	/* comes with the Vis5D Software */
#include "binio.h"	/* comes with the Vis5D Software */

#else
#define MAXVARS     1	/* Vis5D Software not avail, so 	*/
#define MAXTIMES    1   /* set up dummy constants 		*/
#define MAXROWS     1
#define MAXCOLUMNS  1
#define MAXROWS     1
#define MAXCOLUMNS  1
#define MAXLEVELS   1
#define MAXPROJARGS 1
#define MAXVERTARGS 1
#define MISSING   0.0
#endif

/*####################################################################
*  THE NEXT STRUCTS ARE USED BY GRIBSIMP TO MAINTAIN vis5d LEVEL INFO,
*  A LINKED LIST OF PARM, A LINKED LIST OF GRIB MESSAGES THAT 
*  QUALIFY TO BE INCLUDED IN THE V5D FILES;
*####################################################################*/


typedef struct V5D_LVL		/* LEVEL STRUCTURE  */
{
  unsigned short usLevel_id;	/* grib level id */
  unsigned short numheights;	/* #heights for this level, max100 */
  unsigned short height[MAXLEVELS];/* array of height vals for chosen lvl */
				/* Max level defined in 'v5d.h' */
} V5D_LVL;


typedef struct V5D_PARM	/* LINKED LIST OF PARAMETER STRUCTS */
{				/* used for 'VAR' info line in .ctl file */
  char		 abbrv[9];	/* 'axxxyyy[n]' x=parmid, y=lvlid, n=a-z */ 
				/* n= not used for ZDEF level     */
  char  	 varnm[41];	/* unabbreviated name & unit */
  unsigned short usParm_id;     /* parm id for this parm */
  struct V5D_LVL *lvl_ptr;	/* pts to the cell within V5D_LVL linked list 
				   which has the level id same as usLevel_id */
  struct V5D_PARM *next;	/* pts to next cell in list */
  struct V5D_PARM *last;	/* pts to next cell in list */
} V5D_PARM;


/*****************************************************/

typedef struct  	/* Common Info for Vis5D */
{
  long   base_dtg;		/* yyyymmddhh, default is 1st msg  */
  double ebase_dtg;		/* base dtg in Epochal Time (hrs) */
  int    nrows;			/* number of rows */
  int    ncols;			/* number of cols */
  int    zdef_lvl;		/* the User-specified level */
  int    v5d_proj;		/* as defined in Vis5D, different from GRIB */
  unsigned short usModel;       /* constant for entire set */
  unsigned short usGeom;        /* constant for entire set */
  float  proj_args[7];    	/* info on projection */
} V5_INFO;

typedef struct V5_REC 	/* LINKED LIST OF Vis5D MSGS */
{
  long     offset;	    /* Byte Offset from beg. of input file */
  double   valid_etime;	    /* Base Time + Fcst Period : yyyymmddhhmm */
  long     base_dtg;	    /* Base Time only: yyyymmddhhmm */
  unsigned short usP1;	    /* Forecast Period */
  unsigned short usHeight1; /* Height value of this msg */
  V5D_PARM *parm_ptr;	/* pts to parm type in parm list */
  struct V5_REC *next;	/* pts to next cell in list */
  struct V5_REC *last;	/* pts to next cell in list */
} V5_REC;

