/****************************************************************************
*   File:  grads.h
*   GRADS program requires 2 input files:  fn.ctl and fn.gmp
*  - Control File ".ctl" is in ascii, shows info of all fields found in .gmp;
*    --> see make_grad_file() to see how it's created;
*  - ".gmp" is the Binary file which contains info on each Field present;
****************************************************************************/

#include <stdio.h>
#define	MAXHTS	100	/* max hts per level id */

/*####################################################################
*  THE NEXT 4 STRUCTS ARE USED BY GRIBSIMP TO MAINTAIN LINKED LISTS OF 
*  LEVEL/PARM INFO AS WELL AS LINKED-LIST OF GRIB MESSAGES THAT 
*  QUALIFY TO BE INCLUDED IN THE GRADS CONTROL FILES;
*####################################################################*/

typedef struct  	/* Common Info for Grads */
{
  long  base_dtg;		/* yyyymmddhh, default is 1st msg  */
  double ebase_dtg;		/* base dtg in Epochal Time (hrs) */
  long ulGrid_size;		/* size of grid */
  int  zdef_lvl;		/* level to enumerate */
} GRAD_INFO;


typedef struct GRAD_LVL		/* LINKED LIST OF LEVEL STRUCTS */
{
  unsigned short usLevel_id;	/* grib level id */
  unsigned short numheights;	/* #heights for this level, max100 */
  unsigned short height[MAXHTS];/* array of heights vals for this lvl;  if
				   non-ZDEF level then only Height[0] is used */
  struct GRAD_LVL *next;	/* pts to next cell in list */
  struct GRAD_LVL *last;	/* pts to next cell in list */
} GRAD_LVL;


typedef struct GRAD_PARM	/* LINKED LIST OF PARAMETER STRUCTS */
{				/* used for 'VAR' info line in .ctl file */
  char		 abbrv[10];	/* 'Pxxxxyyy[z]' 
				    P  : Main/sub parm table 'm'/'a'-'e'
				    xxx: 3 digit parmid
				    yyy: 3 digit level id
				    z  : 'a'-'z', not used for ZDEF level     */
  char  	 varnm[41];	/* unabbreviated name & unit */
  unsigned short usParm_id;     /* combo of parmid & parmsubid  for this parm */
  unsigned short usLevel_id;    /* level id for this parm.  If level equals
				   the ZDEF level then usHeight is not used;
				   in that case, actual Heights values are 
				   stored in that GRAD_LVL's height[] array; 
				   curr_gradparm->lvl_ptr->height[0:#heights] */
  unsigned short usHeight;      /* only used for non-ZDEF levels, else 65535 */
  struct GRAD_LVL *lvl_ptr;	/* pts to the cell within GRAD_LVL linked list 
				   which has the level id same as usLevel_id */
  struct GRAD_PARM *next;	/* pts to next cell in list */
  struct GRAD_PARM *last;	/* pts to next cell in list */
} GRAD_PARM;


typedef struct GRAD_REC 	/* LINKED LIST OF MSGS READ IN */
{
  long  base_dtg;	     	/* Base Time only: yyyymmddhhmm */
  unsigned short ustau;		/* #of time_units offset from basetm */
  unsigned short usHeight1;	/* Top Height of this level type  */

  GRAD_PARM *parm_ptr;		/* pts to parm type in parm list */
  long   tau_incr;		/* increment of time unit from Base dtg */
  int   dpos;			/* #bytes fr BOF to BDS's bitstream */
  int   bpos;			/* #bytes fr BOF to BMS's bitstream */
  int   bnum;			/* number of bits per datapoint */
  float	fDec_sc_fctr;		/* decimal scale factor */
  float fBin_sc_fctr;		/* binary scale factor */
  float fReference;		/* Reference value */ 
  struct GRAD_REC *next;	/* pts to next cell in list */
  struct GRAD_REC *last;	/* pts to next cell in list */
} GRAD_REC;



/*####################################################################
*   ".gmp" is extension of 1 of the input files required by the GrRADS 
*   software package;  its Fwrite format is= 
*    		- GMP_BLK0, 
*    		- (GMP_BLK0.blk1_elements * GMP_BLK1), 
*    		- (GMP_BLK0.blk3_elements * GMP_BLK3), 
*    		- (GMP_BLK0.blk4_elements * GMP_BLK4);
*####################################################################*/
/*
* BLOCK 0:  header info of file, store only once at beg. of file
*/
typedef struct GMP_BLK0 
{
  int   type;       /* Indexing file type, 1 for GRIB */
  int   blk1_elements;  /* Num of Ints expected in Blk1, (always 4)  */
  int   blk2_elements;  /* Num of Floats expected in Blk2, (always 0) */
  int   blk3_elements;  /* Num of Ints expected in Blk3 (3*times*rescpertime) */
  int   blk4_elements;  /* Num of Floats expected Blk4 (3*times*rescpertime) */
  int   *notused1;      /* unknown usage, leave null */
  float *notused2;      /* unknown usage, leave null */
  int   *notused3;      /* unknown usage, leave null */
  float *notused4;      /* unknown usage, leave null */
} GMP_BLK0;

/* 
* BLOCK 1: num of times to store this is GMP_BLK0.blk1_elements 
*/
typedef struct GMP_BLK1 {   /* only one per entire GMP file */
  int  filetype;	/* 1 means grib */
  int  tdef;		/* number of times  */
  int  recspertime;     /* #recs (all Heights counted) per Time */
  int  usGrid_id;       /* Grid Defn Section identification  */
} GMP_BLK1;

/* 
* BLOCK 3: num of times to store this is GMP_BLK0.blk3_elements 
*/
typedef struct GMP_BLK3  {  /* Block of (times*recspertime) records,*/
		            /* 1 record per Grib message */
   int dpos; 		/* where BDS starts w/respect to beg. of file */
   int bpos;		/* Where PDS starts w/respect to beg. of file  */
   int bnum;		/* number of bits per data point */
} GMP_BLK3;

/* 
* BLOCK 4: num of times to store this is GMP_BLK0.blk4_elements 
*/
typedef struct GMP_BLK4  {  /* Block of (times*recspertime) records,*/
		            /* 1 record per Grib message */
   float fDec_sc_fctr;
   float fBin_sc_fctr;
   float fReference;
} GMP_BLK4;
