/*.............................................................................
encoder_ex1.c

This example shows how to manually initialize the 3 Encoder structures
then use the GRIB library to encode a GRIB message.

INPUT:   $GRIB_ENV/data/IEEE.input    file contain Rows*Cols floating pt values

OUTPUT:  '075_237_1997070100012_011_105.00002.0.grb'
     
recommended GRIB filename convention:
	     'MID_GID_yyyymmddhhttt_PID_LID_lvl1.c.grb' 
     where 	
	MID	= 3-digit model id
	GID	= 3-digit geometry id
	yyyy	= 4-digit year of reference date/time
	mm	= 2-digit month of reference date/time
	dd	= 2-digit day of reference date/time
	hh	= 2-digit hour of reference date/time
	ttt	= 3-digit forecast period relative to reference date/time
	PID	= 3-digit Parameter id
	LID	= 3-digit level id
	lvl1	= 5-digit Level 1 value scaled to integer
	c 	= 1-digit alphanumeric Case Identifier as in USER_INPUT struct

Revisions:
   18jun97 first version;
.............................................................................*/
#include <stdio.h>
#include <stdlib.h>
#include "gribfuncs.h"	    /* Library function propotypes */
#define  INPUT_IEEE_FN	    "data/IEEE.input"

main ()
{
FILE	    *fieee;	    /* file pointer to IEEE data file */
GRIB_HDR    *gh= 0;         /* will hold encoded msg & its info */
DATA_INPUT  data_input;     /* input header structure for Encoder */
GEOM_IN     geom_in;        /* geometry description for Encoder */
USER_INPUT  user_input;     /* user input from input.dat for Encoder */
char	    *grib_env;	    /* working variable */
char	    input_fn[200];  /* full path to the ieee file */
char        errmsg[2000];   /* buffer to hold error message */
char	    out_fn[50];	    /* name of output filename */
float	    *flt_arr=NULL;  /* array to hold float data, init to null */
float	    min, max;	    /* min and max value of data */
int	    cnt,i;	    /* working variables */


    fprintf(stdout,"Clear out the 3 Encoder structures \n");
    init_enc_struct (&data_input, &geom_in, &user_input);

    fprintf(stdout, "Manually fill GEOM_IN structure \n");
    strcpy (geom_in.prjn_name, "spherical");
    geom_in.nx= 61;		/* num of cols */
    geom_in.ny= 51;		/* num of rows */
    geom_in.x_int_dis= 22.;	/* X-dir grid length, in meters */
    geom_in.y_int_dis= 22.;	/* Y-dir grid length, in meters */
    geom_in.parm_1= 0.2;	/* Latitude spacing */
    geom_in.parm_2= 0.2;	/* Longitude spacing */
    geom_in.parm_3= -1.0;   	/* Unused for Spherical */
    geom_in.first_lat= 29.;   	/* Lat of 1st pt, in degrees */
    geom_in.first_lon= -126.5; 	/* Lon of 1st pt, in degrees */
    geom_in.last_lat= 39.;      /* Lat of last pt, in degrees */
    geom_in.last_lon= -114.5;   /* Lon of last pt, in degrees */
    geom_in.scan= 64; /* pts scan in +i, +j, adjcent pts in i-dir consecutive */
    geom_in.usRes_flag= 0;   /* Earth spherical, UV rel. to East/Northerly Dir*/

   fprintf(stdout, "Manually fill DATA_INPUT structure \n");
   data_input.usProc_id= 75;       /* Model/Generating Process ID (TabA) */
   data_input.usGrid_id= 237;      /* Grid Identification (Table B) */
   data_input.usParm_id= 11;        /* GRIB parameter id */
   data_input.usParm_sub_id= 0;    /* GRIB parameter sub-id */
   data_input.usLevel_id= 105;     /* GRIB level id */
   data_input.nLvl_1=   2;         /* 1st level value - scaled to an integer*/
   data_input.nLvl_2=   0;         /* 2nd level value - scaled to an integer*/
   data_input.nYear= 1997;         /* year of data e.g. 1993 */
   data_input.nMonth=  7;          /* month of year e.g. 8 */
   data_input.nDay=    1;          /* day of month e.g. 31 */
   data_input.nHour=   0;          /* hour of day e.g. 0 */
   data_input.nMinute=  0;         /* minute of hour e.g. 0 */
   data_input.nSecond=  0;         /* second of minute e.g. 0 */
   data_input.usFcst_id=   1;      /*  HOURS Forecast time unit id - Table 4 */
   data_input.usFcst_per1= 12;     /* forecast time 1 (tau) e.g. 0. */
   data_input.usFcst_per2= 0;      /* forecast time 2 (tau) e.g. 0. */
   data_input.usTime_range_id=  0; /* Time range indicator - Table 5 */
   data_input.usTime_range_avg= 0;/* Number in average */
   data_input.usTime_range_mis= 0;/* Number missing from average */
   data_input.nDec_sc_fctr=  1;   /* Decimal scale factor */

    fprintf(stdout, "Manually fill USER_INPUT structure\n");
    user_input.chCase_id= '0'; /* User defined Case ID (1 digit alphanumeric)*/
    user_input.usParm_tbl= 2;  /* GRIB Table Version Number             */
    user_input.usSub_tbl= 1;   /* Local Table Version Number            */
    user_input.usCenter_id= 128;  /* ID of Originating Center (Table 0)     */
    user_input.usCenter_sub= 1;   /* Sub-Table Entry for originatingCtr (Tbl 0)*/
    user_input.usTrack_num= 0;    /* Tracking ID for data set               */
    user_input.usBDS_flag=  0;    /* Binary Data Section Flag (Table 11)    */
    user_input.usGds_bms_id= 128; /* GDS present but not BMS */
    user_input.usBit_pack_num= 0; /* No. of bits into which data is packed, */
				  /* or Zero to encode in Least num of bits */
 
  /* 
     MAKE a float array large enough to hold grid size of this geom; 
     Fill array with the #Rows by #Cols elements of Float values  
   */ 
   fprintf(stdout,"Malloc float array %d by %d\n", geom_in.nx, geom_in.ny);
   if (!(flt_arr = (float *) malloc(geom_in.nx * geom_in.ny  *sizeof(float))))
       { fprintf(stderr,"Failed to malloc Float array\n");
	 exit(1);
  	}
    
    grib_env = getenv ("GRIB_ENV");
    if (grib_env == NULL || *grib_env == '\0') {
	fprintf(stderr,"Environment variable GRIB_ENV not defined\n");
	exit(1);
	}

    sprintf (input_fn, "%s/%s", grib_env, INPUT_IEEE_FN);
    fieee = fopen (input_fn, "rb");
    if (fieee == NULL) {
	fprintf(stderr,"Failed to open '%s'\n", input_fn);
	if (flt_arr != NULL) free(flt_arr);
	exit(1);
	}
    else { 
	    fprintf(stdout,
	    "Read in (%d X %d) float data elements from '%s'\n",
	     geom_in.nx, geom_in.ny, input_fn);

	    if (fread((void*)flt_arr, sizeof(float), geom_in.nx*geom_in.ny, 
		fieee) != geom_in.nx*geom_in.ny) 
	    {
	   	fprintf(stderr,"Failed to read 'ieee.input'\n");
        	if (flt_arr != NULL) free(flt_arr);
		fclose (fieee);
        	exit(1);
             }

	/* print out a sample of the data */
        for (i=0; i < 100; i++) 
	    if (i % 5 == 4)
		 fprintf(stdout,"%.*lf\n",
		  (data_input.nDec_sc_fctr>0? data_input.nDec_sc_fctr: 0), 
		  flt_arr[i]); 
	    else fprintf(stdout,"%.*lf ",
		  (data_input.nDec_sc_fctr>0? data_input.nDec_sc_fctr: 0), 
		  flt_arr[i]); 
		  
	/* Just to see what Range of Data is, and also count #Zeros too;  */
        cnt = 0;
        min = max = flt_arr[0];
	for (i=geom_in.nx*geom_in.ny-1; i > 0; i--) {
	    if (min > flt_arr[i])  min = flt_arr[i]; 
	    if (max < flt_arr[i])  max = flt_arr[i]; 
	    if (flt_arr[i] == 0.0) cnt++;
	  }

	fprintf(stdout,"--> MIN value= %lf, MAX= %lf,  %ld Zeros\n", 
	min, max, cnt);
	fclose (fieee); 
	} 
  
  fprintf(stdout,"Allocate storage and initialize GRIB_HDR structure\n");
  if (init_gribhdr (&gh, errmsg))  {
        fprintf(stderr,"Abort; error=%s\n", errmsg);
	if (flt_arr != NULL) free(flt_arr);
	exit(1);
	}

   /* go encode message now
   */
   fprintf(stdout,"Encode message now\n");
   if ( grib_enc(data_input, user_input, geom_in, flt_arr, gh, errmsg) != 0){
	fprintf(stderr,"Abort; error=%s\n",errmsg);
	if (flt_arr != NULL) free(flt_arr);
	free_gribhdr (&gh);
        exit(1);
        }
   fprintf(stdout,"Encoder complete, message is in GRIB_HDR\n\n");

   /* If desired, display content of GRIB_HDR struct 

   	display_gribhdr (gh);
   */

   /* form a filename for the output file which reflects
      the type of message this is
   */
   make_default_grbfn (data_input, user_input, out_fn);

   /* Go save message in GRIB_HDR struct out to Output file
   */
   fprintf(stdout,"Save encoded message in '%s'\n"\
   "output file length should be %ld bytes\n", out_fn,gh->msg_length);

   if ( gribhdr2file (gh, out_fn, errmsg) != 0) 
	{
        fprintf(stderr,"Abort; error= %s\n", errmsg);
	if (flt_arr != NULL) free(flt_arr);
	free_gribhdr (&gh);
        exit(1);
        }

   if (flt_arr != NULL) free(flt_arr);  /* release float array */
   free_gribhdr (&gh);			/* release grib header storage */
   exit(0);
}
