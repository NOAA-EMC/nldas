/*.............................................................................
File:  encoder_ex2.c

Example program to encode one message, loading all information from
external ascii files using library functions.  A GRIB edition 1
message is then created using the grib_enc library function and 
written to file using the default file name.

INPUT FILES:  
	  CONFIG_FN       "$GRIB_ENV/config/encoder.config"
	  GEOM_FN         "$GRIB_ENV/data/encoder_ex2.geom"
	  INFO_FN	  "$GRIB_ENV/data/encoder_ex2.info"
	  IEEE_FN	  "$GRIB_ENV/data/IEEE.input"

OUTPUT:  '075_237_1997070100012_011_105.00002.1.grb'

Recommended Naming convention for GRIB files=
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
	c 	= 1-digit alphanumeric Case Identifier as in USER_INPUT struct;
	.grb    = string, as is

Revisions:
   20jun97 first version;
.............................................................................*/
#include <stdio.h>
#include <stdlib.h>
#include "gribfuncs.h"	  /* Library library function prototypes */

#define  CONFIG_FN        "config/encoder.config"
#define  GEOM_FN          "data/encoder_ex2.geom"
#define  INFO_FN	  "data/encoder_ex2.info"
#define  IEEE_FN	  "data/IEEE.input"

main ()
{
GRIB_HDR    *gh= 0;         /* will hold encoded msg & its info */
DATA_INPUT  data_input;     /* input header structure for Encoder */
GEOM_IN     geom_in;        /* geometry description for Encoder */
USER_INPUT  user_input;     /* user input from input.dat for Encoder */
char        errmsg[2000];   /* buffer to hold error message */
char	    input_config[200];/* name of input filename */
char	    input_geom[200];  /* name of input filename */
char	    input_info[200];  /* name of input filename */
char	    input_ieee[200];  /* name of input filename */
char	    *grib_env;	    /* working variable */
char	    out_fn[50];	    /* name of output filename */
float	    *flt_arr=NULL;  /* array to hold float data, init to null */
float	    min, max;	    /* min and max value of data */
int	    cnt,i;	    /* working variables */

    grib_env = getenv ("GRIB_ENV");
    if (grib_env == NULL || *grib_env == '\0') {
	fprintf(stderr,"Environment variable GRIB_ENV not defined\n");
	exit(1);
	}

    fprintf(stdout,"Clear out the 3 Encoder structures \n");
    init_enc_struct (&data_input, &geom_in, &user_input);
    
    fprintf(stdout,"call ld_enc_geomfile, ld_enc_ffinfo, ld_enc_config\n");
    sprintf (input_config, "%s/%s", grib_env, (char *)CONFIG_FN);
    sprintf (input_geom, "%s/%s", grib_env, (char *)GEOM_FN);
    sprintf (input_info, "%s/%s", grib_env, (char *)INFO_FN);
    sprintf (input_ieee, "%s/%s", grib_env, (char *)IEEE_FN);

    if (ld_enc_geomfile (input_geom, &geom_in, errmsg) != 0 ||
        ld_enc_ffinfo (input_info,  &data_input, errmsg) != 0 ||
        ld_enc_config (input_config, &user_input, errmsg) != 0 )
        {
	fprintf(stderr,"Fatal error= %s\n", errmsg);
	exit(1);
	}
 
    fprintf(stdout,"Malloc float array %d by %d\n", geom_in.nx, geom_in.ny);
    if (!(flt_arr = (float *) malloc(geom_in.nx * geom_in.ny *sizeof(float))))
       { 
	fprintf(stderr,"Failed to malloc Float array\n");
	 exit(1);
  	}

    fprintf(stdout,"call ld_enc_ieeeff\n");
    if (ld_enc_ieeeff (input_ieee, flt_arr, geom_in.nx * geom_in.ny, errmsg))
	{
        fprintf(stderr,"Failed to load IEEE file, ermsg=%s\n", errmsg);
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

