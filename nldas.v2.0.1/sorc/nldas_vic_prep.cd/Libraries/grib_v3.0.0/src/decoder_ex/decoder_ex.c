/**********************************************************************
File:  decoder_ex.c 					June 28, 1997
Example program showing how one can decode an input file which has  
multiple-GRIB messages by calling the GRIB Library functions directly.

Input:  '$GRIB_ENV/data/GRIB0797.tar'
	a file containing 6 GRIB Edition 1 messages, with filler bytes
	before and after each message.

Output: 'decoder_ex.output' 
	ASCII file contains data values for ALL of the messages found in
	input file.

Revisions:
02/02/98 atn:  change parmid/parmsub printing;
*********************************************************************/ 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "gribfuncs.h"  		/* GRIB library function prototypes */

#define TABLE_PATH      "../tables"    /* default dir of Lookup Tables */ 
#define FILL_VALUE      -9999.00       /* default for missing data pts */


#ifdef PROTOTYPE_NEEDED
void main (int argc, char **argv) 
#else
void main (argc, argv) int argc; char **argv; 
#endif
{
FILE    *fout;	 	/* output file */
char    InFile[200];   	/* input file name */
char    *grib_env_dir;   /* holds value of GRIB_ENV environment variable */
char    errmsg[2000];    /* required, error msgs back from Grib library */
char    OutFn[200]; 	 /* name of output file */
int 	i;		 /* working var */
int	D;		 /* Decimal Scale Factor */
int	decoded_cnt;     /* number of messages decoded successfully */
int     nReturn=0;       /* return status from Grib_dec           */
int     Rd_Indexfile=0;  /* Zero to decode all msgs in input file  */
long    offset;   	 /* byte offset from beginning of file, for inputfiles
			    with multiple GRIB messages; */
float   *grib_data;  	 /* pointer to block of decoded data */

BMS_INPUT	 bms;	/* input structure for bitmap section   */
PDS_INPUT        pds;   /* input structure for product definition section */
grid_desc_sec    gds;   /* input structure for grid description section */
BDS_HEAD_INPUT   bds_head; /* input structure for bds section */
GRIB_HDR 	 *gh1;	/* holds Msg and its info, put there by grib Seek */

/*
* Initializes variables
*/
   errmsg[0] = '\0';		/* error buffer */
   offset= 0L;			/* byte offset */
   decoded_cnt = 0;		/* count of messages */
   fout = (FILE *)NULL;		/* output file */
   grib_data= (float *)NULL;  	/* decoder will create storage */
   
/*
* Set up name of file to decode   ! $GRIB_ENV/data/GRIB0797.tar
*/   
   grib_env_dir= getenv("GRIB_ENV");
   if (grib_env_dir==NULL || *grib_env_dir=='\0') {
	fprintf(stdout,"Error: Environment variable GRIB_ENV not defined\n");
	exit(0);
        }
   else sprintf (InFile, "%s/data/GRIB0797.tar", grib_env_dir);

/*
* MAKE storage for Grib Header !exit on error
*/
   if ((nReturn= init_gribhdr (&gh1, errmsg) ))  goto bail_out; 

/*
* Open output file 'decoder_ex.output' to store ALL message's data (in ASCII)
*/
   sprintf (OutFn, "decoder_ex.output");
   fout = fopen (OutFn, "w");
   if (fout == NULL) { 
	fprintf(stderr,"Failed to open output file %s\n", OutFn); 
	nReturn = 1;  /* error stat */
	goto bail_out;
	} 

/*
* loop until no more messages in InFile 
*/
  for  (offset = 0L; nReturn == 0;  offset += gh1->msg_length) 
  {
     if (nReturn= grib_seek(InFile, &offset, Rd_Indexfile, gh1, errmsg))
	{	
	  fprintf(stdout,"Grib_seek returned non zero stat (%d)\n", nReturn);
	  if (nReturn == 2) break;	/* End of file error */
	  else goto bail_out;		/* abort if other error */
	}

     if (errmsg[0] != '\0') 
   	{ /* NO errors but got a Warning msg from seek */
	  fprintf(stdout,"%s; Skip Decoding...\n",errmsg);  
	  errmsg[0]='\0'; 	/* reset error to continue */
	  gh1->msg_length = 1L; 	/* set to 1 to bump offset up */
	  continue;
	}
	
     if (gh1->msg_length < 0) {
	  fprintf(stderr, "Error:  message returned had bad length (%ld)\n",
	  gh1->msg_length);
	  goto bail_out;
	}
     else if (gh1->msg_length == 0) {  
	  fprintf(stdout,"msg_lenth is Zero, set offset to 1\n");
	  gh1->msg_length = 1L;  /* set to 1 to bump offset up */
	  continue; 
	}

/*
* Clear out the Structures 
*/
      init_dec_struct(&pds,&gds,&bms,&bds_head);
      fprintf(stdout,"Decoding message found at %ld bytes ...\n",offset);

/*
* Go decode the message currently in Grib Header;
* Float array MUST be Null, Decoder will allocate storage for it
* If success, PDS, GDS, BDS_HEAD, BMS, GRIB_DATA are returned filled;
*/
      grib_data= (float *) NULL; 
      if (nReturn = grib_dec ((char *)gh1->entire_msg, 
	   &pds, &gds, &bds_head, &bms, &grib_data, errmsg)) goto bail_out;

      decoded_cnt++;

/*
* IF Bitmap Section is present, go apply BMS to Grib Data float array
*/
      if (bms.uslength>0 && 
	 (nReturn=apply_bitmap(&bms, &grib_data, FILL_VALUE, &bds_head,errmsg)))
	  goto bail_out; 

/*
* .... 	  FLOAT DATA ARRAY IS READY TO GO             ....  
* .... 	  YOU MAY CHOOSE TO DO WHATEVER HERE          ....
* .... 	  PUT YOUR CODE FROM HERE TO <<END MARKER>>   .... 
*/

/*
* Here, let's go print it out
*/
   fprintf(stdout, "\n\t>>>  Content of GRIB Msg #%d =\n", decoded_cnt);
   prt_inp_struct(&pds, &gds, &bms, &bds_head, &grib_data);
/*
* Also display entire data array into ASCII output file
*/
     fprintf(fout,
     "decoded GRIB Msg #%d from $GRIB_ENV/data/GRIB0797.tar\n"\
     "  dtg= %02u%02u%02u%02u Fcstper=%03u\n",
      decoded_cnt, pds.usYear, pds.usMonth, pds.usDay, pds.usHour,pds.usP1);

     if (pds.usExt_flag == EXTENSION_FLAG && pds.usParm_sub != 0)
          /* GRIB Extension used */
          fprintf(fout,"  Parmid in Sub-Tbl '%c' =%03u ",
          pds.usParm_id-250 + 'A', pds.usParm_sub );

     else fprintf(fout,"  Parmid=%03u (Main Tbl), ", pds.usParm_id);

     fprintf(fout,"Levelid=%03u, GridId=%03u, lvl1=%05ld\n",
     pds.usLevel_id, pds.usGrid_id, pds.usHeight1);
     fprintf(fout,"  Decimal Scale Factor =  %d\n", pds.sDec_sc_fctr);

    /* Change local 'D' to zero if the Decimal Scale factor is negative
       Only doing this to work around "%10.-nf" format in the print statement 
    */
    D = ((int) pds.sDec_sc_fctr < 0 ?  0 : (int) pds.sDec_sc_fctr );

    for (i=0; i<bds_head.ulGrid_size; i++) {
        if (i%5==0) fprintf(fout,"\n%5d:  ", i);
	fprintf(fout, "%10.*f  ", D, grib_data[i]);
	}
    fprintf(fout, "\n");

/*
*	.....     <<< END MARKER  >>>     .....
*       .....  RETAIN THE REMAINING CODE  .....
*/

/*
* Done, now free up float array;
*/
   if (grib_data!=NULL) { free(grib_data); grib_data = NULL; }
   
 } /* FOR-loop */
 
 nReturn = 0;		/* all done, set to no error */


bail_out:
/*
* IF (there is an error message) , print it out
*/
   if (errmsg[0] != '\0') 
 	fprintf(stderr,"\n***ERROR: %s\n", errmsg); 
/*
*
* CLEAN up before leaving
*       !free data array if defined
*	!close up output files if still opened
*       !free up Grib Header struct
*/
   if (grib_data!=NULL) free(grib_data);
   if (fout != NULL) fclose (fout);
   free_gribhdr (&gh1); 
   if (!errmsg[0] && InFile!=NULL && decoded_cnt==0) 
   fprintf(stdout, "No GRIB messages decoded from '%s'\n", InFile);
   exit ( nReturn);
}
