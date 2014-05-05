/* 08/19/97 by Alice T. Nakajima, SAIC, MRY
*********************************************************************
* FILENAME:  getgribieee.c
*   Used to decode a list of GRIB Edition 1 messages using the GRIB Library 
*   and create an 'IEEE' file per message found.
* INPUT =
*   1)  a List file ($dir/$listfn, whith the $dir being optional) which 
*       contains names of all of the GRIB input files to decode.  
*       Each line entry has format '$dir/$inputfn' with $dir being optional.
*   2)  the files whose names appear in the List file mentioned.
* OUTPUT= '$InFile.IEEE', 1 file per GRIB message found.
*
* NOTES:  
*   -Assume there will only be 1 GRIB message per input file. 
***********************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "dprints.h"		/* for all debug printing */
#include "gribfuncs.h"		/* all GRIB func prototypes */

#define  FILL_VALUE -9999.00	/* Value for missing datapts  */

char        errmsg[2000];       /* leave it big! used by Grib Library */
GRIB_HDR    *gh= NULL;	  	/* grib header block, filled by Encoder*/

/*
************************************************************************
* A.  FUNCTION:  main
************************************************************************
*/
void main (int argc, char **argv) {
int 	        i;
char		*grib_env;	/* pts to Environment Variable */
char 	        InFile[200];	/* current input GRIB filename */
char		ieee_fn[201];	/* current output Ieee filename */
int             nReturn=0;      /* return status from grib_dec           */
int		S_stat;	        /* status from grib seek function */
int		n;
int		Read_Index=0;   /* Zero so Seek will search until get a Msg */
unsigned long   msg_length = 0; /* the total length of the message */
long            offset;   /* offset within the multiple GRIB file, indicating 
		         end of current msg relative to the top of the file*/
float 		*grib_data=NULL;        /* array of decoded data values */

GRIB_HDR 	 *gh1;		/* holds info on msg returned from Seek */
BMS_INPUT	 bms;		/* input structure for bitmap section   */
PDS_INPUT        pds;           /* product definition section */
grid_desc_sec    gds;           /* grid description section */
BDS_HEAD_INPUT   bds_head;      /* input structure for bds section */
FILE		*fpo, *flist;

/*
* A.1       CLEAR out name of input file variable
*/
   InFile[0]='\0';
        
/*
* A.2       PARSE command line arguments !exit if error
*           Expecting List filename as the argument
*/
   if (argc!=2 ||    (argv[1][0]=='\0' || argv[1][0]=='-') ) {
	fprintf(stderr, "Usage:  %s  List_fn \n",argv[0]); exit(0); }

/*
* A.3       OPEN list file for reading !exit if error
*/
   if ((flist = fopen(argv[1], "r")) == NULL)  {
	fprintf(stderr,
	"Unable to open List file %s;\nAbort Program\n", argv[1]); exit(1); }
   fprintf(stdout,"List file is '%s'\n", argv[1]);

/*
*
* A.4       FUNCTION init_gribhdr   !make storage for grib header struct
*           EXIT if fails;
*/
    if (init_gribhdr (&gh1, errmsg))
    {   fprintf(stderr,"%s:  %s\nAbort Program\n", argv[0], errmsg);
	fclose(flist);exit(1);
    }

/*
*
* A.5       LOOP while not end of List file yet & no errors yet
*/
while (!ferror(flist) && !feof(flist)) 
{
/*
* A.5.1       READ next entry - name of GRIB input file
*             QUIT if fails
*/
  if ((n=fscanf (flist, "%s", InFile)) != 1)  break; 
  fprintf(stdout,"\nNext Infile= %s\n", InFile);

/*
* A.5.2       KEEP doing:
*                FUNCTION grib_seek   !find next message in input file
*                !Proceed only if no errors
*/
  for (offset=0L, errmsg[0]='\0'
	; 	! (S_stat= grib_seek(InFile, &offset, Read_Index, gh1, errmsg))
	; offset += (msg_length+ 1L), errmsg[0]='\0')
  {
/*
* A.5.2.1       IF (message has length > 0) THEN
*/
  if ((msg_length = gh1->msg_length) > 0) 
   {
      fprintf(stdout, "Decoding message found at %ld bytes ...\n",offset);

/*
* A.5.2.1.1        FUNCTION init_dec_struct   !initialize decoder structures
*/
      init_dec_struct(&pds,&gds,&bms,&bds_head);

/*
* A.5.2.1.2        FUNCTION grib_dec      !perform GRIB decoding
*                  !results in internal structures holding 
*                  !PDS, GDS , BDS, BMS, and new float data array
*                  IF (error decoding message)
*                      PRINT error
*                      EXIT
*                  ENDIF
*/
      grib_data= (float *) NULL;

      if (nReturn = grib_dec((char*)gh1->entire_msg,
                &pds,&gds,&bds_head,&bms,&grib_data, errmsg))
      {
	 fprintf(stderr, "%s %s: %s\n", argv[0], InFile, errmsg); 
   	 display_gribhdr(gh1); 
         break;
      }

/*
* A.5.2.1.3        IF (applying Bitmap Section to Data AND 
*                    Bitmap Section is present) 
*                  THEN
*                      FUNCTION apply_bitmap !apply BMS to 'grib_data'
*                      EXIT if  error applying bitmap
*                  ENDIF
*/
      if (bms.uslength>0)
	 if (nReturn=apply_bitmap (&bms, &grib_data, FILL_VALUE, &bds_head,
	     errmsg)) 
	 { fprintf(stderr, "%s %s: %s\n", argv[0], InFile, errmsg); break; }

/*
* A.5.2.1.4         !Optional: FUNCTION prt_inp_struct  !show internal structs
*
*
#        fprintf(stdout,
# 	"Main: Decoder Completed;\n\nShow Decoded Message info now=\n");
#        prt_inp_struct(&pds, &gds, &bms, &bds_head, &grib_data); 
*/


/*
* A.5.2.1.5   CREATE  output file to hold the Data, named '$InFile.IEEE'
*             QUIT if fails
*             Assumption:  only 1 GRIB message per input file
*/
	sprintf (ieee_fn, "%s.IEEE",
		(strrchr(InFile,'/') ? strrchr(InFile,'/') + 1 : InFile));

	   if ( ! (fpo = fopen(ieee_fn, "wb"))) {
		    fprintf(stderr,
		    "Cannot open  %s for writing\nAbort Program;\n", ieee_fn);
		    fclose(flist); exit (1); 
	     }

	   if ( fwrite(grib_data, sizeof(float), bds_head.ulGrid_size, fpo)
	        != bds_head.ulGrid_size) 
	     {
		fprintf(stderr, 
		"Failed to write out %d float elements to %s\nAbort Program;\n",
		bds_head.ulGrid_size,ieee_fn);
		fclose(fpo); fclose(flist);exit(1);
	     }

	   fprintf(stdout,"Ieee output= %s\n", ieee_fn);
	   fclose(fpo);

/*
* A.5.2.1.6       IF (float buffer have data in it) THEN
*                  FREE up its storage;
*               ENDIF
*/
      if (grib_data!=NULL) { free(grib_data); grib_data = NULL; }

/*
* A.5.2.1       ENDIF 	!got a message
*/
   }

/*
* A.5.2       ENDLOOP
*/
  } /* FOR */


/*
* A.5.3     IF  last call to seek grib returned 'corrupted len' error THEN
*               print as much of the grib hdr struct as possible
*           ENDIF
*/
   if (S_stat != 0 && errmsg[0]!= '\0') 
	fprintf(stderr, "%s %s: %s\n", argv[0], InFile, errmsg); 

/*
* A.5       END LOOP
*/
}

/*
* A.6       HOUSE clean & exit !free storage, close files
*/
   fclose(flist);   
   if (grib_data!=NULL) { free(grib_data); grib_data = NULL; }
   free_gribhdr (&gh1);  /* alloced in ENcode() */
   exit(0);

/*
*
* END OF FUNCTION
*
*/
}
