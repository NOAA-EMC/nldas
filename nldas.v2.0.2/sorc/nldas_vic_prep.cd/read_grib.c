#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "preproc.h"
#include "gribfuncs.h"

#define DEBUG_GRIB	0
/*  
     
   Read a grib record from file.										

*/

int read_grib(char *fname, int id, int proc, float **grib_data, int reset)
{

  int      nReturn=0;        /* return status from Grib_dec   */
  char     errmsg[BUFSIZ+1]; /* required, error msgs back from Grib library */
  static long offset=0L;   	 /* byte offset from beginning of file, for inputfiles
				with multiple GRIB messages; */
  int      Rd_Indexfile=0;   /* Zero to decode all msgs in input file  */
  int	   decoded_cnt;      /* number of messages decoded successfully */

  BMS_INPUT	   bms;	      /* input structure for bitmap section   */
  PDS_INPUT        pds;       /* input structure for product definition section */
  grid_desc_sec    gds;       /* input structure for grid description section */
  BDS_HEAD_INPUT   bds_head;  /* input structure for bds section */
  GRIB_HDR         *gh;      /* holds Msg and its info, put there by grib Seek */

  int i, match=0, tries=0;

  /* initialise the error message to a null string. 
     Otherwise the test of the string contents 
     below may give incorrect results */
  strcpy(errmsg, "\0");
  /* MAKE storage for Grib Header !exit on error */
  if((nReturn= init_gribhdr (&gh, errmsg) ))  goto bail_out; 

  if(reset) offset=0L;

  /* read thru the grib data file */
//  for(offset=0L, decoded_cnt=0; nReturn == 0;  offset += gh->msg_length) {
  for(decoded_cnt=0; nReturn == 0;  offset += gh->msg_length) {

#if DEBUG_GRIB
    fprintf(stderr, ">>> Grib record %d\n", decoded_cnt);
#endif

    /* move to the next record */
    if( (nReturn = grib_seek(fname, &offset, Rd_Indexfile, gh, errmsg))) {	
#if DEBUG_GRIB
      fprintf(stdout,"Grib_seek returned non zero stat (%d)\n", nReturn);
#endif
#if 0
      if(nReturn == 2) break;	/* End of file error */
#else
      if(nReturn == 2) {	/* End of file error */
        /* rewind file */
        if(tries < 1) {
          nReturn = 0;
          offset = 0L;
          tries++;
        }
        /* we have already rewound the file */
        else {
          goto bail_out;
        }
      }
#endif
      else goto bail_out;		/* abort if other error */
    }
	
    /* NO errors but got a Warning msg from seek */
    if (errmsg[0] != '\0'){ 
#if DEBUG_GRIB
      fprintf(stderr, "Warning: %s; Skip Decoding...\n",errmsg);  
#endif
      errmsg[0]='\0'; 	       /* reset error to continue */
      gh->msg_length = 1L;        /* set to 1 to bump offset up */
      continue;
    }

    if(gh->msg_length < 0) {
      fprintf(stderr, "Error:  message returned had bad length (%ld)\n", gh->msg_length);
      goto bail_out;
    }
    else if(gh->msg_length == 0) {  
#if DEBUG_GRIB
      fprintf(stderr,"msg_lenth is Zero, set offset to 1\n");
#endif
      gh->msg_length = 1L;  /* set to 1 to bump offset up */
      continue; 
    }
	
    /* Clear out the Structures */
    init_dec_struct(&pds, &gds, &bms, &bds_head);
#if DEBUG_GRIB
    fprintf(stderr,"Decoding message found at %ld bytes ...\n",offset);
#endif
	
    /* decode the Grib record */
    if(*grib_data != NULL) { free(*grib_data); *grib_data = (float *)NULL; }  
    if((nReturn = grib_dec ((char *)gh->entire_msg, &pds, &gds, &bds_head, &bms, grib_data, errmsg))) goto bail_out;
#if DEBUG_GRIB
    fprintf(stderr, "grib id\t = %d\n", pds.usParm_id);
#endif
/*    fprintf(stderr, "grib id = %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", pds.uslength, pds.usEd_num, pds.usParm_tbl, pds.usCenter_id, pds.usProc_id, pds.usGrid_id, pds.usGds_bms_id, pds.usParm_id, pds.usLevel_id, pds.usLevel_octets, pds.usHeight1, pds.usHeight2, pds.usYear, pds.usMonth, pds.usDay, pds.usHour, pds.usMinute, pds.usFcst_unit_id, pds.usP1, pds.usP2, pds.usTime_range, pds.usTime_range_avg, pds.usTime_range_mis, pds.usCentury, pds.usCenter_sub, pds.sDec_sc_fctr, pds.ausZero[0], pds.usExt_flag, pds.usSecond, pds.usTrack_num, pds.usParm_sub, pds.usSub_tbl);*/

    decoded_cnt++;

    /* check that we have the correct variable. Note that the proc id is 84 in the retro but 84 or 89 in the realtime */
    if( id == pds.usParm_id && ( (proc == pds.usProc_id) || (proc == 84 && pds.usProc_id == 89) )) {

      match = 1;
	
      /* get the data */
      if(bms.uslength>0 && (nReturn = apply_bitmap(&bms, grib_data, FILL_VALUE, &bds_head, errmsg)))
        goto bail_out; 

#if 0	
      /* check the grib file and mask file are compatible*/
      if( (ldas_index.nr != gds.llg.usNj) || (ldas_index.nc != gds.llg.usNi) ){
        fprintf(stderr, "Incompatible grid dimensions [row cols]:\t %d %d\t%d %d\n", 
                ldas_index.nr, gds.llg.usNj, ldas_index.nc, gds.llg.usNi );
        goto bail_out;
      }
#endif

#if DEBUG_GRIB
      /* write out grib header info and 100 pieces of data */
      fprintf(stdout, "\n\t>>>  Content of GRIB Msg #%d =\n", decoded_cnt);
      prt_inp_struct(&pds, &gds, &bms, &bds_head, grib_data);
#endif

      break;

    }	

  } /* end of for (offset = 0L.. loop */

  /* clean up after every grib file */
  free_gribhdr (&gh);

  if(match)
    return(1);
  else
    return(0);

bail_out:
  /* IF (there is an error message) , print it out */
  if (errmsg[0] != '\0') 
    fprintf(stderr,"\n***ERROR: %s\n", errmsg); 

  return(0);


} /* read_grib */

