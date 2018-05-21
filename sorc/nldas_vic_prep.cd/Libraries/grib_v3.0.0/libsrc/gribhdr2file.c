#include <stdio.h>
#include "dprints.h"		/* for dprints */
#include "gribfuncs.h"		/* prototypes */
/*
*
************************************************************************
* A.  FUNCTION  gribhdr2file
*       write out the Grib message stored in GRIB_HDR struct to external file;
*       if the 'shuffle' flag is set, write each individual section out, else
*       write 'entire_msg' all at once;
*               
*    INTERFACE:
*       int    gribhdr2file (gh, fn, errmsg)
*
*    ARGUMENTS (I=input, O=output, I&O=input and output):
*      (I) GRIB_HDR *gh    holds the GRIB message to be written out
*      (I) char *fn        name of file to write to (includes absolute path)
*      (O) char *errmsg    array returned empty unless error occurred;
*
*     RETURN CODE:
*     0>  no errors,  GRIB file successfully created;
*     1>  error; errmsg is filled;
************************************************************************
*/
#if PROTOTYPE_NEEDED
int    gribhdr2file ( GRIB_HDR *gh, char *fn, char *errmsg)
#else
int    gribhdr2file ( gh, fn, errmsg)
			GRIB_HDR *gh; 
			char *fn; 
			char *errmsg;
#endif
{
/*
*
* A.0     DEFAULT to error status of 1
*/
FILE *f1;
char *func= "gribhdr2file";
int  stat=1;		

/*
*
* A.1     IF (entire msg array is null or msg length is 0)
*         THEN 
*            RETURN error stat !errmsg filled
*         ENDIF
*/
  DPRINT1("Entering %s\n", func);
  if (gh->entire_msg == NULL || gh->msg_length <= 0) {
	DPRINT1 ("%s: GRIB_HDR message buffer is null, OR msg_length=0\n",func);
	sprintf(errmsg,"%s: GRIB_HDR message buffer is null, OR msg_length=0\n",
	func);
	goto BYE;
	}	

/*
*
* A.2     IF (in Shuffle mode)
*         THEN 
*            IF (length of EDS/PDS/BDS/EDS is 0) THEN
*               RETURN error stat !errmsg filled
*            ENDIF
*         ENDIF
*/
  if (gh->shuffled) {
	if (!gh->ids_len|| !gh->pds_len || !gh->bds_len|| !gh->eds_len) {
	   DPRINT1("%s:  Shuffle mode: Zero length encountered, quit\n", func);
	   sprintf(errmsg,
	   "%s:  Shuffle mode: Zero length encountered, quit\n", func);
	   goto BYE; }
	DPRINT1 ("%s:   this mesg is in shuffled mode;\n", func);
	}

/*
*
* A.3     OPEN file for output;
*         IF (failed) THEN
*               RETURN error stat !errmsg filled
*         ENDIF
*/
  if (! (f1 = fopen (fn, "wb"))) {
	DPRINT2("%s: Unable to open %s\n",func,fn); 
	sprintf(errmsg,"%s: Unable to open %s\n",func,fn); 
	goto BYE; 
	}
  
  DPRINT1 ("creating %s\n", fn);
/*
*
* A.4     IF (in shuffled mode) 
* A.4.a   THEN
*/
  if (gh->shuffled) {
/*
* A.4.a.1    IF (fails to write IDS OR fails to write PDS OR
*                (GDS exists AND fails to write GDS) OR 
*                (BMS exists AND fails to write BMS) OR 
*                fails to write BDS or fails to write EDS)
*            THEN
*               RETURN error stat !errmsg filled
*            ENDIF
*/
     if (fwrite (gh->ids_ptr , gh->ids_len, 1, f1) != 1)
        {
          DPRINT1 ("%s:  failed to Fwrite IDS to file\n", func);
          sprintf(errmsg,"%s:  failed to Fwrite IDS to file\n", func);
          goto BYE;
        }
     if (fwrite (gh->pds_ptr , gh->pds_len, 1, f1) != 1)
        {
          DPRINT1 ("%s:  failed to Fwrite PDS to file\n", func);
          sprintf(errmsg,"%s:  failed to Fwrite PDS to file\n", func);
          goto BYE;
        }
     if (gh->gds_len)
     if (fwrite (gh->gds_ptr , gh->gds_len, 1, f1) != 1)
        {
          DPRINT1 ("%s:  failed to Fwrite GDS to file\n", func);
          sprintf(errmsg,"%s:  failed to Fwrite GDS to file\n", func);
          goto BYE;
        }
     if (gh->bms_len)
     if (fwrite (gh->bms_ptr , gh->bms_len, 1, f1) != 1)
        {
          DPRINT1 ("%s:  failed to Fwrite BMS to file\n", func);
          sprintf(errmsg,"%s:  failed to Fwrite BMS to file\n", func);
          goto BYE;
        }
     if (fwrite (gh->bds_ptr , gh->bds_len, 1, f1) != 1)
        {
          DPRINT1 ("%s:  failed to Fwrite BDS to file\n", func);
          sprintf(errmsg,"%s:  failed to Fwrite BDS to file\n", func);
          goto BYE;
        }
     if (fwrite (gh->eds_ptr , gh->eds_len, 1, f1) != 1)
        {
          DPRINT1 ("%s:  failed to Fwrite EDS to file\n", func);
          sprintf(errmsg,"%s:  failed to Fwrite EDS to file\n", func);
          goto BYE;
        }
     DPRINT0 ("ALL Sections to written to file successfully\n");
     }
/*
* A.4.b   ELSE
*/
   else {
        DPRINT0 ("Writing gh->entire_msg (non-shuffled)\n");
/*
* A.4.b.1    IF (fails to write msg_length byte straight from Entire_msg)
*            THEN
*                RETURN error stat !errmsg filled
*            ENDIF
*/
        if (fwrite (gh->entire_msg, gh->msg_length, 1, f1) != 1) {
	    DPRINT1( "%s:  failed to write GH's entire Msg to file\n",func);
	    sprintf(errmsg,
	    "%s:  failed to write GH's entire Msg to file\n",func);
	    goto BYE; 
	   }
  	DPRINT0 ("write GH's entire_msg to file successful\n");
/*
* A.4     ENDIF
*/
        }
/*
*
* A.5     DONE, set status to 0  !no errors
*/
  stat = 0;

BYE:
/*
*
* A.6     CLOSE output file if not already closed
*/
  if (f1!=NULL) fclose(f1);
/*
*
* A.7     RETURN with stat
*/
  DPRINT2 ("Leaving %s, stat=%d;\n", func, stat);
  return stat;
/*
*
* END OF FUNCTION
*/
}
