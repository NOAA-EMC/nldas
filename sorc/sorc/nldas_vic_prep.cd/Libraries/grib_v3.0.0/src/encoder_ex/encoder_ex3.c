/*.............................................................................
File:  encoder_ex3.c

Example program to encode a list of GRIB messages.

Usage:   encoder_ex3 model_name geom_name listfn
where
  model_name is the name of the Model
  geom_name is name of the geometry.  The program expects to load an
 	geometry file by this name with a .GEOM extension.
  listfn is  a text file containing the names of files to process, 1 per line.


One GRIB Edition-1 message is created per IEEE file (names stored in List file)
using the grib_enc library function and written to file using the default 
file name.  In this example, each IEEE filename in the Listfile has the format 
'ieee.parmname.lvltype.l1.l2.yyyymmddhh.tau' and the program will extract
the unique information of each field such as parmameter name, level type, 
height, date time group and forecast time and store in DATA_INPUT structure.

INPUT FILES:  
  $GRIB_ENV/tables/neons2grib.2.1   For mapping the Parameter & Level in the
				    IEEE filename & getting unit conversion info
  $GRIB_ENV/config/encoder.config   Encoder Configuration file
  $GRIB_ENV/data/'$GEOM'.geom       Geometry information file
  $list_fn			    Contains names of IEEE files, one per line.
  ieee.parmname.lvltype.l1.l2.yyyymmddhh.tau
  				    must be located in $GRIB_ENV/data/. Contains
				    float data that is LOCAL units, must be 
				    converted to GRIB unit prior to encoding.
OUTPUT FILES:  
  One Edition 1 GRIB message per entry in List file.
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
   14Aug97 first version;
   23Feb98 set forecast time unit to Hour;
   02Jul98 fix comments;
   23Sep98 Apply the Parameter's scale factor and reference to the data loaded
	   from IEEE file to convert it to GRIB units.
.............................................................................*/

#include <stdio.h>
#include <stdlib.h>
#include "gribfuncs.h"	  /* Library library function prototypes */
#include "grib_lookup.h"  /* Conversion Info */
#include "dprints.h"  	  /* Debug printing info*/

#define  CONFIG_FN        "config/encoder.config"

extern PARM_DEFN   db_parm_tbl[];  /* holds parm conversion info */
extern LVL_DEFN    db_lvl_tbl[];   /* holds level conversion info */
extern MODEL_DEFN  db_mdl_tbl[];   /* holds model conversion info */
extern GEOM_DEFN   db_geom_tbl[];  /* holds Geom conversion info */

main (int argc, char *argv[])
{
FILE	    *flist = NULL;  /* List file pointer */
GRIB_HDR    *gh= 0;         /* will hold encoded msg & its info */
DATA_INPUT  data_input;     /* input header structure for Encoder */
GEOM_IN     geom_in;        /* geometry description for Encoder */
USER_INPUT  user_input;     /* user input from input.dat for Encoder */
char	    *geom_name;     /* cmd line input */
char	    *model_type;     /* cmd line input */
char	    *list_fn;	    /* name of List filename */
char	    *grib_env;	    /* working variable */
char	    line[100],temp[100],dummy[100];  /* working vars */
char	    *p1,*p2;
char        errmsg[2000];   /* buffer to hold error message */
char	    config_fn[200];/* name of input filename */
char	    lookup_fn[100]; /* name of Conversion file */
char	    geom_fn[100];   /* name of Geometry file */
char	    ieee_fn[200];  /* name of input filename */
char	    out_fn[50];	    /* name of output filename */
float	    *flt_arr=NULL;  /* array to hold float data, init to null */
float	    parmscl, parmref; /* to change unit of the parameter */
float	    lvlscl, lvlref; /* to change unit of the level */
float	    min, max; 	    /* min and max value of data */
int	    id,cnt,i;	    /* working variables */
int	    quit;	    /* set to quit */
int	    stat;	    /* working var */
unsigned long uldtg;

    /* Get Command Line  Arguments */
    if (argc != 4) 
    	{ fprintf(stderr,"Usage= %s Modeltype Geomname List_fn\n\n", argv[0]);
	  exit(1); }

    fprintf(stdout,"\nStarting %s\n", argv[0]);
    model_type = argv[1]; geom_name= argv[2];  list_fn= argv[3];

    /* Check Environment Var */
    grib_env = getenv ("GRIB_ENV");
    if (grib_env == NULL || *grib_env == '\0') {
	fprintf(stderr,"Environment variable GRIB_ENV not defined\n");
	exit(1);
	}
    fprintf(stdout,"GRIB_ENV = %s\n", grib_env);

    fprintf(stdout,"Clear out the 3 Encoder structures \n");
    init_enc_struct (&data_input, &geom_in, &user_input);

    /* Load Encoder Configuration */
    sprintf (config_fn, "%s/%s", grib_env, (char *)CONFIG_FN);
    fprintf(stdout,"Loading Config file= %s\n", config_fn);
    if (ld_enc_config (config_fn, &user_input, errmsg) != 0 ) {
	fprintf(stderr,"Fatal error= %s\n", errmsg);
	exit(1);
	}

    /* 	Build Lookup filename to be used;
  	Parameter and Level name are extracted from the IEEE filename then
	used to map to the Parameter ID and Level ID number using info from
 	this lookup file.  The Unit conversion info is also loaded.
     */
    sprintf (lookup_fn, "%s/tables/neons2grib.%d.%d", grib_env,
	user_input.usParm_tbl, user_input.usSub_tbl);
    fprintf(stdout,"Loading Lookup file= %s\n", lookup_fn);

    /*  Go Load Lookup file into the following structure arrays:
	   PARM_DEFN   db_parm_tbl[NPARM]
		attribs ->usParm_id, ->usParm_sub, ->fScale, ->fOffset, ->sDSF;
	   LVL_DEFN    db_lvl_tbl[NLEV];
		attribs ->usLevel_id, ->db_name, ->fScale, ->fOffset;
	   MODEL_DEFN  db_mdl_tbl[NMODEL];
		attribs ->usModel_id, ->db_name;
	   GEOM_DEFN   db_geom_tbl[NCTRS];
		attribs ->usGeom_id, ->db_name;
     */
    if (ld_enc_lookup (lookup_fn, errmsg)) {
	fprintf(stderr,"Fatal error= %s\n", errmsg);
	exit(1);
	}
    
    /* Check input Geom */
    for (id=0; id < NGEOM; id++) 
	if (!strcmp (db_geom_tbl[id].db_name, geom_name)) 
	{   data_input.usGrid_id = id; break;  }	/* found it */

    if (id==NGEOM) { 
	fprintf(stderr,"Error:  invalid Geom %s\n", geom_name); 
	exit(1); 
	}

    /* Check input Model */
    for (id=0; id < NMODEL; id++) 
	if (!strcmp (db_mdl_tbl[id].db_name, model_type)) 
	{   data_input.usProc_id = id; break;  }	/* found it */

    if (id==NMODEL) {
	fprintf(stderr,"Error:  invalid Model %s\n", model_type); 
	exit(1); 
	}

    data_input.usFcst_id= 1;    /* set Fcst Time Unit to Hours */

    /* build the name of Geometry filename */
    sprintf (geom_fn, "%s/data/%s.geom", grib_env, geom_name);
    fprintf(stdout,"call ld_enc_geomfile (%s)\n", geom_fn);

    /* Load into GEOM_IN struct */
    if (ld_enc_geomfile (geom_fn, &geom_in, errmsg) != 0 ) {
	fprintf(stderr,"Fatal error= %s\n", errmsg);
	exit(1);
	}

    fprintf(stdout,"Malloc float array %d by %d\n", geom_in.nx, geom_in.ny);
    if (!(flt_arr = (float *) malloc(geom_in.nx * geom_in.ny *sizeof(float)))) {
	fprintf(stderr,"Failed to malloc Float array\n");
	exit(1);
  	}

    fprintf(stdout,"Prepare List file '%s' for reading\n", list_fn); 
    if ( !(flist = fopen (list_fn, "r")) ) {
	  fprintf(stderr,"Unable to open Listfile %s\n", list_fn);
   	  if (flt_arr != NULL) free(flt_arr);  /* release float array */
   	  exit(1);
	}

     fprintf(stdout,"Allocate storage and initialize GRIB_HDR structure\n");
     if (init_gribhdr (&gh, errmsg))  {
          fprintf(stderr,"Abort; error=%s\n", errmsg);
	  if (flist) fclose (flist);
   	  if (flt_arr != NULL) free(flt_arr);  /* release float array */
   	  exit(1);
   	}
    
    /* Loop: create a Grib msg per entry in List file */
    for (stat=0;  !feof(flist) && !ferror(flist);  ) 
    {
       memset ((void*)line, '\0', sizeof line);
       memset ((void*)dummy, '\0', sizeof dummy);
       memset ((void*)temp, '\0', sizeof temp);
       if (fgets(line, 100, flist) == NULL) break;
       if (sscanf (line, "%s%s", temp, dummy) != 1) {
	  fprintf(stderr,"Invalid List_fn entry:  [%s %s...]\n", temp,dummy);
	  continue; }

        fprintf(stdout,"\n[ %s ]=\n", temp);

	/* 
	    Each message to be encode requires the actual float data and
	    some header information on the field.  In this example, we
	    extract these info from the filename.
	    
	    Extract info from filename:
		'ieee.parmname.lvltype.l1.l2.yyyymmddhh.tau'
	*/
	for (p1=temp, quit=i=0; !quit && i < 7 ; i++, p1=p2+1) 
	{
	   if ((p2=(char *)strchr(p1, '.')) == NULL) p2=temp+strlen(temp);
	   /* Extract */
	   strncpy (dummy, p1, p2-p1);	dummy [ p2-p1] = '\0';

	   switch (i) {
	      case 0:  /* 'ieee.' part */
		if (strcmp (dummy, "ieee")) {
		  fprintf(stderr,"No 'ieee.', drop [%s]\n",temp); quit=1; } 
		break;

	      case 1:  /* PARM part, 
		go fill data_input's usParmid, usParmsub, nDec_sc_fctr;
		Save the scale factor and reference for converting the
		data into GRIB units later */
		if (map_parm (dummy, &data_input, &parmscl, &parmref, errmsg)) {
		  fprintf(stderr,"%s, drop [%s]\n", errmsg,temp); quit=1; } 
		break;
		
	      case 2:  /* LEVEL part, 
		go fill data_INPUT'S usLevelid, SCALE up lvl1 & lvl2 
		which currently have not been read in yet.
		Save the scale factor and reference for scaling the
		level value later. */
		if (map_lvl (dummy, &data_input, &lvlscl, &lvlref, errmsg)) {
		  fprintf(stderr,"%s, drop [%s]\n", errmsg,temp); quit=1; } 
		break;
		
	      case 3:  /* LEVEL_1 part, scale it */
		if (strcspn (dummy, "0123456789") != 0) {
		  fprintf(stderr,"Bad level_1, drop [%s]\n", temp); quit=1; }
		else data_input.nLvl_1 = lvlscl*atoi(dummy) + lvlref; 
		break;
		
	      case 4:  /* LEVEL_2 part, scale it */
		if (strcspn (dummy, "0123456789") != 0) {
		  fprintf(stderr,"Bad level_2, drop [%s]\n", temp); quit=1; }
		else data_input.nLvl_2 = lvlscl*atoi(dummy) + lvlref; 
		break;

	      case 5:  /* DTG part */
		if (strlen(dummy)!=10 || strcspn (dummy,"0123456789")!= 0)
		  { fprintf(stderr,"Bad DTG, drop [%s]\n", temp); quit=1; }
		else {
			uldtg = (unsigned long)atol (dummy);
			data_input.nYear = uldtg / 1000000;
			data_input.nMonth = (uldtg / 10000) % 100;
			data_input.nDay = (uldtg / 100) % 100;
			data_input.nHour = (uldtg % 100); 
			if (*(p2+1) == '\0') { fprintf(stderr,
		  	"Missing Forecast Period, drop [%s]\n", temp); quit=1; }
		      }
		break;

	      case 6:
	        /* Forecast Period part */
		if (strcspn (dummy, "0123456789") != 0) {
		  fprintf(stderr,"Bad Forecast Period, drop [%s]\n", temp); 
		  quit=1; }
		else 
		{  data_input.usFcst_per1 = (unsigned short)atoi(dummy); 
		   if (*(p2+1) != '\0') 
		      { fprintf(stderr,
		       "Invalid info AFTER Fcst_per, drop [%s]\n", temp); 
			quit=1;} 
		}
		break;
	   } /* Switch */
	} /* For */

     if (quit) continue;
     if (i < 6) {
	fprintf(stderr,
	"Proper entry= 'ieee.$parm.$lvl.$lvl1.$lvl2.yyyymmddhh.tau', drop [%s]\n", temp);
	continue; 
	}
	
	/* Go load IEEE file into Float array,  skip if failed */
       sprintf (ieee_fn, "%s/data/%s", grib_env, temp);
       fprintf(stdout,"call ld_enc_ieeeff (%s)\n", ieee_fn);
       if (ld_enc_ieeeff (ieee_fn, flt_arr, geom_in.nx * geom_in.ny, errmsg))
   	   {
           fprintf(stderr,"Failed to load IEEE file %s, ermsg=%s\n", 
	   ieee_fn, errmsg);
           continue;
           }
       
	/* Ensure data is in GRIB units, convert if needed */
	if (parmref != 0.0 || parmscl != 1.0) {
	   fprintf(stdout,"applying Scale Fctor= %f and Reference=%f to data\n",
	   parmscl, parmref);
	   for (i=0; i <  geom_in.nx * geom_in.ny; i++) 
		flt_arr[i] = flt_arr[i] * parmscl + parmref;
	  }
	else fprintf(stdout,
	   "No need to apply Scale Fctor= %f and Reference=%f to data\n",
	   parmscl, parmref);

   	/* print out a sample of the data :
           for (i=0; i < 100; i++) 
   	    if (i % 5 == 4)
   		 fprintf(stdout,"%.*lf\n",
   		  (data_input.nDec_sc_fctr>0? data_input.nDec_sc_fctr: 0), 
   		  flt_arr[i]); 
   	    else fprintf(stdout,"%.*lf ",
   		  (data_input.nDec_sc_fctr>0? data_input.nDec_sc_fctr: 0), 
   		  flt_arr[i]); 
	*/
   		  
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

	/* PRINT OUT CONTENT OF THE 3 STRUCTS */
      user_input.chCase_id='3';  
      fprintf(stdout,
	"Encode message with Forced CASE ID= %c\n",user_input.chCase_id);
	 	  P_CHAR (user_input.chCase_id);
                  P_USHORT (user_input.usParm_tbl);
		  P_USHORT (user_input.usSub_tbl );
                  P_USHORT (user_input.usCenter_id );
                  P_USHORT (user_input.usCenter_sub );
                  P_USHORT (user_input.usTrack_num );
                  P_USHORT (user_input.usGds_bms_id );
                  P_USHORT (user_input.usBit_pack_num );

		  P_USHORT (data_input.usProc_id); 
                  P_USHORT (data_input.usGrid_id); 
                  P_USHORT (data_input.usParm_id); 
                  P_USHORT (data_input.usParm_sub_id); 
                  P_USHORT (data_input.usLevel_id); 
                  P_INT (data_input.nLvl_1); 
                  P_INT (data_input.nLvl_2); 
                  P_INT (data_input.nYear); 
                  P_INT (data_input.nMonth); 
                  P_INT (data_input.nDay); 
                  P_INT (data_input.nHour); 
                  P_USHORT (data_input.usFcst_per1); 
                  P_INT (data_input.nDec_sc_fctr);

		  P_STRING (geom_in.prjn_name); 
                  P_INT (geom_in.nx); 
                  P_INT (geom_in.ny); 
                  P_DOUBLE (geom_in.x_int_dis); 
                  P_DOUBLE (geom_in.y_int_dis); 
                  P_DOUBLE (geom_in.parm_1); 
                  P_DOUBLE (geom_in.parm_2); 
                  P_DOUBLE (geom_in.parm_3); 
                  P_DOUBLE (geom_in.first_lat); 
                  P_DOUBLE (geom_in.first_lon); 
                  P_DOUBLE (geom_in.last_lat); 
                  P_DOUBLE (geom_in.last_lon); 
                  P_USHORT (geom_in.scan); 
                  P_USHORT (geom_in.usRes_flag ); 
                  P_USHORT (geom_in.usRes_flag ); 
   	 
      /* go encode message now
	 with DATA_INPUT struct (unique info on this msg)
	 and USER_INPUT struct (encoder configuration data)
	 and GEOM_IN struct (geometry info)
	 and float data array (actual floating point data)
	 and errmsg buffer (filled by grib_enc when error occurs)
      */
      if (grib_enc (data_input, user_input, geom_in, flt_arr, gh, errmsg)!= 0) {
   	  fprintf(stderr,"Abort; error=%s\n",errmsg);
          stat=1; break;
        }
   
      /* If desired, display content of GRIB_HDR struct 
      	display_gribhdr (gh);
      */
   
      /* form a filename for the output file which reflects
         the type of message this is
      */
      make_default_grbfn (data_input, user_input, out_fn);
   
      /* Go save message in GRIB_HDR struct out to Output file
      */
      fprintf(stdout,"ENCODED FILE= '%s'\n"\
      "output file length should be %ld bytes\n\n", out_fn,gh->msg_length);
   
      if ( gribhdr2file (gh, out_fn, errmsg) != 0) {
           fprintf(stderr,"Abort; error= %s\n", errmsg);
	   stat=1; break;
           }
   } /* loop for each Ieee fn in List_fn */

   if (flist != NULL) fclose (flist);
   if (flt_arr != NULL) free(flt_arr);  /* release float array */
   free_gribhdr (&gh);			/* release grib header storage */
   exit(stat);
}
