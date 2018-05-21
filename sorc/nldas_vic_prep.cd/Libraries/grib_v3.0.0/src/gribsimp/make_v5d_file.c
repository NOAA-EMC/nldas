#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "isdb.h"  /* DATE struct */
#include "gribfuncs.h"
#include "dprints.h"

#include "gsv5d.h"	  /* for all grads structs & links to v5d.h too */

#define  V5_ERRCODE  0	  /* return code of vis5d funcstions upon error */

extern char 	       *InFile;    /* name of input file */
extern BMS_INPUT	bms;	   /* input structure for bitmap section */
extern PDS_INPUT        pds;       /* product definition section */
extern grid_desc_sec    gds;       /* grid description section */
extern BDS_HEAD_INPUT   bds_head;  /* input structure for bds section */
extern GRIB_HDR 	*gh1;      /* already defined by now */
/*
*======================================================================
* A. FUNCTION:  get_htarray_index
*     Locate the cell with 'val' in the HtArray[numhts] and return its
*     array index.
*
*  INTERFACE:   get_htarray_index (val, HtArray, numhts)
*
*  ARGUMENTS (I)=Input, (O)=Output
*      (I) unsigned short val; 	     value to match 1 of array's element
*      (I) unsigned short *HtArray;  array to search in
*      (I) int numhts;		     num of elements in the array
*
*  RETURNS:
*    -1:  error, did not find a match. 
*   else: an index in the range of 0 to numhts-1;
*======================================================================
*/
#if PROTOTYPE_NEEDED
int get_htarray_index (	unsigned short val, 
			unsigned short *HtArray, 
			int numhts)
#else
int 	get_htarray_index (val, HtArray, numhts)

unsigned short val; 		/* value to match 1 of array's element */
unsigned short *HtArray; 	/* array to search in */
int numhts;			/* num of elements in the array */
#endif
{ int x;
  for (x=0;x < numhts;x++) if (HtArray[x]== val) return(x); 
  return(-1);
}

/*
*======================================================================
* B.  FUNCTION:  make_v5d_file
*     Create a vis5d output file from the Message list, the Variable
*     List and the Level structure.
* 
*  INTERFACE:  make_v5d_file  (outfile, NumVars, v5d_msgs_head,
*              v5d_parm_head, v5d_lvl_head, v5d_info, errmsg)
*  
*  ARGUMENTS (I)=Input, (O)=Output
*   (I&O)  char    *outfile;        array for output filename;  Will be
*                                    assigned a filename if Null;
*   (I)  int       NumVars;          number of Variables 
*   (I)  V5_REC    **v5d_msgs_head;  head of V5d Messages list 
*   (I)  V5D_PARM  **v5d_parm_head;  head of V5d Param list 
*   (I)  V5D_LVL   **v5d_lvl_head;   V5d Level info 
*   (I)  V5_INFO   v5d_info;         common info for all messages
* (I&O)  char      *errmsg;          buffer for error message 
*
* RETURN CODE:
*   0:  success, output file has name '$lvlid.v5d'
*   1:  failed, errmsg filled;
*======================================================================
*/
#if PROTOTYPE_NEEDED
int	make_v5d_file ( char		*outfile,
			int 		NumVars, 
			V5_REC   	**v5d_msgs_head, 
			V5D_PARM  	**v5d_parm_head, 
			V5D_LVL   	**v5d_lvl_head,
			V5_INFO  	v5d_info, 
			char 		*errmsg)
#else
int	make_v5d_file  ( outfile, chosen_name, NumVars, 
			v5d_msgs_head, v5d_parm_head, v5d_lvl_head, 
			v5d_info, errmsg)
	char	  *outfile; 		/* output file name */
	int	  NumVars;		/* number of Variables */
	V5_REC    **v5d_msgs_head;  	/* head of V5d Messages list */
	V5D_PARM  **v5d_parm_head;  	/* head of V5d Param list */
	V5D_LVL   **v5d_lvl_head;   	/* head of V5d Level List */
	V5_INFO   v5d_info;		/* common info for Vis5D data */
	char	  *errmsg;		/* array for error mesg */
#endif
{
   char *func="make_v5d_file";

   const char VarName[MAXVARS][10]; /* mapped names of these variables */
   const char Abbrv[MAXVARS][10];   /* 'PpidSsid' names of these variables */

   int LowLev[MAXVARS]; /* index (within LVL->height[]) of Start Heigt per Var*/
   int HiLev[MAXVARS]; /* index (within LVL->height[]) of End Height per Var*/
   int Nht[MAXVARS];    /* Count of Heights per variable */
			/* includes ALL Height steps between LowLev[Var] */
			/* HiLev[var], even if some are not present 	*/
			/*   Nht[0] is #Heights in 1st variable, */
			/*   Nht[1] is #Heights in 2nd variable, etc*/
   int MaxNht;		/* the largest value in Nht[] */

   int NumTimes;                /* Nos. of time steps in ValidETime[] */
   double ValidEtimes[MAXTIMES];/* holds Unique Valid times */
   int TimeStamp[MAXTIMES];     /* real times for each time step */
   int DateStamp[MAXTIMES];     /* real dates for each time step */

   int Projection;              /* a projection number */
   float ProjArgs[100];         /* the projection parameters */

   int Vertical;                /* a vertical coord system number */
   float VertArgs[MAXLEVELS];   /* the vertical coord sys parameters */

   int Nr, Nc;			/* #rows & #cols (constant for all msgs) */
   int CompressMode;            /* number of bytes per grid */
   int  doy;			/* day of year */
   DATE date;			/* Date structure (yr, month, day) */
   double hour;			/* time of the day in Hours */
   double temp;			/* working var */
   float  *farr=NULL;		/* data from a single msg, created by Decoder*/
   float  *v5darr=NULL;		/* data from ALL msgs */
   float  *EmptyGrid=NULL;	/* entire grid holding MISSING value */
   V5D_LVL *Lvl;    		/* local variable holds content from Lvl LIst*/
   V5_REC   *Msg;        	/* temp ptr */
   V5D_PARM  *Parm;       	/* temp ptr */
   int itm,iht,ivar,ind,itmp; /* working var */
   int bytespergrid;		/* num of Bytes per [row * col] */
   int file_open=0;		/* marks a sessions going */
   int stat=1;			/* return code, default is error */

#ifndef V5D_ON
   fprintf(stdout,"%s: cannot proceed, V5D_ON not available\n",func);
   return(0);

#else
  /******************************************************************* 
			** NOTE **
   This code only gets included if -DV5D_ON was used during compilation 
   else go straight to 'bail_out';
  *******************************************************************/
   fprintf(stdout,"\nEntering %s:\n", func);
   memset ((void*)ValidEtimes, '\0', MAXTIMES*sizeof(int));
   memset ((void*)TimeStamp, '\0', MAXTIMES*sizeof(int));
   memset ((void*)DateStamp, '\0', MAXTIMES*sizeof(int));
   memset ((void*)Nht, '\0', MAXVARS * sizeof(int));
   memset ((void*)LowLev, '\0', MAXVARS * sizeof(int));
   memset ((void*)HiLev, '\0', MAXVARS * sizeof(int));

   Lvl = (*v5d_lvl_head); 	/* set local pointer to Level struct */

   /* Fill: 	NumTimes, DateStamp[], TimeStamp[]
   */
   for (NumTimes = 0, Msg=(*v5d_msgs_head) ; Msg != NULL; Msg=Msg->next) 
   {
      /* Loop until either matched, or list is exhausted */
      for (ind=0; ind < MAXTIMES ; ind++) 
      if (ValidEtimes[ind]==Msg->valid_etime || 
	  ValidEtimes[ind] == 0.0) break;

      /* If list is exhausted, Add this Base Time into our list,
	 Calculate its equivalent Date & Time stamp
      */
      if (ValidEtimes[ind]== 0.0) 		/* ADD IT */
      { 
	++NumTimes;
        /* convert to Human time */
	ValidEtimes[ind] = Msg->valid_etime; 
        e_h_time (&(Msg->valid_etime), &date, &hour); 
 	fprintf(stdout,
	"Time#%d (%lf) is %04d/%02d/%02d hr=%10.6lf",
	NumTimes, Msg->valid_etime, date.year, date.month, date.day, hour);

 	/* Calculate Date Stamp 'YYDDD' */
	md_doy (&date, &doy);		/* get Julian date from DATE struct*/
	DateStamp[ind] = ((date.year % 100) * 1000)  + doy;
	/*fprintf(stdout,"Julian day = %d ",doy);*/

 	/* Calculate Time Stamp 'HHMMSS' */
	temp = hour;
	TimeStamp[ind] = (int) floor(temp);   /* start with 'HH' */
 	temp = (hour - (double)TimeStamp[ind]) * 60.; /* move MM up */
	TimeStamp[ind] =  TimeStamp[ind]*100 + (int)floor(temp); /*tackon 'MM'*/
	/*fprintf(stdout,"MM moved up is %lf, TimeStamp is now='%d'\n", 
	temp,TimeStamp[ind]);*/
 	temp = (hour*100. - (double)TimeStamp[ind]) * 60.; /* move SS up */
	TimeStamp[ind] =  TimeStamp[ind]*100 + (int)(temp); /* tackon 'SS' */
	/*fprintf(stdout,"SS moved up is %lf", temp); */

	fprintf(stdout, " => '%05d' '%06d'\n", DateStamp[ind], TimeStamp[ind]);
      }
   } /* Parse unique times */

   /* Fill the VIS5D variables */
   Nr= v5d_info.nrows;
   Nc= v5d_info.ncols;
   CompressMode = 1;  	/* number of bytes per grid point */

   /* Set up VIS5D Projection (different from GRIB's) and its arguments */
   for (ind=0; ind<7; ind++) ProjArgs[ind] = v5d_info.proj_args[ind];
   Projection = v5d_info.v5d_proj;  /* Vis5d projection */

   DPRINT4 ("#Rows= %d, #Cols=%d, CompressMode=%d, Projection=%d\n",
   Nr, Nc, CompressMode , Projection);

   /* Set the vertical coordinate system */
   switch (Lvl->usLevel_id) {
     case 100:  Vertical= 3;  break;
     case 103:  Vertical= 2; break;
     case 105:  Vertical= 2; break;
     default: fprintf(stdout,"%s: Level %d not currently suppored\n", 
	func,Lvl->usLevel_id);
  	return (0);	/* return with no errors */
   }

    /* Vertical 3 means that the Vert_Args[] & data must go 
	from bottom to top.  So must flip the Height[] 
	but only do so if more than 1 heights in array.
    */
    if (Vertical == 3 &&  Lvl->numheights > 1) {
        fprintf(stdout,"Vertical=3,   Must Flip the Height[];\n");
	for (ind=0; ind < Lvl->numheights/2; ind++) {  
	   itmp  = Lvl->height[Lvl->numheights -1 -ind];
	   Lvl->height [Lvl->numheights -1 -ind] = Lvl->height[ind];
	   Lvl->height [ind] = itmp;
	  }
     }

   fprintf(stdout,"Height[%d]= (", Lvl->numheights);
   for (ind=0; ind<Lvl->numheights; ind++) 
	fprintf(stdout," %05d;", Lvl->height[ind]);
   	fprintf(stdout,")\n");

   /* Fill:    VertArgs[] array (Pressure 'mb' of Height)  */
   DPRINT1 ("VertArgs[%d]=(",Lvl->numheights);
   for (ind=0; ind<Lvl->numheights; ind++) {
      	VertArgs[ind]= (float)Lvl->height[ind]; 
	DPRINT1 (" %f;", VertArgs[ind]);
	}
   DPRINT0 (")\n");
  
/* Traverse Parameter Linked List, and fill:
	VarName []  name of this variable
	LowLev []   INDEX of the first Height w/in LVL->HEIGHT[]
	HiLev []    INDEX of the last Height w/in LVL->HEIGHT[]
	Nht []      count of steps w/in in LVL->HEIGHT[] going from 
		    LowLev to HiLev.
Levels:
001 => Surface
100 => Isobaric
102 => Mean Sea Level
105 => Height Level Above Ground

Parameters:
001 => Pressure (Pa) 			=> PRES
002 => Pressure reduced to MSL (Pa)	=> PRES_MSL
007 => Geopotential Height (gpm)	=> GEOP_HGT
011 => Air Temperature (K)		=> AIR_TEMP
031 => Wind Direction (deg)		=> WIND_DIR
032 => Wind Speed (m/s)		=> WIND_SPD
033 => U wind component (m/s)		=> U
034 => V wind component (m/s)		=> V
039 => Vertical Velocity (Pa/s)		=> OMEGA
040 => Vertical Velocity (m/s)		=> W
041 => Absolute Vorticity (s-1)		=> ABS_VORT
052 => Relative Humidity (%)		=> REL_HUM
061 => Total Precipitation		=> TTL_PRCP
083 => Land Cover (%)			=> TERRAIN
133 => Ground/Sea temperature 		=> SFC_TEMP
170 => Evaporative duct height		=> EVAP_DUCT

*/
   fprintf(stdout,"Vis5d summary:\n");
   for (MaxNht=-1, ivar=0, Parm=*v5d_parm_head ; 
		ivar < MAXVARS && Parm != NULL; Parm=Parm->next, ivar++)
    {  
      strcpy ((char *)Abbrv[ivar], Parm->abbrv); 

      /* Map Variable Name, Maximum of 9 chars */
      switch (Parm->usParm_id) 
	{
	case   1:  /* Pressure (Pa) */ 
		   strcpy ((char *)VarName[ivar], "PRES"); break;
	case   2:  /* Pressure reduced to MSL (Pa)*/ 
		   strcpy ((char *)VarName[ivar], "PRES_MSL"); break;
	case   7:  /* Geopotential Height (gpm)*/ 
		   strcpy ((char *)VarName[ivar], "GEOP_HGT"); break;
	case  11:  /* Air Temperature (K)*/ 
		   strcpy ((char *)VarName[ivar], "AIR_TEMP"); break;
	case  31:  /* Wind Direction (deg)*/ 
		   strcpy ((char *)VarName[ivar], "WIND_DIR"); break;
	case  32:  /* Wind Speed (m/s)*/ 
		   strcpy ((char *)VarName[ivar], "WIND_SPD"); break;
	case  33:  /* U wind component (m/s)*/ 
		   strcpy ((char *)VarName[ivar], "U"); break;
	case  34:  /* V wind component (m/s)*/ 
		   strcpy ((char *)VarName[ivar], "V"); break;
	case  39:  /* Vertical Velocity (Pa/s)*/ 
		   strcpy ((char *)VarName[ivar], "OMEGA"); break;
	case  40:  /* Vertical Velocity (m/s)*/ 
		   strcpy ((char *)VarName[ivar], "W"); break;
	case  41:  /* Absolute Vorticity (s-1)*/ 
		   strcpy ((char *)VarName[ivar], "ABS_VORT"); break;
	case  52:  /* Relative Humidity (%)*/ 
		   strcpy ((char *)VarName[ivar], "REL_HUM"); break;
	case  61:  /* Total Precipitation*/ 
		   strcpy ((char *)VarName[ivar], "TTL_PRCP"); break;
	case  83:  /* Land Cover (%)*/ 
		   strcpy ((char *)VarName[ivar], "TERRAIN"); break;
	case 133:  /* Ground/Sea temperature */ 
		   strcpy ((char *)VarName[ivar], "SFC_TEMP"); break;
	case 170:  /* Evaporative duct height*/ 
		   strcpy ((char *)VarName[ivar], "EVAP_DUCT"); break;
	default :  /* Name not avail, leave same as abbreviated IDs */
		   strcpy ((char *)VarName[ivar], (char *)Abbrv[ivar]); 
		   break;
	}
	
	LowLev[ivar]= Lvl->numheights;  /* set to big number for now */
	HiLev[ivar]= -1;  		/* set to low number for now */

	/* check every single message, look for the ones with this
	   ABBRV, and get its Min and Max Height Index
	*/
      	for (ind=1, Msg=(*v5d_msgs_head) ; Msg != NULL; Msg=Msg->next) {
	   /* skip if ABBRV is different */
	   if (strcmp (Msg->parm_ptr->abbrv, Abbrv[ivar])) continue;

	   /* get Array index of this height within the LVL's Height[] */
	   ind = get_htarray_index (Msg->usHeight1, Lvl->height, 
					Lvl->numheights);

	   /* Mark Lowest & Highest Index if needed */
	   if (ind < LowLev[ivar])  LowLev[ivar]= ind;
	   if (ind > HiLev[ivar])   HiLev[ivar]= ind;
	}

	Nht[ivar] = HiLev[ivar] - LowLev[ivar] + 1;   /* num of Ht steps */
	if (MaxNht < Nht[ivar]) MaxNht = Nht[ivar];   /* find max #steps */

	fprintf(stdout,
	"%-9s %s  has Heights [%05d (ht#%d) - %05d (ht#%d)] so Nht=%d\n",
	VarName[ivar], Abbrv[ivar],
	Lvl->height[LowLev[ivar]],  LowLev[ivar], 
	Lvl->height[HiLev[ivar]], HiLev[ivar], Nht[ivar]);
     }

   if (MaxNht == -1) {   /* double check */
	sprintf(errmsg,"%s:  MaxNht is -1\n", func);
	goto bail_out;
	}

   /* Print all message's info */
   DPRINT0("\nVis5d Message List=\n");
   for (ind=1, Msg=(*v5d_msgs_head) ; Msg != NULL; Msg=Msg->next) {
	/* Double Space if the Variable changes:
	   if (Msg->last != NULL && 
	   strcmp (Msg->parm_ptr->abbrv, Msg->last->parm_ptr->abbrv))
	   fprintf(stdout,"\n"); */
	DPRINT6 ( "#%d: Offs=%ld, %ld P1=%u;  Parm='%s', Height=%u\n",
	ind++, Msg->offset, Msg->base_dtg, Msg->usP1, Msg->parm_ptr->abbrv,
	Msg->usHeight1);
    }

   bytespergrid =  v5d_info.nrows * v5d_info.ncols * sizeof(float);

   /* Create an empty grid and fill it with value MISSING */
   if (! (EmptyGrid = (float*)malloc (bytespergrid))) {
	sprintf(errmsg, 
	"%s: failed to make an Empty Grid (R=%d * C=%d floats, or %ld bytes)\n",
	func, v5d_info.nrows, v5d_info.ncols, bytespergrid);
	goto bail_out;
	}
    else {
	for (itmp=v5d_info.nrows*v5d_info.ncols; itmp>=0;) 
	     EmptyGrid[--itmp]= MISSING;
	DPRINT1("Empty Grid has %d elements with value MISSING\n",
	v5d_info.nrows*v5d_info.ncols);
	}

   /* If no output file yet in then make up one: 'YYYYMMDDHH.LID.v5d' */
   if (outfile==NULL ) 
   {
      if (! (outfile = malloc(sizeof(20))))
	{ sprintf(errmsg,"%s: failed to create output filename\n", func);
	  goto bail_out;
	}
      sprintf(outfile, "%010d.%03d.v5d", v5d_info.base_dtg, Lvl->usLevel_id );
   }

   /* Create the v5d file and write the header (Zero means error) */
   if (v5dCreate  (outfile, NumTimes, NumVars, Nr, Nc, Nht,
		   VarName, TimeStamp, DateStamp, CompressMode,
                   Projection, ProjArgs, Vertical, VertArgs ) == V5_ERRCODE) 
    { 
      sprintf(errmsg, "%s:  v5dCreate failed to create %s\n",func,outfile);
      goto bail_out;
    }
   DPRINT1("v5dCreate() just finished, OUTFILE= '%s'\n", outfile);

   file_open = 1;   /* if abort somewhere then this flag will let us know
			to call V5dClose function
		    */

   /* Set the Bottommost Height value for all variables present */
   if (v5dSetLowLev (LowLev)==  V5_ERRCODE)
    { sprintf(errmsg, "%s:  v5dSetLowLev failed\n",func);
      goto bail_out;
    }
   DPRINT0("v5dSetLowLev() just finished\n");

   /*  Create a Data array required by Vis5d
       large enough to hold all grids from LOWLEV to HILEV
       of the Variable with the largest MAXNHT;
    */
   if (! (v5darr = (float*)malloc (bytespergrid * MaxNht))) { 
	sprintf(errmsg, 
	"%s: failed to make V5DARR (R=%d * C=%d * %d htsteps)\n",
	func, v5d_info.nrows, v5d_info.ncols, MaxNht);
	goto bail_out;
	}
    else {
	DPRINT3("\nStuff V5DARR has size (R=%d * C=%d * MaxNht=%d) elements\n",
	 v5d_info.nrows, v5d_info.ncols, MaxNht);
	}

    /* Load the data into the Vis5D data array:
    This Data is sorted by Height, then by Variable, then by Time.
    Each Data array to send Vis5d is of one single Time step only.
    for ex:
	Time X, Var 1, Ht 1 :  stuff data or Nulls
		       Ht 2 :  stuff data or Nulls
		       Ht 3 :  stuff data or Nulls
		Var 2, Ht 1 :  stuff data or Nulls
		       Ht 2 :  stuff data or Nulls
		       Ht 3 :  stuff data or Nulls
    */
    DPRINT0 ("START TO LOAD DATA:\n");
    for (itm=0;itm < NumTimes;itm++) /* For each TIME */
    {
     	DPRINT3 ("\n\n--> NEW TIME [itm=%d]=   YYJJJ=%d  HHMMSS=%d\n",
	itm, DateStamp[itm], TimeStamp[itm]);

	for (ivar=0; ivar < NumVars; ivar++)   /* For each VARIABLE */
	{
	   DPRINT4 ("--------> VAR[%d] '%s' %s has %d HTS (",
	   ivar,VarName[ivar], Abbrv[ivar], Nht[ivar]);
	   for (ind=0;ind<Nht[ivar];ind++) DPRINT1(" %d",Lvl->height[ind]);
	   DPRINT1 (")\n------------> Loop for %d HEIGHTS\n",Nht[ivar]);

	   for (iht=LowLev[ivar]; iht <= HiLev[ivar]; iht++) /*For each HEIGHT*/
	   {
	     /* Locate the Msg cell for curr Variable, curr valid 
		epochal Time and current Heigt
	     */
   	    for (Msg=(*v5d_msgs_head) ; Msg != NULL; Msg=Msg->next) 
	         if ( !strcmp(Msg->parm_ptr->abbrv, Abbrv[ivar]) 
		 && ValidEtimes[itm] == Msg->valid_etime 
		 && Msg->usHeight1 == Lvl->height[iht]) break;

	     /* find index within V5DARR[] where this Height begins */
	     ind = (iht - LowLev[ivar]) * Nc * Nr;

 	     if (Msg == NULL)    /* no match, Copy the Empty Grid*/
	      {
		memcpy ((void *)(v5darr+ind), (void *)EmptyGrid, bytespergrid);
		DPRINT8 ("%s #%d: no (%s Ht=%d), "\
		   "Stuff MISSING in v5darr[%d*%d*%d]= %d\n",
		   func,iht+1- LowLev[ivar], VarName[ivar],Lvl->height[iht],
		   iht-LowLev[ivar], Nc, Nr, ind);
		}   /* no match */
	     else 
		{ 		/* go get the matched msg's data */
    	        farr = (float *)NULL;
	        if ( getV5d_data (Msg->offset, &farr, errmsg)) 
    		  { upd_child_errmsg (func, errmsg); goto bail_out; }
	   	   
	        /* Copy Msg's data into v5d's data array */
	        memcpy ((void *)(v5darr + ind), (void *)farr, bytespergrid);

	        DPRINT9(
		" #%d: Stuff (%s Ht=%d %d P1=%d)in v5darr[%d*%d*%d]= %d\n",
		iht+1- LowLev[ivar], Msg->parm_ptr->abbrv, Msg->usHeight1, 
		Msg->base_dtg, Msg->usP1, iht-LowLev[ivar], Nc, Nr, ind);	

	        /* Release array returned by GetV5d_Data*/
	        if (farr != NULL) { free (farr); farr= NULL; }
	       }  /* get data */
	 }  /* HEIGHT loop */

         /* Write data to v5d file
	    This Data if of  1 time step, ordered by Height then by Variable
	 */
         if (v5dWrite( itm+1, ivar+1, v5darr) == V5_ERRCODE) {  
	    sprintf(errmsg,"%s:  v5dWrite returned error\n", func);
            goto bail_out;
         }
         DPRINT2( "v5dWrite(itm=%d, ivar=%d, v5darr) just finished\n", 
	 itm+1, ivar+1);
      }  /* VAR loop */
   }  /* TME loop */

  /* Close the v5d file */
   if (v5dClose() == V5_ERRCODE) {
	sprintf(errmsg,"%s:  v5dClose returned error\n",func); goto bail_out; }
   DPRINT0 ("v5dClose() just finished\n");

   /* clear out flag */
   file_open = 0;

   /* Change our Return Code */
   stat = 0;
   fprintf(stdout,"Vis5d file created= '%s'\n", outfile);
 
bail_out:

   /* Close up  Vis5d */
   if (file_open) v5dClose();

   if (v5darr != NULL) free(v5darr);
   if (farr != NULL)  free (farr); 
   if (EmptyGrid != NULL) free (EmptyGrid);
   DPRINT1("Exiting make_v5d_file with code=%d\n", stat);
   return(stat);
#endif
}

/*
*
***************************************************************
* B.  FUNCTION:  getV5d_data
*        Only get here if program was compiled with V5D_ON defined;
*        Loads a single GRIB message whose offset is provided
*        and returns the array of data in the storage order
*        as expected by Vis5D.
*
*     INTERFACE:   getV5d_data (Byte_offset, ARRAY, errmsg)
*
*     ARGUMENTS (I=input, O=output, I&O=input and output):
*       (I) long  Byte_offset;  #bytes fr. beginning of file
*       (O) float **ARRAY;      points to NULL upon entering function
*     (I&O) char  *errmsg;      holds err message if occurs
*
*     RETURN CODE:
*      0:  successful, ARRAY is returned pointing to a newly created array;
*          OR, flag V5D_ON is not available so no action was taken;
*      1:  failed, error buffer is filled;
***************************************************************/
int	getV5d_data (long Byte_offset, float **ARRAY, char *errmsg)
{
char *func="getV5d_data";
long offset = Byte_offset;
int  ind,nReturn;
float *grib_data = NULL;

  DPRINT1 ("Entering %s\n", func);
/*
* B.1       FUNCTION grib_seek    !find next GRIB msg starting at Offset
*           ! return if fails
*/
  if  (nReturn= grib_seek(InFile, &offset, 1, gh1, errmsg)) 
    	{ upd_child_errmsg (func, errmsg); goto bail_out; }
/*
* B.2       IF got a Warning from GribSeek THEN 
*                PRINT warning
*                RETURN with  no errors
*           ENDIF
*/
  /* NO errors but got a Warning msg from seek */
  if (errmsg[0] != '\0' || gh1->msg_length <= 0L)
   	{ sprintf(errmsg,"%s: grib_seek returned warning OR msglen<0;\n",func);
	goto bail_out;
	}
/*
* B.3       FUNCTION init_dec_struct   !init decoder structures
*/
      DPRINT1("Decoding message found at %ld bytes ...\n",offset);
      init_dec_struct(&pds,&gds,&bms,&bds_head);
/*
* B.4.      FUNCTION grib_dec      !perform GRIB decoding
*           !return if fails
*/
      grib_data=NULL;        
      if (nReturn = grib_dec ((char *)gh1->entire_msg, 
			&pds, &gds, &bds_head, &bms, &grib_data, errmsg)) 
    	{ upd_child_errmsg (func, errmsg); goto bail_out; }

/*
* B.5       IF (Bitmap Section is present) 
*           THEN
*             FUNCTION apply_bitmap !apply BMS to 'grib_data'
*             !return if fails
*           ENDIF
*/
      if (bms.uslength > 0  && 
	(nReturn= apply_bitmap(&bms,&grib_data, MISSING, &bds_head,errmsg)))
    	{ upd_child_errmsg (func, errmsg); goto bail_out; }

     /* print the first 3 rows & last 3 rows :
      { int ncols, rows_to_print=3;
	DPRINT1("Decoder returns (ulgridsize=%ld) :\n", bds_head.ulGrid_size);
	switch (gds.head.usData_type)  {
	  case LATLON_PRJ: ncols= gds.llg.usNj; break;
          case LAMB_PRJ:   ncols= gds.lam.iNy;  break;
	  }
        for (ind=0;  ind < ncols*rows_to_print; ind=ind+ ncols)
           DPRINT7 ("Row %d [%5d-]:  %f  %f  %f  %f  %f\n", 
	   ind/ncols, ind, grib_data[ind], grib_data[ind+1], grib_data[ind+2], 
	   grib_data[ind+3], grib_data[ind+4] );
        for (ind=bds_head.ulGrid_size-(ncols*rows_to_print);  
		ind < bds_head.ulGrid_size-1; ind=ind+ncols)
           DPRINT7 ("Row %d [%5d-]:  %f  %f  %f  %f  %f\n", 
	   ind/ncols, ind, grib_data[ind], grib_data[ind+1], grib_data[ind+2], 
	   grib_data[ind+3], grib_data[ind+4] );
     } */

/*
* B.6       FLIP the data to the orientation expected by Vis5D
*/

      switch ((int)gds.head.usData_type)
 	{
	 case LATLON_PRJ:  /* spherical grid */
            if (flip_data (gds.llg.usNi, gds.llg.usNj, 
			gds.llg.usScan_mode, grib_data, errmsg))
	        { upd_child_errmsg (func, errmsg); goto bail_out; } break;

         case LAMB_PRJ:  /* Lambert */
	    if (flip_data (gds.lam.iNx, gds.lam.iNy, 
			gds.lam.usScan_mode, grib_data, errmsg))
	        { upd_child_errmsg (func, errmsg); goto bail_out; } break;
	}
            

     /* print the first 3 rows & last 3 rows #/
     { int ncols, rows_to_print=3;
	DPRINT1 ("Flipped Data (ulgridsize=%ld) :\n", bds_head.ulGrid_size);
	switch (gds.head.usData_type)  {
	  case LATLON_PRJ: /# Spherical #/ ncols= gds.llg.usNj; break;
          case LAMB_PRJ: /# Lambert #/ ncols= gds.lam.iNy;  break;
	  }
        for (ind=0;  ind < ncols*rows_to_print; ind=ind+ ncols)
           DPRINT7 ("Row %d [%5d-]:  %f  %f  %f  %f  %f\n", 
	   ind/ncols, ind, grib_data[ind], grib_data[ind+1], grib_data[ind+2], 
	   grib_data[ind+3], grib_data[ind+4] );
        for (ind=bds_head.ulGrid_size-(ncols*rows_to_print);  
		ind < bds_head.ulGrid_size-1; ind=ind+ncols)
           DPRINT7 ("Row %d [%5d-]:  %f  %f  %f  %f  %f\n", 
	   ind/ncols, ind, grib_data[ind], grib_data[ind+1], grib_data[ind+2], 
	   grib_data[ind+3], grib_data[ind+4] );
     } */

/*
* B.7       ASSIGN Data array to user's ptr
*           RETURN addres of data array 	! successful exit
*/
      *ARRAY = grib_data;
      return (0);

bail_out:
/*
* B.8       !BAIL_OUT:
*           RETURN Null pointer; 	! exit with error
*/
      if (grib_data!=NULL)  free(grib_data); 
      *ARRAY = NULL;
      return (1);
/*
*  END of function
*/
}

/*
**************************************************************
* C. FUNCTION:   flip_data
*     To turn the mapping of the GRIB data to the orientation as
*     expected by the Vis5D program.   
*     Result grid will have (-y in +x) storage description=
*     first point at top left, data scanning top to bottom, left to right;
*
*   INTERFACE: 	flip_data (ncols, nrows, origfarr, errmsg)
*   ARGUMENTS (I=input, O=output):
*    (I)  int ncols;	  num of rows in grid 
*    (I)  int nrows;	  num of rows in grid 
*    (I)  unsigned short Scan_mode; scan mode of the data
*    (I)  float *farr;	  array with data to flip 
*    (I)  char *errmsg;   error buffer
* 
*   RETURN CODE:
*     0>  success, data has been flipped in array;
*     1>  Malloc error, see errmsg buffer.
**************************************************************
*/
#if PROTOTYPE_NEEDED
int    flip_data (int ncols, int nrows, unsigned short Scan_mode, 
			float *oldbuff, char *errmsg)
#else
int    flip_data (ncols, nrows, Scan_mode, oldbuff, errmsg)
int	 ncols;		/* num of cols in grid */
int	 nrows;		/* num of rows in grid */
unsigned short  Scan_mode;	/* scan mode of the data */
float 	*oldbuff;	/* array with data to flip */
char	*errmsg;	/* error buffer */
#endif
{
/*
* C.1        INITIALIZE variables
*/
char    *func="flip_data";
float    *newbuff = NULL;		/* temp newbuffer */
int	ind, r, c, size;		/* working vars */

    DPRINT1("Entering %s\n", func);
/*
* C.2        CREATE a temp array (#Rows by #Cols)
*            !return if failed
*/
    size = nrows * ncols * sizeof(float);   
    if ((newbuff = (float *)malloc (size)) == NULL) {
	sprintf(errmsg,"%s: failed to create temp newbuffer\n",func);
        DPRINT1("Leaving %s with malloc errors\n", func);
	return (1);
	}
/* original existing code (64)
* C.3        DEPENDING on the Scan Mode of the data
*               case   0: 
*                 !switch from GRIB's 1st point is top Left (+i,-j, Horiz first)
*                 !to V5D's 1st point is top Left (+i,-j, Vertic first)
*               case  32: 
*                 !do nothing, GRIB's 1st point is same as V5D's
*                 !1st point is already at top Left (+i,-j, Vertic first)
*               case  64:     
*                 !switch from GRIB's 1st point bottom left (+i, +j, Horz first)
*                 !to V5D's 1st point is top Left (+i,-j, Vertic first)
*               case  96: 
*                 !switch from GRIB's 1st point bottom Left (+i,+j,Vertic first)
*                 !to V5D's 1st point is top Left (+i,-j, Vertic first)
*               case  128: 
*                 !switch from GRIB's 1st point Top Right (-i,-j, Horiz first)
*                 !to V5D's 1st point is top Left (+i,-j, Vertic first)
*               case  192: 
*                 !switch from GRIB's 1st point Top Right (-i,-j, Vertic first)
*                 !to V5D's 1st point is top Left (+i,-j, Vertic first)
*               case  160: 
*                 !switch from GRIB's 1st point Bottom Right (-i,+j,Horiz first)
*                 !to V5D's 1st point is top Left (+i,-j, Vertic first)
*               case  224: 
*                 !switch from GRIB's 1st point Bottom Rt (-i,+j,Vertic first)
*                 !to V5D's 1st point is top Left (+i,-j, Vertic first)
*            ENDTEST
*/

/* 
   newbuff[] holds the result data after the flip the way Vis5d expects
   ( 1st point is top Left (+i,-j, Vertic first) )
*/
    switch ((int) Scan_mode) {
	case  0:  /* GRIB's 1st point is top Left (+i,-j, Horiz first) */
	    for (ind=0, c=0; c < ncols; c++)
            for ( r=0; r < nrows; r++)
            newbuff[ind++] = oldbuff [(r * ncols) + c];
            break;

	case 32:  /* GRIB's 1st point is top Left (+i,-j, Vertic first) */
            break; /* do nothing !!! */

	case 64: /* GRIB's 1st point is bottom left (+i, +j, Horz first) */ 
		 /* ***ORIGINAL flip algorithm*** "x in +y" to "-y in +x" */
	    for (ind=0, c=0; c < ncols; c++) 
	    for ( r=0; r < nrows; r++)  
	    newbuff[ind++] = oldbuff [((nrows-r-1) * ncols) + c];
	    break;

	case 96:  /* GRIB's 1st point is bottom Left (+i,+j, Vertic first) */
	    for (ind=0, c=0; c < ncols; c++)
            for ( r=0; r < nrows; r++)
            newbuff[ind++] = oldbuff [(nrows-r-1) + (c * nrows)];
            break;

	case 128: /* GRIB's 1st point is Top Right (-i,-j, Horiz first) */
	    for (ind=0, c=0; c < ncols; c++)
            for ( r=0; r < nrows; r++)
            newbuff[ind++] = oldbuff [(r * ncols) + (ncols-c-1)];
            break;

	case 192: /* GRIB's 1st point is top Right (-i,-j, Vertic first) */
	    for (ind=0, c=0; c < ncols; c++)
            for ( r=0; r < nrows; r++)
            newbuff[ind++] = oldbuff [ r + (ncols-c-1) * nrows];
            break;

	case 160:  /* GRIB's 1st point is Bottom Right (-i,+j, Horiz first) */
	    for (ind=0, c=0; c < ncols; c++)
            for ( r=0; r < nrows; r++)
            newbuff[ind++] = oldbuff [(nrows-r-1) * ncols + (ncols-c-1)];
            break;

	case 224:  /* GRIB's 1st point is Bottom Right (-i,+j, Vertic first) */
	    for (ind=0, c=0; c < ncols; c++)
            for ( r=0; r < nrows; r++)
            newbuff[ind++] = oldbuff [(nrows-r-1) + (ncols-c-1)*nrows];
            break;
    }
	
/*
* C.4         PUT the content of the flipped array into the old array 
*/
    memcpy ((void*)oldbuff, (void*)newbuff, size);

/*
* C.5        RELEASE temporary space
*/
    if (newbuff!=NULL)  free (newbuff);

/*
* C.6        RETURN with no errors
*/
    DPRINT1("Leaving %s with no errors\n", func);
    return (0);
}

/*
OLD>**************************************************************
OLD>* D. FUNCTION:   old_flip_data
OLD>*     To turn the mapping of the Grid to the orientation 
OLD>*     from "+x in -y" to "x in +y"
OLD>**************************************************************
OLD>#/
OLD>#if PROTOTYPE_NEEDED
OLD>int  old_flip_data (int ncols, int nrows, float *farr, char *errmsg)
OLD>#else
OLD>int  old_flip_data (ncols, nrows, origfarr, errmsg)
OLD>int	 ncols;		/# num of rows in grid #/
OLD>int	 nrows;		/# num of rows in grid #/
OLD>float 	*farr;		/# array with data to flip #/
OLD>char	*errmsg;	/# error buffer #/
OLD>#endif
OLD>{
OLD>char    *func="flip_data";
OLD>float    *buff = NULL;		/# temp buffer #/
OLD>int	top,bot, r, c, size;		/# working vars #/
OLD>
OLD> 	DPRINT1("Entering %s\n", func);
OLD>/#
OLD>	>>> flipping [left to right, top to bottom] to
OLD>	[left to right, bottom to top],  a row at a time.
OLD>* D.2        FLIP the data, "+x in -y" to "x in +y"
OLD>#/
OLD>    size = ncols * sizeof(float);   /# number of bytes per ROW #/
OLD>    if ((buff = (float *)malloc (size)) == NULL) {
OLD>	sprintf(errmsg,"%s: failed to create temp buffer\n",func);return(1);}
OLD>
OLD>    /#  put Top half of grid into Bottom half of grid, 1 line at a time; 
OLD>	stop at middle else will end up with original grid! #/
OLD>    for (r=0; r < nrows/2; r++) {
OLD>	/# copy entire row to temp buffer #/
OLD>	top = r *ncols; bot = ((nrows-1-r) *ncols);
OLD> 	/# if (r < 5) {
OLD>	   fprintf(stdout,"TOP = r * ncols= (%d * %d) = %d\n",r,ncols,top);
OLD>	   fprintf(stdout,"BOT = (nrows-1-r)*ncols= (%d-1-%d)*%d= %d\n",
OLD>		nrows,r,ncols,bot);
OLD>   	   fprintf(stdout,"switching line %d with %d\n",top,bot);
OLD>	   fprintf(stdout,"#%5d: %f %f %f %f ...\n",
OLD>	   top,*(farr+top), *(farr+top+1), *(farr+top+2), *(farr+top+3));
OLD>	   fprintf(stdout,"#%5d: %f %f %f %f ...\n",
OLD>	   bot,*(farr+bot), *(farr+bot+1), *(farr+bot+2), *(farr+bot+3));
OLD>	  } #/
OLD>	memcpy ((void*)buff, (void*)(farr+top), size);
OLD>	memcpy ((void*)(farr+top), (void*)(farr+bot), size);
OLD>	memcpy ((void*)(farr+bot), (void*)buff, size);
OLD>    /#if (r < 5) {
OLD>	   fprintf(stdout,"After MEMCPY,\n");
OLD>	   fprintf(stdout,"#%5d: %f %f %f %f ...\n",
OLD>	   top,*(farr+top), *(farr+top+1), *(farr+top+2), *(farr+top+3));
OLD>	   fprintf(stdout,"#%5d: %f %f %f %f ...\n",
OLD>	   bot,*(farr+bot), *(farr+bot+1), *(farr+bot+2), *(farr+bot+3));
OLD>	  } #/
OLD>	}
OLD>/#
OLD>* D.3        RELEASE temporary space
OLD>#/
OLD>    if (buff!=NULL)  free (buff);
OLD>/#
OLD>* D.4        RETURN with no errors
OLD>#/
OLD>    DPRINT1("Leaving %s with no errors\n", func);
OLD>    return (0);
OLD>}
*/
