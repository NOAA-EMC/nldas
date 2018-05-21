/*........................................................
Filename:  ld_v5d_msg.c
Author  :  Alice T. Nakajima, SAIC MRY
Date    :  11/18/96
*/
#include <stdio.h>
#include <math.h>
#include "grib_lookup.h"  /* for all lookup structs */
#include "isdb.h"  	  /* struct date */
#include "dprints.h"	/* for dprints */
#include "gribfuncs.h"	/* prototypes */
#include "gsv5d.h"	  /* local vis 5d structs & links to v5d.h too */

static int  unit_offset[5]= 
 { 1, 100, 10000, 1000000, 100000000 }; /* hr, month, day, year, cent */

/*
*
* ======================================================================
* A.  FUNCTION:  ld_v5d_msg
*
*     PURPOSE: 
*     Called on to store info of the Msg just read into the arrays
*     that will be used to create Vis5d output file;
*     - does not allow duplicates of messages;
*     - gives preference to msg of certain parm&Lvl&Ht of Earlier DTGs
*
*     ARGUMENTS (I)=input, (O)=output, (I&O)=input and output:
*       (I) long       GRIB_offs;    
*           num of bytes from beginning of file to the start of this GRIB msg
*       (I) PDS_INPUT  pds;          
*           Product description section of this message
*       (I) grid_desc_sec  gds;      
*           Grid descr section of this message
*       (I) BDS_HEAD_INPUT bds_head; 
*           Binary Data section of this message
*       (I) BMS_INPUT  bms;          
*           Bitmap section of this message
*     (I&O) int *tot_parms;   
*           count of unique parameter types, may get updated;
*     (I&O) int *tot_msgs;    
*           count of unique Messages, may get updated;
*     (I&O) V5D_LVL **lvl_head;   
*           head of Level linked list, may get updated;
*     (I&O) V5D_PARM **parm_head;  
*           head of Parameter linked list, may get updated;
*     (I&O) V5D_PARM **parm_tail;  
*           tail of Parameter linked list, may get updated;
*     (I&O) V5_REC **msgs_head;  
*           head of Grib Msg linked list, may get updated;
*     (I&O) V5_REC **msgs_tail;  
*           tail of Grib Msg linked list, may get updated;
*     (I&O) V5_INFO  *grad_info;   
*           holds common info, gets updated if this is 1st msg in v5d list
*       (O) char *errmsg;      
*           holds err message if occurs
*
*     RETURN CODE:
*       0>  no errors, may/may not have stored this msg;
*             (Duplicate message, out of Exntensions, invalid forecast unitid)
*          -message is inserted in Msg linked list sorted by Base time,
*           Forecast Period, Parameter and Leveltype;
*          -if new Level type then add it to Level Linked list sorted
*           by level_id;
*          -if new Height then add it to Height array of this Levelid,
*           sorted in ascending order;
*          -if new parameter type then add it to Parameter Linked list
*           sorted by parm_id;
*       1>  error, errmsg gets filled;
* ======================================================================
*/
extern int 		UseTables;
extern PARM_DEFN        db_parm_tbl[];

#if PROTOTYPE_NEEDED
   int  ld_v5d_msg (
		long 		GRIB_offs,
		PDS_INPUT 	pds, 
		grid_desc_sec 	gds, 
		BDS_HEAD_INPUT 	bds_head,
		BMS_INPUT 	bms,
		int 		*tot_parms, 
		int 		*tot_msgs,
		V5D_LVL  	**lvl_head, 
		V5D_PARM 	**parm_head, V5D_PARM 	**parm_tail,
		V5_REC  	**msgs_head, V5_REC 	**msgs_tail,
		V5_INFO 	*v5_info,
		char 		*errmsg)
#else
   int  ld_v5d_msg ( GRIB_offs, pds, gds, bds_head, bms,
		tot_parms,  tot_msgs, lvl_head, parm_head, parm_tail, 
		msgs_head, msgs_tail, v5_info, errmsg)
		
		long 		GRIB_offs;
		PDS_INPUT 	pds; 
		grid_desc_sec 	gds; 
		BDS_HEAD_INPUT 	bds_head;
		BMS_INPUT 	bms;
		int 		*tot_parms, *tot_msgs;
		V5D_LVL  	**lvl_head; 
		V5D_PARM 	**parm_head, **parm_tail;
		V5_REC  	**msgs_head, **msgs_tail;
		V5_INFO 	*v5_info;
		char 		*errmsg;
#endif
{
char    *func="ld_v5d_msg";
char 	abbrv[9];	/* work var */
int	Parmtbl_indx;   /* within the Parameter Conversion Array */
int     Recycling_cell = 0;  /* default of not recycling Cell */
int	lvl_indx, n,i,x, unit_incr;
int     newparm;	/* set if need to add new parm */
int	newlvl=0;	/* set if need to add new Level */
int	newheight=0;	/* set if need to add new height */
int	ht_index=-1;	/* indx within lvl's ht array */
int     stat=0;
long	base_dtg;	/* work var */
V5_REC   *msg_ptr; 	/* temp pointer */
V5_REC   *NewMsg=NULL;   /* cell for current message */
V5D_PARM  *parm_ptr;    /* temp pointer */
V5D_PARM  *insrt_parm;  /* working ptr */
V5D_PARM  *TmpParm;	 /* is Parm Cell that matched/got created */
V5D_LVL   *Lptr;        /*  is Lvl Cell that matched/got created */
DATE    date; 
double  valid_etime;    /* (Ref.time + FcstPer) in Epochal Time */
double  etime;   	/* Reference Time in epochal time */
double	hour;		/* the Hour, Min, Sec converted to Hours */

/*
* A.1          IF (vis5d option is not available) THEN
*                 PRINT message
*                 RETURN with no errors
*              ENDIF
*/
#ifndef V5D_ON
	fprintf(stdout,"Inside %s:  skip because V5D_ON is not available\n",
	func);
	return (0);
#else

/*
*
* A.2          INIT variables
*/
      DPRINT6 ("Inside ld_v5d_msg\nLoad: %02d%02d%02d%02d%02d fcst=%02d ",
                  pds.usCentury-1, pds.usYear,pds.usMonth, pds.usDay, 
                 pds.usHour, pds.usP1);
      DPRINT4 ("P=%03d (subP=%03d) L=%03d Ht=%03d;", 
	pds.usParm_id, pds.usParm_sub, pds.usLevel_id, pds.usHeight1);

     /*..........................................................
	Since Vis5d only works with one single Leveltype,
	we do not have a Level Linked list.
	Instead, we have just 1 record with all its info, and
	the Height values are sorted in ascending order.

        if Head Ptr is not NULL then this LEVEL already exists,  
              -check if curr ht is already in height[] array;
		add to array if not.
        else if Head Ptr is NULL:  add new lvl type;
     ..........................................................*/
/*
*
* A.3          IF (level_id already exists) 
*/
      Lptr	= *lvl_head;	/* default */
      if (Lptr != NULL)
	{    
	/* LEVL EXISTS, ADD HEIGHT IF NEW  
	*/
/*
* A.3.a        THEN
*                 DEFAULT as an old level, new height
*/
	   newlvl=0;
	   newheight= 1;
/*
* A.3.a.1         FOR (each Height of this Level) DO
*/
           for (ht_index=-1, x=0, i=0; i < (int)Lptr->numheights; i++)
	    {
/*
* A.3.a.1.a           IF (height is greater than currmsg's height) BREAK;
*/
		if (Lptr->height[i] > pds.usHeight1) break;

/*
* A.3.a.1.b           ELSE IF (height is smaller than currmsg's height) 
*                         MARK ht_index to curr index; !insrt after this index
*/
		else if (Lptr->height[i] < pds.usHeight1) 
		     	ht_index= i;  /* insrt AFTER this index */

/*
* A.3.a.1.c           ELSE IF (height matches currmsg's height) THEN
*                         SET newheight flag to 0
*                         BREAK out of loop
*/
		else if (Lptr->height[i]==pds.usHeight1)
                     {  newheight= 0;
			break;      /* Ht already defined, do nothing*/
		     }
/*
* A.3.a.1.c           ENDIF
*/
		/**** 
		DONOT CHANGE 'HT_INDEX' AFTER THIS POINT; 
		*****
		Note:  will only add new height later if we're keeping 
		this message
		*/
/*
* A.3.a.1         ENDFOR
*/
	    }
        }
/*
* A.3.b        ELSE
*/
      else 			/* First time around  */
         {    
	   /* 	 ADD NEW LEVEL (& HEIGHT)
		 -- Make Lptr the new cell to be added;
	    */
/*
* A.3.b.1         SET flags, newlvl=1 and newheight=0
*/
	   newlvl=1;
	   newheight= 0;
/*
* A.3.b.2         ALLOCATE storage for new level cell
*                 IF (failed) THEN
*                    PRINT message
*                    RETURN 1		!malloc error
*                 ENDIF
*/
           if ((Lptr=(V5D_LVL *)malloc(sizeof(V5D_LVL)))==NULL) { 
		 sprintf(errmsg,"ld_v5d_msg: Failed to malloc new level \n");
                 stat=(1); goto BYE;
		}

/*
* A.3.b.3         FILL Temp level cell  !leave unlinked to List for now
*/
	   if (v5_info->zdef_lvl == 105 && 
		(pds.usLevel_id==1 || pds.usLevel_id==102)) 
	   {
	     DPRINT2("%s: Level %03d gets treated as Level 105 at Ht=0\n",
	     func, pds.usLevel_id);
             Lptr->usLevel_id     = 105;
             Lptr->height[0]      = 0;
	   }
	  else {
	     Lptr->usLevel_id     = pds.usLevel_id;
             Lptr->height[0]      = pds.usHeight1;
	   }

           Lptr->numheights     = 1;
	   /* 
		Note:  link Lptr to list later if we're definately keeping
		this message 
	   */
/*
* A.3.b        ENDIF
*/
        }
     /**** 	
      DONOT CHANGE 'LPTR' AFTER THIS POINT;
      ****/


/*
*
* A.4          CHECK if parm already exists;
*              ! SET Parameter flags and ptrs to default values
*/
   sprintf (abbrv, "P%03dS%03d", pds.usParm_id, pds.usParm_sub);

   insrt_parm=NULL;	/* default */
   TmpParm = NULL;	/* default */
   newparm = 1;		/* default */

/*
* A.5          FOR (each cell in Parameter Linked List) DO
*/
     for (parm_ptr=*parm_head;  newparm && parm_ptr!=NULL;  
						parm_ptr=parm_ptr->next)
     { /* check if Parm already exists */
/*
* A.5.1             IF (cell's ABBRV matches current ABBRV)
*                   THEN
*                      MARK not a new parameter
*                      DROP out of For-loop
*                   ELSE
*                      MARK insrt_parm as place to insert after
*                      LOOP, check next cell
*                   ENDIF
*/
	/* Vis5d:  Keep Unique by ABBRV, 'PpidSsid' 
	   Keeping it in Ascending order
	*/
	  x=  (int)strcmp(parm_ptr->abbrv, abbrv);

	  if (x==0) 
	  /* matched, not a new parm */  
		{ newparm=0; break; }
	  else if (x < 0) 
	  /* Cell's Abbrv < curr Abbrv, so Insert AFTER Insrt_parm */
		insrt_parm=parm_ptr; 
	  else break;		
	  /* Cell's Abbrv > curr Abbrv, insert after PREVIOUS insrt_parm */ 
/*
* A.5          ENDFOR  ! traverse list
*/
     }  /* for */


     /*....................................................
        create linked list of unique VAriable;
        .. if lvl is Zdef levl then each var is a unique
           combo of (parmid, lvlid);
        .. if lvl is non-Zdef lvl then each var is a unique
           combo of (parmid, lvlid, height);
        Only create if desired combo is not in list;
      *...................................................*/
/*
*
* A.6         IF (this is an old parameter) THEN
*                SET TmpParm point to that location in parm list
*/
      if (! newparm)	
	{			/* TmpPtr is the one that matched */
	  TmpParm= parm_ptr;	/* already exists, no set up needed */
	}
/*
* A.6.a       ELSE
*/
      else   			/* Set up new param AFTER Insrt_Ptr */
	{		/* Here, TmpPtr is new parmcell just created */

/*
* A.6.a.0        IF reached Vis5d's max number of Variables
*                    DROP msg, fill errmsg
*                    RETURN
*                ENDIF
*/
	 if (*tot_parms + 1 >= MAXVARS) {
	    sprintf(errmsg,"%s: cannot add '%s', excced VIS5D's MAXVARS\n",
	    func, abbrv);
	    goto BYE;
	    }
/*
* A.6.a.1        ALLOCATE storage for the new parm
*                IF (error) RETURN 1
*/
          if ((TmpParm=(V5D_PARM *)malloc(sizeof(V5D_PARM)))==NULL)
          { sprintf(errmsg,"ld_v5d_msg: Failed to malloc V5D_PARM\n");
	    stat=(1); goto BYE; }
/*
  Not counting until actually linking it to the Parm List
GRADS>* A.6.a.2        INCREMENT parameter counter
GRADS>-/
GRADS>	  *tot_parms += 1;
*/
	DPRINT1("ld_v5d:  make new cell for new parm '%s'\n",abbrv);
/*
* A.6.a.3        FILL this temp Parameter cell
*                !leaving it unlinked to Parameter list for now
*/
	  strcpy (TmpParm->abbrv, abbrv);
          TmpParm->usParm_id	= pds.usParm_id;

          if (UseTables) 
	    {   /* sub_tab2[pds.usParm_sub]. */

	    Parmtbl_indx =  PARMTBL_INDX(pds.usParm_id, pds.usParm_sub);

	    if (Parmtbl_indx && Parmtbl_indx % 256 == 0) 
		strcpy (TmpParm->varnm, "[Reserved Code, not used]");

	    else if ( ! db_parm_tbl[Parmtbl_indx].grib_dsc[0] ||
		      ! db_parm_tbl[Parmtbl_indx].grib_unit_dsc[0])
		strcpy (TmpParm->varnm, "[Undefined Variable]");

	    else {
                /* 40 chars max:  'ParameterName [unit]' */
                n= strlen(db_parm_tbl[Parmtbl_indx].grib_unit_dsc);
                strncpy (TmpParm->varnm, 
			db_parm_tbl[Parmtbl_indx].grib_dsc, 
			40-n-3);
                strcat (TmpParm->varnm," [");
                strcat (TmpParm->varnm, 
			db_parm_tbl[Parmtbl_indx].grib_unit_dsc);
                strcat (TmpParm->varnm, "]");
    		}
            }
          else strcpy (TmpParm->varnm, "[No Lookup File]");
	  TmpParm->lvl_ptr = Lptr;  /* pts to Level cell ABOVE */
	  /* Note:
		Link cellto list later only if keeping this Msg
	  */
/*
* A.6.a       ENDIF
*/
        }  /* NewParm */
 	/*************************************
	DONOT CHANGE 'INSRT_PARM' OR 'TmpParm' FROM HERE ON
	************************************/

/*
* A.7         CALCULATE this msg's base time (in hours)
*/
  base_dtg =  (int)(pds.usCentury-1)*unit_offset[4] 
		+ (int)pds.usYear*unit_offset[3]
		+ (int)pds.usMonth*unit_offset[2] 
		+ (int)pds.usDay*unit_offset[1] 
		+ (int)pds.usHour;

  
  /* 
   Find time when field is Valid in Epochal time (Ref time + Fcst Period) 
  */
  date.year  = (int) ((pds.usCentury-1)*100+pds.usYear);
  date.month = (int) pds.usMonth;
  date.day   = (int) pds.usDay;
  if (pds.usSecond != 0 && pds.usSecond <= 60)
           hour = (double)pds.usHour+(double)pds.usMinute/60.+
		  (double)pds.usSecond/3600.;
      else hour   = (double)pds.usHour + (double)pds.usMinute/60.;
   h_e_time (&date, &hour, &etime); /* Reference Time in Epochal time */

   /* Now add the Fcst Period to (DATE) structure , Hour
   */
   if  (pds.usFcst_unit_id != 1) 
	fprintf(stdout,"%s warning:  Forecast Time Unit %d is not in HOURS\n",
	func,pds.usFcst_unit_id);

   DPRINT1("Reference time in Epochal Time = %lf ",etime);
   switch  (pds.usFcst_unit_id)
        {
        case 0: DPRINT1("Plus usP1=%d  MINUTES \n",pds.usP1);
                 hour += (double)pds.usP1 / 60.; break;
        case 1: DPRINT1("Plus usP1=%d  HOURS \n",pds.usP1);
                 hour += (double)pds.usP1; break;
        case 2: DPRINT1("Plus usP1=%d  DAYS \n",pds.usP1);
                 date.day += (int)pds.usP1; break;
        case 3: DPRINT1("Plus usP1=%d  MONTHS \n",pds.usP1);
                 date.month += (int)pds.usP1;
                while (date.month > 12) { ++date.year; date.month -= 12; }
                break;
        case 4: DPRINT1("Plus usP1=%d  YEARS \n",pds.usP1);
                 date.year += (int)pds.usP1; break;
        case 5: DPRINT1("Plus usP1=%d  DECADES (10 years) \n",pds.usP1);
                date.year += (10 * (int)pds.usP1); break;
        case 6: DPRINT1("Plus usP1=%d  NORMAL (30 years) \n",pds.usP1);
                date.year += (30 * (int)pds.usP1); break;
        case 7: DPRINT1("Plus usP1=%d  CENTURY (100 years) \n",pds.usP1);
                date.year += (100 * (int)pds.usP1); break;
        case 254: DPRINT1("Plus usP1=%d  SECONDS\n",pds.usP1);
                hour += ((float)pds.usP1 / 3600.); break;
        default: fprintf (stdout,
               "ld_v5d_msg warning: Invalid Forecast Unit Id (%d), drop msg\n"
                ,pds.usFcst_unit_id);
                break;
        }

   /* 	get (RefTime + Fcstper) in Epochal time 
      NOTE:
      The values in DATE, hour may not be in proper range (24 hrs/day, 
      60 mins/hr, etc...).  It depends on which the Fcst Unit Time is in.
      H_E_TIME will sort it out and give back the Valid time in epochal.
      So from here on, it may not be a good idea to use the DATE struc & HOUR
   */
      h_e_time (&date, &hour, &valid_etime); 


/*
* A.8.b.5        FOR (each cell in the Message list) DO
*/
      for (msg_ptr=*msgs_head;  msg_ptr!= NULL;   msg_ptr=msg_ptr->next)
      {
/*
* A.8.b.5.1        !Compare the new msg with the current Cell in Msg list
*                  IF (either Valid time OR Abbrv OR Height value 
*                      is different)
*                  THEN
*                     GO check the next cell in msg list
*                  ENDIF
*/
         if ( valid_etime != msg_ptr->valid_etime 
	     ||  strcmp (abbrv, msg_ptr->parm_ptr->abbrv)
	     ||  pds.usHeight1 != msg_ptr->usHeight1 ) 	  
 	 continue;

/*  Same (Valid time, Abbrv, and Height value)  so far: 
 *  Now choose the earlier DTG...  
*/

/*
* A.8.b.5.2        IF (newmsg's basetime is same as cell's basetime)
* A.8.b.5.2.a      THEN   !this is a duplicate
*                     FREE temp Parm cell if there is one
*                     FREE temp Level cell if there is one
*                     PRINT message
*                     RETURN 2 !drop msg
*/
         if (base_dtg == msg_ptr->base_dtg)  /* same Base, same FcstPer */
	  {
	    if (newparm) free(TmpParm);
	    if (newlvl) free (Lptr);
	    fprintf (stdout,
   	    "ld_v5d_msg:  drop Duplicate [%d tau=%d Parm=%u ht=%u]\n", 
	    base_dtg, pds.usP1, pds.usParm_id, pds.usHeight1);
	    goto BYE;
	  }
/*
* A.8.b.5.2.b      ELSE IF (newmsg's basetime is earlier than cell's time)
*                  THEN
*/
         else if (base_dtg < msg_ptr->base_dtg) /* curr Msg has smaller DTG */
	  {
/*
* A.8.b.5.2.b.1       DEBUG print  
*/
	    DPRINT3("%s: Remove old '%ld' valid_etime=%lf\n",
	    func, msg_ptr->base_dtg,msg_ptr->valid_etime);
	   /*
            "\nCurrent Msg [DTG=%ld, FcstPer=%d , '%s' ht=%d] "
	    "has smaller DTG than:\n"
            "stored Msg [DTG=%ld, FcstPer=%d , '%s' ht=%d] \n"
	    "REMOVING the Previously Stored Msg & keeping current one\n"
		    ,base_dtg, pds.usP1,  abbrv, pds.usHeight1,
		    msg_ptr->base_dtg,msg_ptr->ustau,
		    msg_ptr->parm_ptr->abbrv, msg_ptr->usHeight1);
	    */
/*
* A.8.b.5.2.b.2       UNLINK cell with bigger basetime from Msg list
*/
	    /* remove the previously stored msg w/ bigger dtg */
	    if (msg_ptr != *msgs_head && msg_ptr!= *msgs_tail)
	         {   
		   /* remove cell between head and tail */
		   msg_ptr->next->last= msg_ptr->last;
		   msg_ptr->last->next= msg_ptr->next;
	         }
	    else { 
		   if (msg_ptr==*msgs_head)  /* rm 1st cell */
		   { *msgs_head=msg_ptr->next;
		     if (*msgs_head !=NULL) (*msgs_head)->last= NULL;
		   }
		   if (msg_ptr== *msgs_tail) /* rm last cell */
		   { *msgs_tail= msg_ptr->last;
		     if (*msgs_tail !=NULL) (*msgs_tail)->next= NULL;
		   }
		}
/*
GRADS>* A.8.b.5.2.b.3       FREE its storage  !Insrt new Cell here later
GRADS>*
GRADS>            free (msg_ptr); /- remove old cell , insert new cell later here-/
*/

/*
* A.8.b.5.2.b.3       SET pointer to cell to be recycled & set Recycle flag;
*                     !fill it with New Msg's info later
*/
            NewMsg = msg_ptr;	   /* set pointer to cell to be recycled */
	    Recycling_cell = 1;  
/*
* A.8.b.5.2.b.4       BREAK out of For loop !all done
*/
	    break;		/* all done */
	  }
/*
* A.8.b.5.2.c      ELSE  ! Current msg's has bigger DTG
*/
         else	/* curr Msg has bigger DTG */
	  {
/*
* A.8.b.5.2.c.1       PRINT message
*/
	    DPRINT1 ("%s: bigger dtg, Ignored\n", func);

	    fprintf(stdout,
            "\nCurrent Msg [DTG=%ld, FcstPer=%d (ValidEtm=%lf), '%s' ht=%d]\n"
	    "has bigger DTG than previously stored message;\n"
            "Previous Msg [DTG=%ld, FcstPer=(n/a) (ValidEtm=%lf), '%s' Ht=%d] "
	    "\nvis5d not storing current message...\n",
	    base_dtg, pds.usP1, valid_etime, abbrv, pds.usHeight1,
	    msg_ptr->base_dtg, msg_ptr->valid_etime, 
	    msg_ptr->parm_ptr->abbrv, msg_ptr->usHeight1);

/*
* A.8.b.5.2.c.2       FREE temp Parm cell if there is one
*                     FREE temp Level cell if there is one
*/
	    if (newparm) free (TmpParm);
	    if (newlvl)  free (Lptr);
/*
* A.8.b.5.2.c.3       RETURN 0 !not going to store later DTG
*/
	    goto BYE;  /* skip msg */
/*
* A.8.b.5.2.c      ENDIF
*/
	  }
/*
* A.8.b.5        ENDFOR !traverse list
*/
    } /* traverse list */
/*
* A.8.b        ENDIF
*/
	    
/*    
*              ! Not a duplicate msg, ok to keep it;
*/

/* ###*#*#*#*#  *#*#*#*#  #*#*#*#*#  *#*#*#*#  
	USED TO INSERT PARM FIRST, THEN HEIGHT, THEN LVL 
	TRY TO INSERT PARM 'AFTER' HEIGHT, LEVEL
 ##*#*#*#*#  #*#*#*#*#  *#*#*#*#  #*#*#*#*#  
*/

/*--- new insert new height in our Level rec if needed */
/*
*
* A.10         IF (this is a new Height)
*              THEN
*/
    if (newheight) {
	  DPRINT0 (" +Ht");
/*
* NEW for VIS5D>>>
* A.10.0           IF (reached max number of Heights defined by Vis5d) THEN
*                     PRINT message
*                     RETURN with errors
*                  ENDIF
*/
	if (Lptr->numheights+1 == MAXLEVELS) {
	   sprintf(errmsg,
	   "%s: cannot add Height %d, exceeded VIS5D's MAXLEVELS\n",
	   func, pds.usHeight1);
	   goto BYE;
	}
/*
* A.10.1           IF (Inserting Height at beginning of list) 
* A.10.1.a          THEN
*                     SHIFT existing Heights right once to make
*                           room for new height
*                     INSERT new height into 1st cell of array
*/
	  if (ht_index == -1)	/* Insert before 1st cell */
	   {
	        /* shift all existing ht right once to make
		   room for new height */
	        for (i= (int)Lptr->numheights; i >= 1 ; i--)
			Lptr->height[i]= Lptr->height[i-1];
		Lptr->height[0]= pds.usHeight1;
	   }
/*
* A.10.1.b         ELSE
*                    SHIFT all cells after Ht_Index right once
*                          to make room for new height
*                    INSERT new height into Ht_Index cell of array
*/
       	  else 			/* Insert AFTER 'Ht_Index' */
	   { 
		/* shift all cells from Ht_Index+1 right once
		   to make room for new height */
	        for (i= (int)Lptr->numheights; i > ht_index; i--)
			Lptr->height[i]= Lptr->height[i-1];
                Lptr->height[ht_index+1]=pds.usHeight1;
/*
* A.10.1           ENDIF
*/
	   } 

/*
* A.10.2        INCREMENT Height counter by one
*/
	  Lptr->numheights += 1;	/* Increment Height counter */

	   DPRINT3 ("%s: Added new Height[%d] => %d, list=[", 
	   func, Lptr->numheights, pds.usHeight1);
	   for (i=0;i<Lptr->numheights; i++) 
	   DPRINT1(" %d;",(*lvl_head)->height[i]);
	   DPRINT0("]\n");
/*
* A.10         ENDIF !new height
*/
	}

/*--- now accept new Level 
  --- For Vis5d, we should only go into this loop Once (on the very 1st
  --- Vis5D message) because Vis5D option only allows 1 level type.
*/
/*
*
* A.11         IF (this is a new Level)
*              THEN
*/
    if (newlvl) {
/*
* A.11.1          SET level head pointer to Level cell
*/
	  *lvl_head = Lptr;

	  DPRINT0(" +L"); /*Insert new Level;\n"); */
	  DPRINT3("%s: Added new Level=%d, Height[0]=%d\n",
	  func,Lptr->usLevel_id, Lptr->height[0]);
/*
* A.11         ENDIF !new Level
*/
	}



/*
*
* A.9          IF (this is a new Parameter type)
* A.9.a        THEN
*/
/*--- now insert new Parm ---*/
   if (newparm) {
	  /*DPRINT0("Linking new parm;\n");*/
/*
* A.9.a.1           IF (Insert at the beginning)
* A.9.a.1.a         THEN
*                     LINK temp parm cell to beginning of linked list
*                     UPDATE Head and Tail pointer as needed
*/
          if (insrt_parm ==NULL)   	/* insrt at Beginnning */
		{
		  TmpParm->last= NULL;
		  TmpParm->next= *parm_head;
		  if (*parm_head!=NULL) (*parm_head)->last=TmpParm;
		  *parm_head=TmpParm;
          	  if (*parm_tail==NULL) *parm_tail= TmpParm; 

		}
/*
* A.9.a.1.b         ELSE
*                     LINK temp parm cell to linked list
*                     UPDATE Head and tail pointer as needed
*/
           else {                       /*-- Insert AFTER Insrt_parm;*/
                  if (insrt_parm->next != NULL) /*-- betw/ 2 cells */
                        {               
                        insrt_parm->next->last= TmpParm;
                        TmpParm->next = insrt_parm->next;
                        }
                  else  {               	/*-- at End of list */
                        TmpParm->next = NULL;
                        *parm_tail = TmpParm;
                        }
                  insrt_parm->next = TmpParm;
                  TmpParm->last = insrt_parm;
/*
* A.9.a.1.b         ENDIF
*/
                }
/*
* A.9.a.2        INCREMENT parameter counter
*/
    	   DPRINT2("%s: <<<%s>>>", func,TmpParm->abbrv);
	   *tot_parms += 1;
		   
    	   DPRINT4 ("%s: Added new Parm[%d]: '%s' (height=%d)\n",
	   func, *tot_parms, TmpParm->abbrv, pds.usHeight1);

	}  /* New Parm */
    else DPRINT1(" '%s'", TmpParm->abbrv);


/*
 --- Keeping this msg, store Info into either an existing cell (recycle mode)
 --- or a brand new cell which then needs to be inserted somewhere in
 --- the Msg Linked list .
*/

/*
* A.12         IF not recycling an old cell in List list THEN
*                  ALLOCATE storage for new cell for Msg List
*                  IF (error) RETURN(1)    !error
*              ENDIF
*/
   if (! Recycling_cell && !(NewMsg= (V5_REC *) malloc(sizeof(V5_REC))))
      { sprintf(errmsg,"ld_v5d_msg: Failed to malloc Grad_cell\n" );
  	stat=(1); goto BYE;
      }

/*
* A.13         FILL cell with info from current msg
*              !leaving it unlinked to Msg List for now
*/
      NewMsg->parm_ptr	   = TmpParm;	/* from  ABOVE */
      NewMsg->base_dtg 	   = (long)base_dtg;
      NewMsg->usP1 	   = pds.usP1;
      NewMsg->valid_etime    = (double)valid_etime; 
      NewMsg->usHeight1    = pds.usHeight1;	/* actual ht */
      NewMsg->offset       = GRIB_offs;


/* -- ONLY ENTER IF USING BRAND NEW CELL & NEED TO INSERT INTO LINK LIST */
/*
* A.14         IF not recycling an old cell (using a new Cell) 
*              THEN
*                 !must now link this new cell to Linked List
*/
    if ( ! Recycling_cell )    
    {
      NewMsg->next         = (V5_REC*)NULL; /* default */
      NewMsg->last         = (V5_REC*)NULL; /* default */

   /* .................................................................
     SORT MSGS BY BASE_TIME, TAU, PARMID, LVLID 
   .................................................................*/
/*
* A.14.1           IF (Msg List is empty)
* A.14.1.a         THEN
*                    SET head and tail of list point to new cell
*/
        if (*msgs_head==NULL) 	/* 1st msg */
	 {
	    *msgs_head= NewMsg; 
	    *msgs_tail= NewMsg; 
	  }
/*
* A.14.1.b         ELSE   !sort msgs by Basetime, Tau, Parmid, Lvlid
*/
	else 
	 {
/*
* A.14.1.b.1          FOR (each cell in Msg List)
*                   DO
*/
	    for (msg_ptr=(*msgs_head); msg_ptr!=NULL; msg_ptr=msg_ptr->next) 
	    {
/*
* A.14.1.b.1.1           IF (cell's VALID time is earlier than new msg's)
*                        THEN  keep looping
*                        ELSE  if (cell's basetime is bigger than newmsg's)
*                        THEN  break;
*/
	       	if (msg_ptr->valid_etime > NewMsg->valid_etime ) break;
	        if (msg_ptr->valid_etime < NewMsg->valid_etime) continue;

/*
* A.14.1.b.1.3           IF (cell's Parm_id is earlier than new msg's)
*                      THEN  keep looping
*                      ELSE  if (cell's Parm_id is bigger than newmsg's)
*                      THEN  break;
*/
	       	if (msg_ptr->parm_ptr->usParm_id < 
			NewMsg->parm_ptr->usParm_id) continue;
	       	else if (msg_ptr->parm_ptr->usParm_id > 
			NewMsg->parm_ptr->usParm_id) break;

/*
* A.14.1.b.1.5           IF (cell's Height is earlier than new msg's)
*                      THEN  keep looping
*                      ELSE  if (cell's Height is bigger than newmsg's)
*                            THEN  break
*                                  ELSE  PRINT message
*                                  RETURN 2 !drop duplicate message
*                            ENDIF
*/
		if (msg_ptr->usHeight1 < NewMsg->usHeight1) continue;
		else if (msg_ptr->usHeight1 > NewMsg->usHeight1) break;
		else {
		     fprintf(stdout,
		     "ld_v5d_msg:  duplicate msg not caught previously, "
		     "Dropping current msg anyway;\n");
		     free (NewMsg);
	    	     goto BYE;
		     }
/*
* A.14.1.b.1          ENDFOR !traverse Msg List
*/
	    } /* traverse Msg list */

/*
* A.14.1.b.2          IF (inserting at end)
* A.14.1.b.2.a        THEN
*                      LINK msg cell to Msg list
*                      UPDATE Head and Tail as needed
*/
	    if (msg_ptr==NULL) /* insert at end */
	     {
	       NewMsg->last 	= *msgs_tail;
   	       (*msgs_tail)->next	= NewMsg;
   	       *msgs_tail	= NewMsg;	
	     }
/*
* A.14.1.b.2.b        ELSE
*                      LINK msg to Msg List before msg_ptr
*                      UPDATE head and tail as needed
*/
	    else  	/* insert BEFORE msg_ptr, betw/ Head and Tail */
	     {
		NewMsg->next = msg_ptr;
		NewMsg->last = msg_ptr->last;
		if (msg_ptr->last!=NULL) msg_ptr->last->next= NewMsg;
	       	msg_ptr->last = NewMsg;
	        if (msg_ptr == *msgs_head) *msgs_head = NewMsg;
/*
* A.14.1.b.2.b        ENDIF
*/
	     }
/*
* A.14.1.b         ENDIF
*/
	  }
/*
* A.14         ENDIF
*/
    }  /* Link new Msg into Linked List */

      /*DPRINT4 ( "%s Prm=%d, tau=%d, ht1=%d, \n",
        NewMsg->parm_ptr->abbrv, NewMsg->parm_ptr->usParm_id,
        NewMsg->ustau);
	DPRINT7(
	"Dpos=%ld, Bpos=%ld, %d bits, D=%.1lf, BinScl=%.2lf, Ref=%.2lf\n\n",
	NewMsg->usHeight1, NewMsg->dpos, NewMsg->bpos, NewMsg->bnum,
        NewMsg->fDec_sc_fctr, NewMsg->fBin_sc_fctr, NewMsg->fReference);
	*/

/*
* 
* A.17           INCREMENT message counter
*/
  *tot_msgs+= 1;
  DPRINT1("Bump totmsgs to %ld\n", *tot_msgs);

BYE:
/*
* 
* A.18           RETURN stat
*/
  DPRINT1("Leaving ld_v5d_msg, Stat=%d\n", stat);
  return(stat);
/*
*
* END OF FUNCTION
*/
#endif
}
