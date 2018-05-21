/*........................................................
Filename:  ld_grad_msg.c
Author  :  Alice Nakajima, SAIC
Date    :  07/17/96

12/05/96 atn:  -Stderr, +errmsg, -(table2 *parmtab), +(PARM_DEFN);
06/26/97 atn:  +prototypes;
07/07/97 atn:  fixed order of GMP's BLK3&4 offset; altered GRADS script file to
	only loop for ZDEF level;
08/30/98 atn:  
	+ comments to ctl file and script file;
	* let GRAD_PARM's usParm_id be the combo of Parmid & ParmSub to 
	  support grib extensions.  
	* Abbrv now has format:
	 'aPPPLLL[a-z]' if parmid <= 999
	 'hxxxLLL[a-z]' if parmid >= 1000 (xxx is ID in Hex)
*/
#include <stdio.h>
#include <math.h>
#include "grib_lookup.h"  /* for all lookup structs */
#include "grads.h"	  /* for all grads structs */
#include "isdb.h"  	  /* struct date */
#include "dprints.h"	/* for dprints */
#include "gribfuncs.h"	/* prototypes */

static int  unit_offset[5]= 
 { 1, 100, 10000, 1000000, 100000000 }; /* hr, month, day, year, cent */
/*
*
* ======================================================================
* A.  FUNCTION:  ld_grad_msg
*
*     PURPOSE: 
*     Called on to store info of the Msg just read into the arrays
*     that will be used to create grads control file & mapping file;
*     - does not allow duplicates of messages;
*     - gives preference to msg of certain parm&Lvl&Ht of Earlier DTGs
*
*     INPUT:
*        int        *grad_full;	  set if all out of extensions 
*        long       GRIB_offs;    #bytes fr. beginning of file
*        PDS_INPUT  pds;          Product descr section
*        grid_desc_sec  gds;      Grid descr section
*        BDS_HEAD_INPUT bds_head; Binary Data section
*        BMS_INPUT  bms;          Bitmap section
*        int        *tot_parms;   count of unique parameter types
*        int        *tot_levels;  count of unique level types
*        int        *tot_msgs;    count of unique Messages
*        GRAD_PARM  **parm_head;  head of Parameter linked list
*        GRAD_PARM  **parm_tail;  tail of Parameter linked list
*        GRAD_REC   **msgs_head;  head of Grib Msg linked list
*        GRAD_REC   **msgs_tail;  tail of Grib Msg linked list
*        GRAD_LVL   **lvl_head;   head of Level linked list
*        GRAD_LVL   **lvl_tail;   tail of Level linked list
*        GRAD_INFO  *grad_info;   holds common info
*        char	    *errmsg;      holds err message if occurs
*
*     OUTPUT:
*         If no errors: 
*          -message is inserted in Msg linked list sorted by Base time,
*           Forecast Period, Parameter and Leveltype;
*          -if new Level type then add it to Level Linked list sorted
*           by level_id;
*          -if new Height then add it to Height array of this Levelid,
*           sorted in ascending order;
*          -if new parameter type then add it to Parameter Linked list
*           sorted by parm_id;
*
*     RETURN CODE:
*         0:  no errors, may/may not have stored this msg;
*             (Duplicate message, out of Exntensions, invalid forecast unitid)
*         1:  Malloc error;
* ======================================================================
*/
extern int 		UseTables;
extern PARM_DEFN        db_parm_tbl[];

#if PROTOTYPE_NEEDED
   int  ld_grad_msg (
		int 		*grad_full,
		long 		GRIB_offs,
		PDS_INPUT 	pds, 
		grid_desc_sec 	gds, 
		BDS_HEAD_INPUT 	bds_head,
		BMS_INPUT 	bms,
		int 		*tot_parms, 
		int 		*tot_levels, 
		int 		*tot_msgs,
		GRAD_PARM 	**parm_head, 
		GRAD_PARM 	**parm_tail,
		GRAD_REC  	**msgs_head, 
		GRAD_REC 	**msgs_tail,
		GRAD_LVL  	**lvl_head, 
		GRAD_LVL 	**lvl_tail,
		GRAD_INFO 	*grad_info,
		char 		*errmsg)
#else
   int  ld_grad_msg ( grad_full, GRIB_offs, pds, gds, bds_head, bms,
		tot_parms, tot_levels, tot_msgs, parm_head, parm_tail, 
		msgs_head, msgs_tail, lvl_head, lvl_tail, grad_info, errmsg)

		int 		*grad_full;
		long 		GRIB_offs;
		PDS_INPUT 	pds; 
		grid_desc_sec 	gds; 
		BDS_HEAD_INPUT 	bds_head;
		BMS_INPUT 	bms;
		int 		*tot_parms, *tot_levels, *tot_msgs;
		GRAD_PARM 	**parm_head, **parm_tail;
		GRAD_REC  	**msgs_head, **msgs_tail;
		GRAD_LVL  	**lvl_head, **lvl_tail;
		GRAD_INFO 	*grad_info;
		char 		*errmsg;
#endif
{
char  	pidStr[5];	/* string representation of the parmid+parmsub */
char 	abbrv[10];	/* work var */
int	Parmtbl_indx;   /* within the Parameter Conversion Array */
int	lvl_indx, n,i,j,x, unit_incr;
int     newparm=1;	/* set if need to add new parm */
int	newlvl=0;	/* set if need to add new Level */
int	newheight=0;	/* set if need to add new height */
int	ht_index=-1;	/* indx within lvl's ht array */
int     tau_incr; 
int     stat=0;
int	ipid;		/* parmtbl_indx(parmid, parmsub) value */
long	base_dtg;	/* work var */
GRAD_REC   *msg_ptr; 	/* temp pointer */
GRAD_REC   *NewMsg;      /* cell for current message */
GRAD_PARM  *parm_ptr;    /* temp pointer */
GRAD_PARM  *insrt_parm;  /* working ptr */
GRAD_PARM  *TmpParm;	 /* is Parm Cell that matched/got created */
GRAD_LVL   *insrt_lvl;	/* marks where to insert cell */
GRAD_LVL   *Lptr;        /*  is Lvl Cell that matched/got created */
DATE    date; 
double  hour;
double  etime;   /* epochal time */
void	h_e_time();

/*
*
* A.1          INIT variables
*/
      DPRINT6 ("Inside ld_grad_msg\nLoad: %02d%02d%02d%02d%02d fcst=%02d ",
                  pds.usCentury-1, pds.usYear,pds.usMonth, pds.usDay, 
                 pds.usHour, pds.usP1);

      DPRINT3 ("P=%03d L=%03d Ht=%03d;", 
		pds.usParm_id, pds.usLevel_id, pds.usHeight1);

     /*..........................................................
        Maintain a Level linked list where level ids are sorted
	in ascending order;
        if usLevel_id is new:  add new lvl type;
        if usLevel_id exists:  only add ht to height[] array if it
                           hasn't already been defined;
     ..........................................................*/
      insrt_lvl	= NULL; 	/* default */
      Lptr	= *lvl_head;	/* default */

/*
*
* A.2          WHILE (more Level cell in Level list AND
*                     msg's level_id is greater than cell's levelid)
*              DO
*                 MARK insrt_lvl as cell to insert before
*                 LOOP to check next cell
*              ENDWHILE
*/
      /* loop until (usLevel_id > curr usLevel_id) or til NULL */
      while (Lptr!=NULL && pds.usLevel_id > Lptr->usLevel_id)
	{ insrt_lvl= Lptr;  Lptr= Lptr->next; }
	/**** 	
	DONOT CHANGE 'INSRT_LVL' AFTER THIS POINT;
	***/

/*
*
* A.3          IF (level_id already exists) 
*/
      if (Lptr != NULL && Lptr->usLevel_id == pds.usLevel_id) 
	{    
	/* LEVL EXISTS, ADD HEIGHT IF NEW  
	-- Lptr is the cell that matched
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
      else 
         {    /* NO MATCH, ADD NEW LEVEL (& HEIGHT)
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
           if ((Lptr=(GRAD_LVL *)calloc(sizeof(GRAD_LVL),1))==NULL) { 
		 sprintf(errmsg,"ld_grad_msg: Failed to malloc new level \n");
                 stat=(1); goto BYE;
		}

/*
* A.3.b.3         FILL Temp level cell  !leave unlinked to List for now
*/
           Lptr->usLevel_id     = pds.usLevel_id;
           Lptr->numheights     = 1;
           Lptr->height[0]      = pds.usHeight1;
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
   insrt_parm=NULL;	/* default */
   TmpParm = NULL;	/* default */
   newparm = 1;		/* default */
    
   /* 	Grads variables are only allowed 8 characters long.
	format:  "xPPPLLL[a]", 
	where  
		x= m/a/b/c/d/e (Main or SubTable A-E);
		PPP= 3 digit parm id if Main table,  or
		     3 digit parm Sub id if Sub table A-E;
		LLL= 3 digit level id;
		[a]= optional char to represent non-zdef variables.
   */
   ipid =  PARMTBL_INDX(pds.usParm_id, pds.usParm_sub);
   if (ipid < 256) 
   	sprintf(pidStr , "m%03d", pds.usParm_id);
   else sprintf(pidStr , "%c%03d", 'a' + (pds.usParm_id - 250), pds.usParm_sub);

   if (pds.usLevel_id==grad_info->zdef_lvl 
       		|| grad_info->zdef_lvl==999)	/* default */ 
       	 	sprintf (abbrv, "%s%03d", pidStr, pds.usLevel_id);
   	else 	sprintf (abbrv, "%s%03da", pidStr, pds.usLevel_id);
    
/*
* A.5          FOR (each cell in Parameter Linked List) DO
*/
     for (parm_ptr=*parm_head;  parm_ptr!=NULL;  parm_ptr=parm_ptr->next)
     { /* check if Parm already exists */
/*
* A.5.1             IF (msg's lvl_id and parm_id differ from 
*                      that of current cell in list)
*                   THEN
*                      MARK insrt_parm as place to insert after
*                      LOOP, check next cell
*                   ENDIF
*/
	if (pds.usLevel_id != parm_ptr->usLevel_id ||
	    PARMTBL_INDX(pds.usParm_id, pds.usParm_sub) != parm_ptr->usParm_id)
	  {	insrt_parm=parm_ptr;	/* Insert AFTER Insrt_parm */
		continue;
	  }
/*
* A.5.2             IF (this is 1st msg    !zdef_lvl still 999
*                       OR  msg's level is the ZDEF level
*                       OR  msg's Height matches that of current non-zdef level)
*                   THEN   
*                      SET newparm to 0     !this parm already exists
*                      COPY abbrv name to fill new cell later
*                      BREAK out of loop
*                   ENDIF
*/
	if (grad_info->zdef_lvl==999 ||
	    pds.usLevel_id == grad_info->zdef_lvl ||
	    pds.usHeight1 == parm_ptr->usHeight )
	  {
		newparm=0;
		strcpy (abbrv, parm_ptr->abbrv);
		break;
	  }


	/* GET HERE IF:  
		ZDEF has already been defined (zdef != 999)
		curr level is not ZDEF level  (&& lvlid!=zdefid )
		curr Height has not been defined (&& ht != nonzht)
	*/
/*
* A.5.3             IF (levelid is not the ZDEF level) THEN
*                      COPY curr cell's abbrv to variable abbrv
*                      IF (extension used is 'z') 
*                      THEN bump new extension to 'A'
*                      ELSE if (extension used is 'Z')
*                      THEN 
*                          SET grad_full flag to 1   !no more room
*                          PRINT message
*                          RETURN 0 !drop msg but no errors
*                      ELSE BUMP extension up by 1
*
*                      INCREMENT the extension by 1
*                      MARK insrt_parm as place to insert after
* A.5.3                ENDIF
*/
	if (pds.usLevel_id != grad_info->zdef_lvl)
	   if (strcmp (abbrv, parm_ptr->abbrv) <= 0) 
	   {
		strcpy (abbrv, parm_ptr->abbrv);
		/* bump extension up for new variable  */
		if (abbrv[7] == 'z') abbrv[7]= 'A';
		else if (abbrv[7] == 'Z') {
		   *grad_full= 1;	/* set global flag */
		   fprintf (stdout,
   "ld_grad_msg warning: out of extensions for GrRADS VarName, dropmsg\n");
		   goto BYE;
		   }
		else abbrv[7]= abbrv[7]+1;	
	   	insrt_parm=parm_ptr;	/* Insert AFTER Insrt_parm */
	   }
	   
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
* A.6.a.1        ALLOCATE storage for the new parm
*                IF (error) RETURN 1
*/
          if ((TmpParm=(GRAD_PARM *)calloc(sizeof(GRAD_PARM),1))==NULL)
          { sprintf(errmsg,"ld_grad_msg: Failed to calloc GRAD_PARM\n");
	    stat=(1); goto BYE; }
/*
* A.6.a.2        INCREMENT parameter counter
*/
	  *tot_parms += 1;
/*
* A.6.a.3        FILL this temp Parameter cell
*                !leaving it unlinked to Parameter list for now
*/
	  strcpy (TmpParm->abbrv, abbrv);
          TmpParm->usLevel_id = pds.usLevel_id;
          TmpParm->usParm_id  = PARMTBL_INDX(pds.usParm_id, pds.usParm_sub);
	  TmpParm->usHeight   = pds.usHeight1;
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
* A.8          IF (the common Time Unit hasn't been defined)
*              THEN
*                 LET tau_incr be 1st msg's forecst time
*/
  if (grad_info->ebase_dtg == 999)
   {
	/* Undefined time unit means this is 1st msg; */
	tau_incr= (int)pds.usP1;
   }
/*
* A.8.b        ELSE
*                !check for duplicates (see if already have msg of same
*                !time from an earlier or later or same dtg)
*/
  else
   {	
/*
* A.8.b.1        SET up base time of curr message
*/
      date.year  = (int) ((pds.usCentury-1)*100+pds.usYear);
      date.month = (int) pds.usMonth;
      date.day   = (int) pds.usDay;
      if (pds.usSecond != 0 && pds.usSecond <= 60)
           hour = (double)pds.usHour + (double)pds.usMinute/60.
		   + (double)pds.usSecond/3600.;
      else hour	= (double)pds.usHour + (double)pds.usMinute/60.;

/*
* A.8.b.2        COMPUTE time offset based on forecast unit id
*/
      switch  (pds.usFcst_unit_id)
	{
	case 0: /* in minutes */ hour += (double)pds.usP1 / 60.; break;
	case 1: /* in hours */ hour += (double)pds.usP1; break;
	case 2: /* in days */ date.day += (int)pds.usP1; break;
	case 3: /* in months */ date.month += (int)pds.usP1; break;
	case 4: /* in years */ date.year += (int)pds.usP1; break;
	case 5: /* in decades (10 years) */ 
		date.year += (10 * (int)pds.usP1); break;
	case 6: /* in normal (30 years) */ 
		date.year += (30 * (int)pds.usP1); break;
	case 7: /* in century (100 years) */ 
		date.year += (100 * (int)pds.usP1); break;
	case 254: /* in seconds */
		hour += ((float)pds.usP1 / 3600.); break;
	default: fprintf (stdout,
	       "ld_grad_msg warning: Invalid Forecast Unit Id (%d), drop msg\n"
		,pds.usFcst_unit_id);
		goto BYE;
	}

/*
* A.8.b.3        FUNCTION h_e_time  !convert to Epochal time
*/
      h_e_time (&date, &hour, &etime);

/*
* A.8.b.4        LET tau_incr be difference between this msg's basetime
*                    and the common basetime
*/
      tau_incr = (int)(etime - grad_info->ebase_dtg);	/* in hours so far */ 
      
/*
* A.8.b.5        FOR (each cell in the Message list) DO
*/
      for (msg_ptr=*msgs_head; msg_ptr!= NULL; )
      {
/*
* A.8.b.5.1        IF (newmsg's tau_incr differs from cell's tau_incr
*                     OR if their abbrv differ
*                     OR if their Height level differ)
*                  THEN
*                     GO check the next cell in msg list
*                  ENDIF
*/
         if (tau_incr!=msg_ptr->tau_incr ||
	     strcmp (abbrv, msg_ptr->parm_ptr->abbrv) ||
	     pds.usHeight1 != msg_ptr->usHeight1) 
	  {  
		msg_ptr=msg_ptr->next; continue;
	   }

         /* Same TauIncr, Same Abbrv so far: 
	    Now choose the earlier DTG...  
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
 "ld_grad_msg:  drop Duplicate msg [%d tau=%d '%s' Pid=%d Sub=%d Lvl=%d @%d]\n",
	    base_dtg, pds.usP1, abbrv,pds.usParm_id,pds.usParm_sub,
	    pds.usLevel_id, pds.usHeight1);
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
	    DPRINT2(
	    " Remove old '%ld' incr=%d\n",msg_ptr->base_dtg,msg_ptr->tau_incr);
	   /*
            "\nCurrent Msg [DTG=%ld, FcstPer=%d (TauIncr=%d), '%s' ht=%d] "
	    "has smaller DTG than:\n"
            "stored Msg [DTG=%ld, FcstPer=%d (TauIncr=%d), '%s' ht=%d] \n"
	    "REMOVING the Previously Stored Msg & keeping current one\n"
		    ,base_dtg, pds.usP1, tau_incr, abbrv, pds.usHeight1,
		    msg_ptr->base_dtg,msg_ptr->ustau,msg_ptr->tau_incr, 
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
* A.8.b.5.2.b.3       FREE its storage
*/
            free (msg_ptr);	/* remove cell */
/*
* A.8.b.5.2.b.4       BREAK !all done
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
	    DPRINT0 (" bigger dtg, Ignored\n");

	    fprintf(stdout,
            "\nCurrent Msg [DTG=%ld, FcstPer=%d (TauIncr=%d), '%s' ht=%d]\n"
	    "has bigger DTG than previously stored message;\n"
            "Previous Msg [DTG=%ld, FcstPer=%d (TauIncr=%d), '%s' Ht=%d] "
	    "\nGrads not storing current message...\n",
	    base_dtg, pds.usP1, tau_incr, abbrv, pds.usHeight1,
	    msg_ptr->base_dtg,msg_ptr->ustau,msg_ptr->tau_incr, 
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
  }
	    
/*    
*              ! Not a duplicate msg, ok to keep it;
*/

/*
*
* A.9          IF (this is a new Parameter type)
*              THEN
*/
/*--- now insert new Parm ---*/
   if (newparm) {
	  /*DPRINT0("Linking new parm;\n");*/
/*
* A.9.1           IF (Insert at the beginning)
* A.9.1.a         THEN
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
* A.9.1.b         ELSE
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
* A.9.1.b         ENDIF
*/
                }
    	   DPRINT1(" <<<%s>>>", TmpParm->abbrv);
/*
* A.9          ENDIF
*/
	}
    else DPRINT1(" '%s'", TmpParm->abbrv);

/*--- new insert new height */
/*
*
* A.10         IF (this is a new Height)
*              THEN
*/
    if (newheight) {
	  DPRINT0 (" +Ht");
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
/*
* A.10         ENDIF !new height
*/
	}

/*--- now insert new level */
/*
*
* A.11         IF (this is a new Level)
*              THEN
*/
    if (newlvl) {
	   DPRINT0(" +L"); /*Insert new Level;\n"); */
/*
* A.11.1          IF (Insert at beginning of list)
* A.11.1.a        THEN
*                    INSERT new level cell to begining of Lvl list
*                    UPDATE Head and Tail pointers as needed
*/
	   if (insrt_lvl == NULL)  	/*-- Insert at beg. of list  */
		{
		  Lptr->last	= NULL;
		  Lptr->next	= *lvl_head;
		  if (*lvl_head!= NULL) (*lvl_head)->last= *lvl_head;
		  *lvl_head		= Lptr;
		  if (*lvl_tail==NULL) *lvl_tail= Lptr;
		}
/*
* A.11.1.b        ELSE
*                    INSERT new level cell after Insrt_lvl within list
*                    UPDATE Head and Tail pointers as needed
*/
	   else {			/*-- Insert AFTER Insrt_Lvl;*/

		  if (insrt_lvl->next != NULL) 
			{		/*-- Insert betw/ 2 cells */
			insrt_lvl->next->last= Lptr;
			Lptr->next	= insrt_lvl->next;
			}
		  else  { 		/*-- Insrt at End of list */
			Lptr->next	= NULL;
			*lvl_tail	= Lptr;
			}
		  insrt_lvl->next	= Lptr;
		  Lptr->last	= insrt_lvl;
/*
* A.11.1.b        ENDIF
*/
		}
/*
* A.11.2          INCREMENT level counter by one
*/
	  *tot_levels += 1;
/*
* A.11         ENDIF !new Level
*/
	}

   /*.........................................................
        store Msg into linked list of messages
	 ABBRV		: set later after knowing Extension;
	 BASE_DTG       : base time of Msg (yyyymmddhh);
         USTAU          : forecast period of Msg;
         USHEIGHT1      :  top level value of Msg;
         BNUM           : number of bits to store values into;
         FDEC_SC_FCTR   : 10**D
         FBIN_SC_FCTR   : 2**SF
         FREFERENCE     : reference val
         BPOS           : byte offs fr beg. of file to where
                          bitstr of BMS begins;
         DPOS           : byte offs fr beg. of file to where
                          bitstr of BDS begins;
    ..........................................................*/

/*
* A.12         ALLOCATE storage for new cell for Msg List
*              IF (error) RETURN(1)    !error
*/
      NewMsg= (GRAD_REC *) calloc(sizeof(GRAD_REC),1);
      if (NewMsg==(GRAD_REC *)NULL)
      { sprintf(errmsg,"ld_grad_msg: Failed to malloc Grad_cell\n" );
  	stat=(1); goto BYE;
      }
/*
* A.13         FILL cell with info from current msg
*              !leaving it unlinked to Msg List for now
*/
      NewMsg->parm_ptr	   = TmpParm;	/* from  ABOVE */
      NewMsg->base_dtg 	   = (long)base_dtg;
      NewMsg->tau_incr	   = (long)tau_incr;
      NewMsg->ustau        = pds.usP1;
      NewMsg->usHeight1    = pds.usHeight1;	/* actual ht */
      NewMsg->dpos         = (bds_head.length <= 11) ? -999:
			     GRIB_offs+sizeof(IDS_GRIB)+pds.uslength
			     + gds.head.uslength+ bms.uslength + 11;
      NewMsg->bpos         = (bms.uslength <= 6) ?  -999:
			     GRIB_offs+ sizeof(IDS_GRIB)+ pds.uslength
			     + gds.head.uslength+6;
      NewMsg->bnum         = bds_head.usBit_pack_num;
      NewMsg->fDec_sc_fctr = (float)pow(10.,(int)pds.sDec_sc_fctr);
      NewMsg->fBin_sc_fctr = (float)pow(2.,(int)bds_head.Bin_sc_fctr);
      NewMsg->fReference   = (float)bds_head.fReference;
      NewMsg->next         = (GRAD_REC*)NULL; /* default */
      NewMsg->last         = (GRAD_REC*)NULL; /* default */

   /* .................................................................
     SORT MSGS BY BASE_TIME, TAU, PARMID, LVLID 
   .................................................................*/
/*
* A.14           IF (Msg List is empty)
* A.14.a         THEN
*                    SET head and tail of list point to new cell
*/
        if (*msgs_head==NULL) 	/* 1st msg */
	 {
	    *msgs_head= NewMsg; 
	    *msgs_tail= NewMsg; 
	  }
/*
* A.14.b         ELSE   !sort msgs by Basetime, Tau, Parmid, Lvlid
*/
	else 
	 {
/*
* A.14.b.1          FOR (each cell in Msg List)
*                   DO
*/
	    for (msg_ptr=*msgs_head; msg_ptr!=NULL; msg_ptr=msg_ptr->next) 
	    {
/*
* A.14.b.1.1           IF (cell's basetime is earlier than new msg's)
*                      THEN  keep looping
*                      ELSE  if (cell's basetime is bigger than newmsg's)
*                      THEN  break;
*/
	        if (msg_ptr->base_dtg < NewMsg->base_dtg) continue;
	       	else if (msg_ptr->base_dtg > NewMsg->base_dtg ) break;

/*
* A.14.b.1.2           IF (cell's forecast time is earlier than new msg's)
*                      THEN  keep looping
*                      ELSE  if (cell's forecast time is bigger than newmsg's)
*                      THEN  break;
*/
	       	if (msg_ptr->ustau < NewMsg->ustau) continue;
		else if (msg_ptr->ustau > NewMsg->ustau) break;

/*
* A.14.b.1.3           IF (cell's Parm_id is earlier than new msg's)
*                      THEN  keep looping
*                      ELSE  if (cell's Parm_id is bigger than newmsg's)
*                      THEN  break;
*/
	       	if (msg_ptr->parm_ptr->usParm_id < 
			NewMsg->parm_ptr->usParm_id) continue;
	       	else if (msg_ptr->parm_ptr->usParm_id > 
			NewMsg->parm_ptr->usParm_id) break;

/*
* A.14.b.1.4           IF (cell's Level_id is earlier than new msg's)
*                      THEN  keep looping
*                      ELSE  if (cell's Level_id is bigger than newmsg's)
*                      THEN  break;
*/
	       	if (msg_ptr->parm_ptr->usLevel_id < 
			NewMsg->parm_ptr->usLevel_id) continue;
	       	else if (msg_ptr->parm_ptr->usLevel_id > 
			NewMsg->parm_ptr->usLevel_id) break;


/*
* A.14.b.1.5           IF (cell's Height is earlier than new msg's)
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
		     "ld_grad_msg:  duplicate msg not caught previously, "
		     "Dropping current msg anyway;\n");
		     free (NewMsg);
	    	     goto BYE;
		     }
/*
* A.14.b.1          ENDFOR !traverse Msg List
*/
	    } /* traverse Msg list */

/*
* A.14.b.2          IF (inserting at end)
* A.14.b.2.a        THEN
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
* A.14.b.2.b        ELSE
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
* A.14.b.2.b        ENDIF
*/
	     }
/*
* A.14.b         ENDIF
*/
	  }

    DPRINT1 (" incr=%d\n", NewMsg->tau_incr);

      /*DPRINT4 ( "%s Prm=%d, lvl=%d, tau=%d, ht1=%d, \n",
        NewMsg->parm_ptr->abbrv, NewMsg->parm_ptr->usParm_id,
        NewMsg->parm_ptr->usLevel_id, NewMsg->ustau);
	DPRINT7(
	"Dpos=%ld, Bpos=%ld, %d bits, D=%.1lf, BinScl=%.2lf, Ref=%.2lf\n\n",
	NewMsg->usHeight1, NewMsg->dpos, NewMsg->bpos, NewMsg->bnum,
        NewMsg->fDec_sc_fctr, NewMsg->fBin_sc_fctr, NewMsg->fReference);
	*/

/*
*
* A.15           IF (common basetime is undefined)
*                THEN
*                    SET the "common" DTG to this msg's dtg
*                    SET the ZDEF level to this msg's level
*                    DEBUG print
*/
    if (grad_info->base_dtg==999)
    {
      /* Since already accepted this msg, allow its characteristics
	 to be the "common" info if "common" info has not yet been
	 defined;
	*/
	grad_info->zdef_lvl = pds.usLevel_id;
        grad_info->base_dtg = base_dtg;
        date.year  = (int) ((pds.usCentury-1)*100+pds.usYear);
        date.month = (int) pds.usMonth;
        date.day   = (int) pds.usDay;
        hour = (double)pds.usHour;
        h_e_time (&date, &hour, &grad_info->ebase_dtg);

        DPRINT3("SETTING Common Base=%d, Zdef_Level=%d, Ebase=%.4lf\n", 
	grad_info->base_dtg, grad_info->zdef_lvl, grad_info->ebase_dtg);
/*
* A.15           ENDIF
*/
    }
  DPRINT1 ("grad_info->base_dtg is %ld\n", grad_info->base_dtg);

/*
*
* A.16           IF (common Grid Size is undefined)
*                THEN
*                    SET them equal to this msg's info
*                    DEBUG print
*/
    if (grad_info->ulGrid_size==999)
	grad_info->ulGrid_size= (long)bds_head.ulGrid_size;
    DPRINT1 ("grad_info->ulGrid_size is %ld\n",grad_info->ulGrid_size);

/*
*
* A.16           IF (this is the ZDEF level)
*                THEN  
*                    MARK attribute usHeight as not used (65535)
*                    !ZDEF level's Heights are actually stored in the
*                    !the GRAD_LVL's Height[] array; 
*                ENDIF
*/
    if (pds.usLevel_id == grad_info->zdef_lvl) {
	TmpParm->usHeight = 65535;
        DPRINT1 ("Zdef Level, NOT using TmpParm->usHeight (set to %d)\n"
        ,TmpParm->usHeight);
	}
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
  DPRINT1("Leaving ld_grad_msg, Stat=%d\n", stat);
  return(stat);
/*
*
* END OF FUNCTION
*/
}
