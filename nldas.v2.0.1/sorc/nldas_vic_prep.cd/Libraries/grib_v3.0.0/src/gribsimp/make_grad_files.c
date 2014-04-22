/*........................................................
Filename:  make_grad_files.c
Author  :  Alice Nakajima, SAIC
Date    :  07/17/96
06/26/97 atn:  +prototypes;
07/10/97 atn:  only loop in GRADS script if level is ZDEF level;
07/02/98 atn:  free xdefn,ydefn,options,pdefn;
10/22/98 atn:  prevent reclosing of closed files;
*/
#include <stdio.h> 
#include <string.h>
#include "isdb.h"    	/* date struct  */
#include "grads.h"	/* grads structs */
#include "dprints.h"	/* for dprints */
#include "gribfuncs.h"	/* prototypes */

#define ISLEAP(year)    (!((year) % 4) && ((year) % 100) || !((year) % 400))
#define GO_HOME		goto BAIL_OUT;

extern char *InFile;	/* name of input filename */
static char *month_name[]=
{ "jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec" };

/*
*
*
* ======================================================================
* A.  FUNCTION:   make_grad_files
*
*     PURPOSE :
*       Create the two control files to be used with the GRADS software,
*       (".gmp" binary file and ".ctl" ascii control file) and also a
*       GRADS script file called "draw_all.gs" ('run draw_all.gs')
*
*     INPUT   :
*       PDS_INPUT *pds,         Production Defn Section
*       grid_desc_sec *gds,     Grid Defn Section
*       GRAD_INFO *common,   Common Grad info
*       int 	parmcnt,         count of parameters
*       GRAD_PARM **parm_head,  head of Parameter Linked list
*       GRAD_REC **msgs_head,   head of Message Linked list
*       GRAD_LVL **lvl_head)    head of Level Linked list
*
*     OUTPUT  :
*       If successfull, ./fn.ctl, ./fn.gmp and 'draw_all.gs'
*       where 'fn' is the prefix of the input Grib file
*       .ctl is an ascii file and .gmp is a binary file
*       draw_all.gs is the grads script file to display all fields
*
*     RETURN CODE:
*       0 for no errors, 
*       1 for errors; see errmsg;
* ======================================================================
*/
#if PROTOTYPE_NEEDED
    int   make_grad_files (
    PDS_INPUT 		*pds,
    grid_desc_sec 	*gds,
    GRAD_INFO 		*common,
    int 		parmcnt,
    GRAD_PARM 		**parm_head,
    GRAD_REC 		**msgs_head,
    GRAD_LVL 		**lvl_head,
    char		*errmsg)
#else
    int   make_grad_files (pds, gds, common, parmcnt, 
			   parm_head, msgs_head, lvl_head, errmsg)

    PDS_INPUT 		*pds;	        /* input */
    grid_desc_sec 	*gds;	        /* input */
    GRAD_INFO 		*common;/* input= common info shared by all messages */
    int 		parmcnt;        /* input */
    GRAD_PARM 		**parm_head; 	/* input */
    GRAD_REC 		**msgs_head; 	/* input */
    GRAD_LVL 		**lvl_head; 	/* input */
    char		*errmsg;	/* output */
#endif

{
  FILE  *fctl= (FILE *)NULL;  		/* output .ctl file */
  FILE  *fgmp= (FILE *)NULL;  		/* output .gmp file */
  FILE  *fscript= (FILE *)NULL; 	/* output script file */
  GRAD_REC *pmesg, *msg;    		/* working pointer */
  GRAD_PARM *parm_ptr;	    		/* working pointer */
  GRAD_LVL *lvl_ptr;	    		/* working pointer */
  GMP_BLK0 *blk0 = (GMP_BLK0 *)NULL;	/* for the .gmp file */
  GMP_BLK1 *blk1 = (GMP_BLK1 *)NULL;	/* for the .gmp file */
  GMP_BLK3 *blk3 = (GMP_BLK3 *)NULL;	/* for the .gmp file */
  GMP_BLK4 *blk4 = (GMP_BLK4 *)NULL; 	/* for the .gmp file */
  DATE  date;		      	        /* for Date routines  */
  char  *func="make_grad_files";
  char  *cp;		    		/* working var */
  char  temp[200];
  char  *ctl_fn= (char*)NULL;    	/* fname for .CTL file */
  char  *gmp_fn= (char*)NULL;    	/* fname for .GMP file */
  char *Options= (char*)NULL; 		/* for the .CTL file */
  char *Pdefn= (char*)NULL; 		/* for the .CTL file */
  char *Xdefn= (char*)NULL; 		/* for the .CTL file */
  char *Ydefn= (char*)NULL; 		/* for the .CTL file */
  int   i, j;		    /* working vars */
  int	newtime;	    /* set if this is new time */
  int	time;		    /* working var */
  int   Tindex;	    	    /* index for Time array */
  int	time_cnt;	    /* count of #unique Times */
  int   recspertime;	    /* number of msgs per TIme 		*/
  int   varht_offs;	    /* offs of each Variable group in each Time */
  int   offset;		    /* offs of each Msg w/in grad group */
  int   elements;	    /* max number of records for .GMP file */
  int	smallest, largest;  /* to compute Tdef line */
  int   daysperyear;	    /* depends on leap yr or not */
  int   doy1, doy2;         /* day of year, compute tdef */
  int   day_cnt;            /* to compute Tdef line */
  int   yr, yr1, yr2;       /* year, compute tdef line */
  int   Ecode = 1;	    /* default to bad exit code */
  long	Time[100]; 	    /* holds unique Tau Increments  */
  long	tau_incr=99999999;  /* Time tau_increrence between msgs */

  DPRINT1("Entering %s\n",func);
  /* Find smallest Tau Increment from base time => TAU_INCR; */
  /* store all the Increments into array Time[];*/

  DPRINT0("Unique Taus are= ");
/*
*
* A.1       FOR (each message in Msg list )
*           DO
*/
  for (Tindex=0, msg=*msgs_head; msg != NULL; msg=msg->next) 
     {
  	/* find smallest Tau increment from base time */
/*
* A.1.1        IF (msg is head of list) 
* A.1.1.a      THEN
*                  IF (this tau_incr is not zero)
*                  SET smallest increment to this tau_incr
*/
	if (msg->last==NULL) {
		if (msg->tau_incr>0) 
		tau_incr= (long) msg->tau_incr; 
	   }
/*
* A.1.1.b      ELSE
*                  IF (difference between tau_incr of current and
*                      previous message is smaller than tau_incr)
*                  SET tau_incr to this difference
*/
	else {
	        if (msg->tau_incr - msg->last->tau_incr > 0 &&
	            msg->tau_incr - msg->last->tau_incr < tau_incr)
		tau_incr= (msg->tau_incr - msg->last->tau_incr);
/*
* A.1.1.b      ENDIF
*/
	   }

/*
* A.1.2        !default to newtime flag=1
*              FOR (each stored tau_incr) DO
*                 IF (see one that matches current tau_incr)
*                  THEN
*                      SET newtime flag to 0
*                      BREAK
*                  ENDIF
*              ENDDO
*/
	for (newtime=1, i=0; i < Tindex; i++) 
	    if (msg->tau_incr == Time[i]) { newtime=0; break; }

/*
* A.1.3        IF (newtime flag is set) THEN
* A.1.3.1          DEBUG print
*/
 	if (newtime)  {
		DPRINT1 ("%d ", msg->tau_incr);
/*
* A.1.3.2          IF (exceeded max number of tau_incr) THEN 
*                      PRINT error message
*                      RETURN 1   !error
*/
		if (Tindex== 100) 
		 {
		  sprintf(errmsg,"%s: Out of room for Time array\n", func);
		  GO_HOME;
/*
* A.1.3.2          ENDIF
*/
		 }
		Time[Tindex++]=  msg->tau_incr;
/*
* A.1.3        ENDIF
*/
		} /* Newtime */
/*
* A.1       ENDFOR
*/
     } /* for */

     DPRINT1 ("; Smallest incr=%d\n", tau_incr);

/*
* 
* A.2       IF (tau_incr is undefined or zero) 
* A.2.a     THEN
*               LET tau_incr be 1
*/
  if (tau_incr== 99999999 || tau_incr==0) 
	{ tau_incr= 1;  }
/* 
* A.2.b     ELSE
*/
  else {
/*
* A.2.b.1       FOR (each message in Msg list)
*               DO
*                   IF (msg's tau_incr is divisible by tau_incr) CONTINUE;
*                   DEBUG print
*                   RESET tau_incr to 1
*                   BREAK
*/
     for (msg=(*msgs_head); msg != NULL; msg= msg->next) {
        if ((int)msg->tau_incr % (int)tau_incr == 0) continue;
	DPRINT4("\nWARNING:  DTG %ld Tau %d (incr=%d) not divisible by %d;"
	" Reset Increment to '1 hr';\n",
	msg->base_dtg, msg->ustau, msg->tau_incr, tau_incr);
	tau_incr=1;
	break;
/*
* A.2.b.1       ENDFOR
*/
	}
/* 
* A.2.b     ENDIF
*/
    }

/*
*
* A.3       CALCULATE number of Times 
*           IF (using Verbose Library) DEBUG print
*/ 
   time_cnt=  (Time[Tindex-1] / tau_incr) + 1;
   DPRINT3 ("Time count= (%d / %d) + 1= %d   [",
		Time[Tindex-1], tau_incr, time_cnt);
   for (i=0; i < time_cnt; i++) DPRINT1("%i ",i*tau_incr);
   DPRINT0 ("]\n");
/*
*
* A.4       FOR (each cell in Parameter list) 
*           DO
*/
  /* sum up #levels used by each variables */
  for (recspertime=0, parm_ptr=*parm_head; parm_ptr!=NULL;
		parm_ptr=parm_ptr->next) 
    {
/*
* A.4.1        IF (cell's level is zdef level) THEN
*                 ADD number of heights it has to recspertime
*              ELSE
*                 ADD 1 to recspertime
*              ENDIF
*/
	if (parm_ptr->lvl_ptr->usLevel_id == common->zdef_lvl)
	{
	DPRINT3 ("\t+(%s  Lvlid=%d numhts=%d)\n",parm_ptr->abbrv,
	parm_ptr->lvl_ptr->usLevel_id, parm_ptr->lvl_ptr->numheights);
	recspertime += parm_ptr->lvl_ptr->numheights;
	}
	else 
	{
	DPRINT2("\t+(%s  Lvlid=%d, only 1 ht)\n",parm_ptr->abbrv,
	parm_ptr->lvl_ptr->usLevel_id);
	recspertime += 1;
	}
/*
* A.4       ENDFOR
*/
    }
   DPRINT1("Number of records per time = %d\n", recspertime);

/*
*
* A.5       FUNCTION grad_boundary_box  !calculate coordinates
*           IF (error) THEN
*               PRINT msg
*               RETURN error
*           ENDIF
*/
   if (grad_boundary_box (&Options, &Pdefn, &Xdefn, &Ydefn, gds))
        { 
	  upd_child_errmsg (func, errmsg); GO_HOME;
	}

/*   CREATE GRADS CONTROL FILES  ".gmp' and '.ctl' extensions
 *   Note that Path is chopped off, files are always created
 *   in current directory no matter where input file is;
*/  
/*
*
* A.6       DETERMINE names for the .ctl and .gmp files based
*              on input Grib file's name
*/
   ctl_fn= (char *)malloc(strlen(InFile)+4);
   gmp_fn=(char *)malloc(strlen(InFile)+4);
   if ((cp=strrchr (InFile, '/')) != NULL)
	strcpy (ctl_fn, cp+1);
   else strcpy (ctl_fn, InFile);
   if ((cp=strrchr (ctl_fn, '.')) == NULL) 
	{
        sprintf (gmp_fn, "%s.gmp", ctl_fn); 
        strcat (ctl_fn, ".ctl");
	}
   else {
        ctl_fn[cp-ctl_fn]= '\0';
	strcpy (gmp_fn, ctl_fn);
        strcat(ctl_fn,".ctl");
        strcat(gmp_fn,".gmp");
        }
   fprintf(stdout,"GRADS Control Files= '%s' & '%s'\n", ctl_fn, gmp_fn);
   /*unlink (ctl_fn); unlink (gmp_fn); */

/*
*
* A.7       OPEN files .ctl and 'draw_all.gs' for writing
*           IF (error) THEN
*              PRINT message
*              RETURN error
*           ENDIF
*/
   if ((fctl= fopen (ctl_fn, "w"))==NULL)
        { sprintf(errmsg,"%s:  Open '%s' failed\n",func, ctl_fn); GO_HOME;
        }
   if ((fscript= fopen ("draw_all.gs", "w"))==NULL)  /* Grads Script file */
	{ sprintf(errmsg,"%s: Open 'draw_all.gs' failed\n", func); GO_HOME;
	}

   DPRINT1("\n-- Creating '%s', grad script file 'draw_all.gs'\n", ctl_fn);
   fprintf(fscript,
	"'open %s'\n"
	"my_ctl=\"%s\"   \n"
	"my_grib=\"%s\"   \n"
	"'query file'\n"
	"line1= sublin (result, 1)\n"
	"grib_fn=subwrd (line1, 4)\n"
	"line2= sublin (result, 2)\n"
	"ctl_fn=subwrd (line2, 2)\n"
	"if (ctl_fn != my_ctl | grib_fn != my_grib)\n"
	"   say '**Error** '\n"
	"   say 'File  draw_all.gs  was made for 'my_grib' & 'my_ctl\n"
	"   say 'not for 'grib_fn ' & 'ctl_fn\n"
	"   exit\n"
	"endif\n",
	ctl_fn, ctl_fn, InFile);
	
/*
*
* A.8       PREPARE internal variables
*/
   /* if !null, add <CR> at end of string */
   if (Options[0]!='\0') 
	{ i=strlen(Options); Options[i]='\n'; Options[i+1]='\0'; }
   if (Pdefn[0] != '\0') 	
	{ i=strlen(Pdefn); Pdefn[i]='\n'; Pdefn[i+1]='\0'; }

/*
*
* A.9        WRITE part 1 to .ctl file
*/
   fprintf (fctl , "dset ^%s\ndtype grib\n%s"
        "index ^%s\nundef -9.99E+33\ntitle %s\n"
        "* pdef isz jsz LCC reflat reflon iref jref "
        "stdlat1 stdlat2 stdlon delx dely \n"
        "%s%s\n%s\n"
	,InFile, Options, gmp_fn, InFile, Pdefn, Xdefn, Ydefn);

/*
*
* A.10       FOR (each Level cell in Level List) DO
*               IF (it's the  zdef level) 
*               THEN
*                  IF (it has more than 1 height) 
*                  WRITE all its heights to .ctl file
*                  ELSE write its height and a dummy one to .ctl file
*                  !grads software does not accept levels with 1 height
*               END
*            ENDDO
*/
	/* List all Hieghts for the ZDEF level; */
        for (lvl_ptr=*lvl_head; lvl_ptr!=NULL; lvl_ptr=lvl_ptr->next)
        if (lvl_ptr->usLevel_id==common->zdef_lvl) 
	{
	   if (lvl_ptr->numheights > 1) 
	   {
	     fprintf (fctl,"* Zdef Level is (%03d)\nzdef %d levels\n",
	     common->zdef_lvl, lvl_ptr->numheights);
             for (j=0; j < (int)lvl_ptr->numheights; j++) 
                fprintf (fctl, "%05d ", lvl_ptr->height[j]);
	    }
	   else 
	    {
	     /* GRADS complains if ZDEF is 1, so fake it out with
		a dummy level 
	     */
	     DPRINT0("Adding Dummy Height in Zdef defn \n");
	     fprintf (stdout,
		"* Warning: Zdef level (%03d) only has one Height (%05d)\n"
		, common->zdef_lvl, lvl_ptr->height[0]);
	     fprintf (fctl, 
		"* Level (%03d) only has 1 Height (%05d), add dummy 2nd Height"\
		"\nzdef 2 levels\n%05d %05d", 
		common->zdef_lvl, lvl_ptr->height[0], 
		lvl_ptr->height[0], lvl_ptr->height[0]+1);
            }
           break;
        }

/*
*
* A.11       !Build tdef line for .ctl line
*            IF (there is only 1 forecast time)
* A.11.a     THEN
*               WRITE Tdef line to file with dummy '1mo' Time increment
*               !grads software expects a non-zero Time increment
*/
	/* TDEF format:  taus, hr, day, monthname, year, parmcnt */
        if (time_cnt==1)    /* GRADS expects non-zero increment */
	  {
 	     i= (int)common->base_dtg;
             fprintf(fctl, 
             "\ntdef   %5d linear %02dZ%02d%3s%02d 1mo\n",
	      time_cnt, i%100, (i/100) % 100,
	      month_name[((i /10000)%100) - 1], (i/1000000) % 100);
	  }
/*
* A.11.b     ELSE       !more than one forecast time
*/
        else 
	 {
/*
* A.11.b.1      LET temp be the Time Increment phrase for tdef line
*               ! defaulting to hours
*/
             sprintf (temp, "%dhr", (int)tau_incr);  /* DEFAULT IN HOURS */

	     /* Now see if we can reduce the Tau Increments 
		ie:  increment of 48hrs to  increment of 2days
	     */

	     smallest=99999999; largest=-99999999;
/*
*               ! Try to reduce the Time Incr (ie. 24hr to 1dy)
* A.11.b.2      IF (the smallest time increment is greater or equal to 
*                   #hours in a year)
* A.11.b.2.a    THEN 
*/
             if ((int)tau_incr >= (365*24)) 
		/* increments of Years 
		   requirement is that MMDDHH must always be same;
		*/
		{
/*
* A.11.b.2.a.1      FIND the common MMDDHH 
*/
		  i = (int)(*msgs_head)->base_dtg % 100000; /* common MMDDHH */
/*
* A.11.b.2.a.2      FOR (each mesg in the Mesgs List) DO
*/
		  for (pmesg=(*msgs_head); pmesg; pmesg=pmesg->next)
		  {
/*
* A.11.b.2.a.2.1        FIND the largest and smallest Year
*/
		      /* find largest and smallest YYYY */
		      if ((int)pmesg->base_dtg/1000000 < smallest) 
			  smallest= (int)pmesg->base_dtg/1000000;
		      if ((int)pmesg->base_dtg/1000000 > largest) 
			  largest= (int)pmesg->base_dtg/1000000;
/*
* A.11.b.2.a.2.2        IF (month day and hour differs from the common MMDDHH)
*                       THEN  quit
*                       ENDIF
*/
		      /* quit if  MMDDHH not same */
		      if ((int)pmesg->base_dtg % 1000000 != i) break;
/*
* A.11.b.2.a.2      ENDFOR 
*/
		  } /* for */

/*
* A.11.b.2.a.3      IF (successfully went through entire List) THEN
*                       LET temp be 'Xyr' where X is the difference
*                       in the largest and smallest year
*                   ENDIF
*/
		  if (pmesg == NULL)  /* iff all mmddhh are same */
     	          sprintf (temp, "%dyr", largest-smallest+1);
		}

/*
*
* A.11.b.2.b    ELSE IF (the smallest time increment is greater or equal to 
*                    #hours in a 28-day month) THEN
*/
             else if ((int)tau_incr >= (28*24)) 
		/* increments of months (always 12 months per year)
		   requirement is that DDHH must always be same;
		*/
		{ 
/*
* A.11.b.2.b.1        FIND the common DDHH
*/
		  i = (int)(*msgs_head)->base_dtg % 1000; /* common DDHH */
/*
* A.11.b.2.b.2        FOR (each message in List)
*                     DO
*/
		  for (pmesg=(*msgs_head); pmesg; pmesg=pmesg->next)
		  {
/*
* A.11.b.2.b.2.1         FIND the smallest and largest YYYYMM
*/
		      /* find largest and smallest YYYYMM */
		      if ((int)pmesg->base_dtg/10000 < smallest) 
			  smallest= (int)pmesg->base_dtg/10000;
		      if ((int)pmesg->base_dtg/10000 > largest) 
			  largest= (int)pmesg->base_dtg/10000;
/*
* A.11.b.2.b.2.2         IF (the DDHH is not same as the common DDHH) quit 
*/
		      /* quit if  DDHH not same */
		      if ((int)pmesg->base_dtg % 10000 != i) break;
/*
* A.11.b.2.b.2        ENDDO
*/
		  } /* for */


/*
* A.11.b.2.b.3        IF (went thru entire list successfully) 
*                     THEN
*                         LET temp by 'Xmo' where X is the difference
*                         in months of the largest and smallest YYYYMM
*                     ENDIF
*/
		  /* 	If & Only if  all mmddhh are same then #months 
			between largest yy2mm2 and smallest yy1mm1 dtg is= 
			((yy2-yy1)*12 + mm2) - mm1;
		   */
		  if (pmesg == NULL) sprintf (temp, "%dmo", 
			12 *((largest/100)-(smallest/100)) 
			+ largest%100 - smallest%100); 
		}
/*
*
* A.11.b.2.c    ELSE IF (the smallest time increment is greater or equal to 
*                    #hours in a day) THEN
*/
             else if ((int)tau_incr >= 24) 
	      {  
		 /* increments of days (365/366 days per year)
		    requirement is that HH must always be same;
		 */
/*
* A.11.b.2.c.1        FIND the common HH
*/
		i = (int)(*msgs_head)->base_dtg % 10; /* common 'HH' */
/*
* A.11.b.2.c.2        FOR (each message in List) DO
*/
		for (pmesg=(*msgs_head); pmesg; pmesg=pmesg->next)
		  {
/*
* A.11.b.2.c.2.1         FIND the smallest and largest YYYYMMDD
*/
		      /* find largest and smallest 'yyyymmdd' */
		      if ((int)pmesg->base_dtg/100 < smallest) 
			  smallest= (int)pmesg->base_dtg/100;
		      if ((int)pmesg->base_dtg/100 > largest) 
			  largest= (int)pmesg->base_dtg/100;
/*
* A.11.b.2.c.2.2         IF (the HH is not same as the common HH) quit 
*/
		      /* quit if  'HH' not same */
		      if ((int)pmesg->base_dtg % 100 != i) break;
/*
* A.11.b.2.c.2        ENDDO
*/
		  } /* for */

/*
* A.11.b.2.c.3        IF (went thru entire list successfully) THEN
*/
		/* 	If & Only if  all 'hh' are same then #days 
			between largest yy2mm2dd2 and smallest yy1mm1dd2 is= 
			((yy2-yy1)*12 + mm2) - mm1;
			format for Smalles/Largest= 'yyyymmdd'
		 */
		if (pmesg == NULL) 
		  {
/*
* A.11.b.2.c.3.1         CALCULATE the Year of the smallest and largest dtg
*/
		     yr1= (smallest/10000);
		     yr2= (largest/10000);
/*
* A.11.b.2.c.3.2         IF (they both are the same year)
* A.11.b.2.c.3.2.a       THEN
*                           FIND day_of_year of the smallest dtg
*                           FIND day_of_year of the largest dtg
*                           LET day_cnt be the difference of the 2 day_of_year
*/
		     if (yr1 == yr2)  /* same year */
		      {
		        date.year = yr1; date.month= (smallest/100) % 100;
		        date.day  = smallest % 100; md_doy (&date, &doy1);
		        date.year = yr2; date.month= (largest/100) % 100;
		        date.day  = largest % 100; md_doy (&date, &doy2);
			day_cnt =  doy2 - doy1;
		      }
/*
* A.11.b.2.c.3.2.b       ELSE   
*                           !the 2 dtgs are in different years
*                           !have to count number of days from 1st dtg to
*                           !last dtg, taking Leap yrs into consideration
*/
		     else 	/* dtg1 and dtg2 are in different years */
		      {
			/* Count number of days from first dtg to last dtg,
			   taking Leap year into consideration too
			*/
/*
* A.11.b.2.c.3.2.b.1        CALCULTE the common MM and DD
*/
			date.month= (smallest/100) % 100; /* common month*/
		        date.day  = smallest % 100;  /* common day */
/*
* A.11.b.2.c.3.2.b.2        FOR (each year from 1st to last year)
*                           DO
*                               FIND day_of_year of the common MM, DD using
*                               current year
*                               CALCULATE number of days for this year
*                               IF (curr year is the smallest dtg's year)
*                                    LET day_cnt is the number of days left
*                                    in this year after MM/DD/curryr
*                               ELSE if (curr year is largest dtg's year)
*                                    ADD the day_of_year of MM/DD/curryr
*                                    to day_cnt
*                               ELSE !in between 1st & last yr
*                                    ADD number of days in curryr to day_cnt
*                               ENDIF
*                            ENDDO
*/

			for (yr= yr1; yr <= yr2; yr++)
			{
			  date.year= yr; 
			  md_doy (&date, &doy1);
			  daysperyear= (ISLEAP(yr) ? 366: 355);
			  if (yr == yr1) day_cnt= daysperyear - doy1;
			  else if (yr == yr2) day_cnt += doy1;
			  else day_cnt += daysperyear;
			}
/*
* A.11.b.2.c.3.2.b       ENDIF !different years
*/
		      }
/*
* A.11.b.2.c.3.3         LET temp by 'Xdy' where X is the computed day_cnt
*/
		     sprintf (temp, "%ddy", day_cnt);
/*
* A.11.b.2.c.3        ENDIF  !went thru list successfully
*/
		 }  /* pmesg is null */
/*
* A.11.b.2.c    ENDIF   !check for increments of days
*/
	       } /* increm. of days */

/*
* A.11.b.3      WRITE tdef line to ctl file
*/
 	    i= (int)common->base_dtg;
	     fprintf (fctl, "\ntdef   %5d linear %02dZ%02d%3s%02d %s\n",
	    time_cnt, i%100, (i/100)%100, month_name[((i/10000)%100) - 1],
	    (i / 1000000) % 100, temp); 
/*
* A.11.b     ENDIF    !more than 1 forecasts
*/
	  } /* more than 1 forecast time */

/*
* 
* A.12        WRITE 'vars' line to .ctl file
*/
	fprintf(fctl,
	"* cccccccccccccccccccccccccccccccccccccccccc\n"\
	"* VAR   HTS Pid,Lid[,Ht]    INFO\n"\
	"* cccccccccccccccccccccccccccccccccccccccccc\n"\
	"vars %d\n", parmcnt);
/*
*
* A.13        FOR (each Parm cel in Parameter List) DO
*                 IF (level is non zdef)
*                 THEN 
*                   LIST its abbrv, parmid, lvlid, ht, varnm to .ctl file
*                 ELSE
*                   LIST its abbrv, heightcount, parmid, lvlid, ht, varnm 
*                   to .ctl file
*                 ENDIF
*             DONE
*/


/* 
* 
* A.14        FOR (each parameter in the Parameter list )
*             DO
*/

        for (parm_ptr=*parm_head; parm_ptr!=NULL; parm_ptr=parm_ptr->next)
         {
/*
* A.14.1          IF this is the Level to enumerate (zdef level)
* A.14.1.a        THEN
*/
            if (parm_ptr->usLevel_id == common->zdef_lvl){  /* ZDEF lvl */
/*
* A.14.1.a.1        BUILD variable definition line for .ctl file 
*                   ! Varname #enum_hts  Parmid, Lvlid            Description
*/
              sprintf (temp, "%-8s %2d   %03d,%03d      %s",
              parm_ptr->abbrv, parm_ptr->lvl_ptr->numheights,
              parm_ptr->usParm_id, parm_ptr->usLevel_id, parm_ptr->varnm);
	      fprintf (fctl, "%s\n", temp);

/*
* A.14.1.a.2        WRITE sequence in GRADS script to draw this Variable
*                   for each time, for each Height;
*/
             fprintf (fscript,"say 'ZDEF Level has %d heights,  Z = {",
		parm_ptr->lvl_ptr->numheights);
                for (j=0; j < (int)parm_ptr->lvl_ptr->numheights; j++) 
                fprintf (fscript, " %05d;", parm_ptr->lvl_ptr->height[j]);

	     fprintf (fscript, 
		" }'\ncurr_var=%s\n"
		"heights=%d\n"
		"max_times=%d\n"
		"z=1\n"
		"while (z <= heights)\n"
		"   'set z 'z \n"
		"   'set t 1 'max_times\n"
		"   say '   VARNAME  ENUM PID,LID[,HT]'\n"
		"   say '=> %s'\n" 
		"   say 'draw 'curr_var' at Z='z' for (T=1,'max_times')'  \n"
		"   'display tloop ('curr_var')'\n"
		"   if (rc!=0); break; endif;\n"
		"   z=z+1\n"
		"endwhile\n"
		"say '-------------------------------'\n\n",
	    	parm_ptr->abbrv, parm_ptr->lvl_ptr->numheights, 
		time_cnt, temp);
	      }
/*
* A.14.1.b        THEN
*/
            else {					/* other levels */

/*
* A.14.1.b.1        BUILD variable definition line for .ctl file 
*                   ! Varname 0          Parmid, Lvlid, HEIGHT    Description
*/
              sprintf (temp, "%-8s  0   %03d,%03d,%03d  %s",
              parm_ptr->abbrv, parm_ptr->usParm_id, parm_ptr->usLevel_id,
              parm_ptr->usHeight, parm_ptr->varnm);
	      fprintf (fctl, "%s\n", temp);

/*
* A.14.1.b.2        WRITE sequence in GRADS script to draw this Variable
*/
	      fprintf (fscript, 
		"curr_var=%s\n"
		"heights=%d\n"
		"max_times=%d\n"
		"   'set LEV %d'\n"
		"   'set t 1 'max_times\n"
		"   say '   VARNAME  ENUM PID,LID[,HT]'\n"
		"   say '=> %s'\n" 
		"   say 'drawing 'curr_var' at  LEV=%d (T=1,'max_times')'  \n"
		"   'display tloop ('curr_var')'\n"
		"say '-------------------------------'\n\n",
	    	parm_ptr->abbrv, parm_ptr->lvl_ptr->numheights, 
		time_cnt, parm_ptr->usHeight, 
		temp, parm_ptr->usHeight);
		
/*
* A.14.1.b        ENDIF
*/
	      }

/*
* A.14        ENDFOR !Parameter list
*/
          } /* FOR */

/*
*
* A.15        WRITE last line 'endvars' to .ctl file
*
*/
   fprintf (fctl, "endvars\n");


/*.....  create GRADS .gmp file ......*/
/*
*
* A.16        OPEN .gmp file for writing
*             IF (error) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   fgmp= fopen (gmp_fn, "wb+");
   if (fgmp==NULL) 
	{sprintf(errmsg,"%s:  Unable to open '.gmp' file\n", func); GO_HOME; }

/*
*
* A.17        ALLOCATE storage for GMP_BLK0
*             IF (fail) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   if ((blk0= (GMP_BLK0 *)malloc(sizeof(GMP_BLK0))) == NULL)
        { sprintf(errmsg, "%s: Failed to Malloc BLK1\n", func);GO_HOME; }

   elements= recspertime * time_cnt;
/*
*
* A.18        FILL blk0
*/
   blk0->type=1;      /* for grib */
   blk0->blk1_elements=4;   /* #ints in Hipnt[] */
   blk0->blk2_elements=0;   /* hard coded to 0 */
   blk0->blk3_elements= 3* elements;  /* 3 ints stored perMsg*/
   blk0->blk4_elements= 3* elements;  /* 3 floats stored per msg */
   blk0->notused1= NULL;   /* unknown usage, leave null */
   blk0->notused2= NULL;   /* unknown usage, leave null */
   blk0->notused3= NULL;   /* unknown usage, leave null */
   blk0->notused4= NULL;   /* unknown usage, leave null */
/*
*
* A.19        DEBUG print
*/
   DPRINT2 ("--Creating %s:\nGMP_BLK0 starts at %ld\n",gmp_fn,ftell(fgmp));
/*
* 
* A.20        WRITE blk0 to .gmp file
*             IF (fail) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   if (fwrite (blk0, sizeof(GMP_BLK0), 1, fgmp) != 1) {
	sprintf(errmsg,"%s: Error writing GMP_BLK0 to .gmp file\n", func);
	GO_HOME;
	}
/*
*  
* A.21        FREE storage for blk0
*/
   free(blk0);

/*
* 
* A.22        ALLOCATE storage for blk1 of .gmp file
*             IF (fail) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   if ((blk1= (GMP_BLK1 *)malloc(sizeof(GMP_BLK1))) == NULL)
	 { sprintf(errmsg,"%s: Failed to Malloc BLK1\n", func);exit(0);}
/*
*
* A.23        FILL blk1 block
*/
   blk1->filetype = 1;                /* means grib */
   blk1->tdef = time_cnt;         /* Tdef */
   blk1->recspertime = recspertime;      /* each var at all levels */
   blk1->usGrid_id = pds->usGrid_id;   /* gds ident */

/*
*
* A.24        DEBUG print
*/
   DPRINT1 ("GMP_BLK1 starts at %ld\n", ftell(fgmp));
/*
* 
* A.25        WRITE blk0 to .gmp file
*             IF (fail) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   if (fwrite (blk1, sizeof(GMP_BLK1), 1, fgmp) != 1) {
	sprintf(errmsg,"%s: Error writing GMP_BLK1 to .gmp file\n", func);
	GO_HOME; }
/*
*  
* A.26        FREE storage for blk1
*/
   free(blk1);

   DPRINT9("*** Recspertau= %d * #times= %d) -> %d Elements\n"
        "*** blk3 has [%d elements] * (%d bytes/elem)= total %d bytes\n"
        "*** blk4 has [%d elements] * (%d bytes/elem)= total %d bytes\n",
	recspertime,time_cnt, elements,
        elements, sizeof(GMP_BLK3), elements,
        elements, sizeof(GMP_BLK4),  elements);

/*
*
*             ! blk 2 is left out on purpose
*
* A.27        ALLOCATE space for blk3
*             IF (fail) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   blk3= (GMP_BLK3 *)malloc (elements* sizeof(GMP_BLK3));
   if (blk3 ==NULL )
   { sprintf(errmsg,"%s: Failed to Malloc BLK3\n",func); GO_HOME; }
/*
*
* A.28        ALLOCATE space for blk3
*             IF (fail) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   blk4= (GMP_BLK4 *)malloc (elements* sizeof(GMP_BLK4));
   if (blk4 ==NULL)
   { sprintf(errmsg,"%s: Failed to Malloc BLK4\n",func); GO_HOME; }

/*
*
* A.29        INIT all elements of blk3 and blk4 to missing
*/
   for (i=0; i < elements; i++) /* init all to missing */
       { 
	(blk3+i)->dpos= -999; 
	(blk3+i)->bpos= -999; 
	(blk3+i)->bnum= -999; 
  	(blk4+i)->fDec_sc_fctr=0.0; 
	(blk4+i)->fBin_sc_fctr=0.0; 
  	(blk4+i)->fReference=0.0; 
  	}

/*
*
* A.30        FOR (each grad message) DO
*                 IF (base_dtg is zero OR  parm_ptr is null OR abbreviated
*                     name is undefined) 
*                 THEN
*                     RETURN error code   !null record
*                 ENDIF
*             ENDDO
*             DEBUG print
*/
  for (i=1, pmesg=(*msgs_head); pmesg!=NULL; pmesg=pmesg->next, ++i)
      if (!pmesg->base_dtg || !pmesg->parm_ptr || !pmesg->parm_ptr->abbrv[0]) 
	{ 
	  sprintf(errmsg, "%s:  NULL RECORD Record#%d\n", func,i); GO_HOME;
	}
   DPRINT0 ("Building Blk3 and Blk4 now;\n");
/*
*
*             !Building Blk3 and Blk4, Info of each Field present
*             !for (var=firstvar (T=1 (Ht=1, maxhts), numtimes), lastvar)
*             !loop order:  innermost  to outermost
*
* A.31        FOR (each mesg cell in Msg list) 
*             DO
*/
   for (pmesg= (*msgs_head), i=1; pmesg != NULL; pmesg=pmesg->next) 
   {
     DPRINT1 ("\n-->  LOOKING FOR=  %s\n", pmesg->parm_ptr->abbrv);

/*
* A.31.1          IF (curr mesg's basetime is zero OR parm_ptr is
*                     NULL or abbrv is NULL)
*                 THEN
*                      PRINT message
*                      RETURN 1   !error
*                 ELSE INCREMENT counter
*/
     if (pmesg->base_dtg==0|| pmesg->parm_ptr==NULL||
	pmesg->parm_ptr->abbrv[0]=='\0') {
	sprintf(errmsg, "%s: ERROR: RECORD #%d (from 1) is null\n",
	func, pmesg-*msgs_head);
	GO_HOME;
	}
	else ++i;

/*
* A.31.2          CALCULATE the grad Time for this message
*/
     time= (int)(pmesg->tau_incr / tau_incr);  /* Time Range= [0-time_cnt] */
     DPRINT1("Current Time=%d\n", time);

/*
*                 ! 'VarHt_Offs' is the accrual of #recs (each is a Height)
*                 !  from ALL previous Parameters ;
*                 ! It is used as Index where to the information in the
*                 ! .gmp file (within in Int/Flt block)
*
* A.31.3          FOR (each parameter in Parameter list) 
*                 DO
*/
     for (varht_offs=0, parm_ptr=*parm_head; 
			parm_ptr!=NULL; parm_ptr=parm_ptr->next)
      { 
/*
* A.31.3.1            IF (cell's abbrv matches mesg's abbrv) 
* A.31.3.1.a          THEN 
*                        BREAK   !found the correct Parm cell in list
* A.31.3.1.b          ELSE 
*                        !Calculate offset of this var within Parameter List
*                        IF (this level is the ZDEF level) 
*                        THEN
*                            ADD number of heights of this parm to varht_offs
*                        ELSE
*                            JUST add 1 to varht_offs
*                        ENDIF
*                     ENDIF
*/
	if (!strcmp (pmesg->parm_ptr->abbrv, parm_ptr->abbrv)) break;

	else if (parm_ptr->lvl_ptr->usLevel_id == common->zdef_lvl)
	     {
	     	varht_offs += parm_ptr->lvl_ptr->numheights;
	        DPRINT3 (
	         "bump varht_offs past %d hts of var '%s', new varhtofs=%d\n",
	         parm_ptr->lvl_ptr->numheights, parm_ptr->abbrv, varht_offs);
	     }
        else {
		varht_offs += 1;
		DPRINT2 ("bump varht_offs past var '%s', new varhtofs=%d\n",
	        parm_ptr->abbrv, varht_offs);
	     }
		
/*
* A.31.3          ENDDO
*/
      }

/*
* A.31.4          FOR (each Height of this level)
*                 DO
*                     IF (it matches the height of curr mesg)
*                     THEN break   !found correct height
*                     ELSE
*                        BUMP varht_offs by 1  !keep looping
*                     ENDIF
*                 ENDDO
*/
     for (i=0; i < (int)parm_ptr->lvl_ptr->numheights; i++)
	if (pmesg->usHeight1 == parm_ptr->lvl_ptr->height[i]) {
	   DPRINT3 ("Found '%s' ht=%d, use varht_ht=%d\n",
	   parm_ptr->abbrv, parm_ptr->lvl_ptr->height[i], varht_offs);
	   break;
	   }
	else if (pmesg->parm_ptr->usLevel_id == common->zdef_lvl)
	   {
	   varht_offs += 1;
	   DPRINT2 ("bump varht_offs past Height %d, varhtoffs=%d\n",
	   parm_ptr->lvl_ptr->height[i], varht_offs); 
	   }

/*
* A.31.5          CALCULATE offset where message should go
*                 !offset = (time * recspertime + varht_offs);
*                 IF (out of limit) THEN
*                     PRINT message
*                     RETURN 1 !error
*                 ENDIF
*/
     offset= (time * recspertime + varht_offs);
     if (offset<0 || offset>= elements) {
	sprintf(errmsg,
		"%s WARNING:  BAD OFFSET=%d > max_off=%d\n"
		"Time=%d, recspertime=%d, varht=%d\n"
	   	"MESG's %s dtg=%ld tau=%d parm=%d lvl=%d ht=%d TAUINCR=%d\n",
		func,
		offset,elements, time, recspertime, varht_offs,
	  	pmesg->parm_ptr->abbrv, 
	  	pmesg->base_dtg, pmesg->ustau, pmesg->parm_ptr->usParm_id,
	  	pmesg->parm_ptr->usLevel_id, pmesg->usHeight1,
	  	pmesg->tau_incr);
	GO_HOME;
	}

/*
* A.31.6          FILL blk3 and blk4
*/
     blk3[offset].dpos= (int)pmesg->dpos;
     blk3[offset].bpos= (int)pmesg->bpos;
     blk3[offset].bnum= (int)pmesg->bnum;

     blk4[offset].fDec_sc_fctr= pmesg->fDec_sc_fctr;
     blk4[offset].fBin_sc_fctr= pmesg->fBin_sc_fctr;
     blk4[offset].fReference  = pmesg->fReference;
/*
	DPRINT3 ("OFFS=%d, REFENRECE= %f, blk4[offs].fRef= %f\n", 
	offset, pmesg->fReference, blk4[offset].fReference);
*/

/*
* A.31.7          DEBUG print
*/
     DPRINT7("Make: T=%d: Varht=%d <OFFS=%d> %010d tau=%d '%s' Ht=%d\n",
	time+1, varht_offs, offset, pmesg->base_dtg, pmesg->ustau,
  	pmesg->parm_ptr->abbrv, pmesg->usHeight1);

	DPRINT7 ("OFS=%d (BDS=%d,BMS=%d,#bits=%d : D=%f,Scl=%f,Ref=%f)\n",
	offset,
	(blk3+offset)->dpos,(blk3+offset)->bpos,(blk3+offset)->bnum,
	(blk4+offset)->fDec_sc_fctr,(blk4+offset)->fBin_sc_fctr,
	(blk4+offset)->fReference);
	
/*
* A.31        ENDFOR 
*/
   }

/*
*
* A.32        DEBUG print
*/
   DPRINT0("\n***** FINAL LIST:\n");
   for (offset=0; offset < elements; offset++) 
     if (blk3[offset].dpos== -999 && blk3[offset].bpos== -999 &&
	    blk3[offset].bnum== -999 && blk4[offset].fDec_sc_fctr==0. &&
	    blk4[offset].fBin_sc_fctr==0. && blk4[offset].fReference==0.)
	{
        DPRINT2("OFS=%02d: T=%d > missing;\n", 
	offset, (offset / recspertime) + 1); }
     else  {
	DPRINT8("OFS=%02d: T=%d > (%d,%d,%d : %f,%f,%f)\n",
	offset,
  	(offset / recspertime) + 1, 
	blk3[offset].dpos,blk3[offset].bpos,blk3[offset].bnum,
	blk4[offset].fDec_sc_fctr,blk4[offset].fBin_sc_fctr,
	blk4[offset].fReference);
	}

/*
*
* A.33        WRITE blk3 to .gmp file
*             IF (fail) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   DPRINT1 ("GMP_BLK2 has no bytes\nGMP_BLK3 starts at %ld\n", ftell(fgmp));
   if (fwrite (blk3, sizeof(GMP_BLK3), elements, fgmp) != elements) {
	sprintf(errmsg,"%s: Error writing GMP_BLK3 to .gmp file\n", func);
	GO_HOME;
	}

/*
*
* A.34        WRITE blk4 to .gmp file
*             IF (fail) THEN
*                 PRINT message
*                 RETURN error
*             ENDIF
*/
   DPRINT1 ("GMP_BLK4 starts at %ld\n", ftell(fgmp));
   if (fwrite (blk4, sizeof(GMP_BLK4), elements, fgmp) != elements) {
        sprintf(errmsg,"%s: Error writing GMP_BLK4 to .gmp file\n", func);
	GO_HOME;
	}
   DPRINT2 ("--end of %s: %ld bytes total--\n",  gmp_fn,ftell(fgmp));

/*
*
* A.35        SET no error status
*/
   Ecode = 0;


BAIL_OUT:
  fprintf(stdout,"at BAIL_OUT:  errmsg=%s\n", errmsg);
/*
*
* A.36        CLOSE files if still opened
*/
   if (fctl!= NULL) fclose (fctl);  
   if (fgmp!= NULL) fclose (fgmp); 
   if (fscript!= NULL) fclose (fscript); 

   if (ctl_fn!=NULL) free(ctl_fn); 
   if (gmp_fn!=NULL) free(gmp_fn); 
   if (blk0 != NULL) free(blk0);  
   if (blk1 != NULL) free(blk1); 
   if (blk3 != NULL) free(blk3);  
   if (blk4 != NULL) free(blk4);
   if (Pdefn != NULL) free(Pdefn);
   if (Xdefn != NULL) free(Xdefn);
   if (Ydefn != NULL) free(Ydefn);
   if (Options != NULL) free(Options);

/*
*
* A.37        SET no error status & RETURN
*/
   DPRINT2 ("Leaving %s, Stat %d\n", func, Ecode); 
   return(Ecode); 
/*
*
* END OF FUNCTION
*/
}
