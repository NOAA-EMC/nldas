#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: open_state_file.c,v 4.2 2000/05/16 21:57:54 vicadmin Exp vicadmin $";

#if SAVE_STATE

FILE *open_state_file(global_param_struct *global,
		      int                  Nlayer,
		      int                  Nnodes) 
/*********************************************************************
  open_state_file      Keith Cherkauer           April 15, 2000

  This subroutine opens the model state file if used.

  Modifications:
  07/20/00 Removed the date from the state file name to ease its 
           manipulation when VIC is run in real-time at NCEP      JS
*********************************************************************/
{
  extern option_struct options;

  FILE   *statefile;
  char    filename[MAXSTRING];
  double  Nsum;

  /* open state file */
/** Start of changes - JS - removed date from state file name */
  sprintf(filename,"%s", global->statename);
/*  sprintf(filename,"%s_%02i%02i%04i", global->statename, 
	  global->stateday, global->statemonth, global->stateyear);*/
/* End of changes - JS */
  statefile = open_file(filename,"w");

  /* Write save state date information */

  /** Start of changes - NGMOD **/
  /* write out '#' in header info - auto skip in open_file */
    fprintf(statefile,"# %i %i %i\n", global->stateday, 
	    global->statemonth, global->stateyear);
    fprintf(statefile,"# %i %i\n", Nlayer, Nnodes);
  /** End of changes - NGMOD **/
  /* Write simulation flags */

  return(statefile);

}

#endif
