#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: check_files.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void check_files(infiles_struct   *infp, 
		 filenames_struct *fnames)
/**********************************************************************
	check_files		Dag Lohmann		January 1996

  This routine opens files for soil, vegitation, and global parameters.

**********************************************************************/
{
  extern option_struct  options;
  extern FILE          *open_file(char string[], char type[]);

  infp->soilparam   = open_file(fnames->soil, "r");
  infp->veglib      = open_file(fnames->veglib, "r");
  infp->vegparam    = open_file(fnames->veg, "r");
  if(options.SNOW_BAND>1)
    infp->snowband    = open_file(fnames->snow_band, "r");
}



