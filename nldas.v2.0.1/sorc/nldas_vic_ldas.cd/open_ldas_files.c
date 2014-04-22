#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"
#include "vicNl_ldas.h"

/***********************************************************************
  Author: Greg O'Donnell [tempgd@hydro.washington.edu], June 2000

  open_ldas_files.c

  This routine opens the ldas met files and the output files.
************************************************************************/

/**** call routines to open data files ****/
void open_ldas_files( ldas_metfiles_struct *ldas_metfiles, 
		      ldas_outfiles_struct *ldas_outfiles,
		      filenames_struct     *fname,
		      dmy_struct           *dmy )
{
  /* open the met files */
  open_ldas_metfiles(ldas_metfiles, fname, dmy, "r");
  open_ldas_outfiles(ldas_outfiles, fname, dmy, "w");
}

