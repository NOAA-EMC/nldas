#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "preproc.h"
#include "vicNl_def.h"

/***********************************************************************
Greg O'Donnell [tempgd@hydro.washington.edu], June, 2000
 
Make filenames to open the grib files  

***********************************************************************/

char * make_fname(char *path, char *prefix, char *suffix, dmy_struct *dmy )
{

  /* strcat and strcpy require space for NULL */
  char cyear[5];
  char cmonth[3];
  char cday[3];
  char chour[3];

  static char filename[BUFSIZ+1];

  /* write the dmy timestamp to characters */
  sprintf(cyear, "%4.4d", dmy->year);
  sprintf(cmonth,"%2.2d", dmy->month);
  sprintf(cday,  "%2.2d", dmy->day);
  sprintf(chour, "%2.2d", dmy->hour);

  /* piece together the filename */
  strcpy(filename,path);
  strcat(filename,prefix);
  strcat(filename,cyear);
  strcat(filename,cmonth);
  strcat(filename,cday);
  strcat(filename,chour);
  strcat(filename,suffix);

  return filename;
}
