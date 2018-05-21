#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "preproc.h"
 
/**********************************************************************
  get_params	Greg O'Donnell	            May 2000

  This routine reads the global control file
  [based on vic].
  
**********************************************************************/
void global_params ( global_param_struct *global, FILE *gp )
{
  
  char cmdstr[BUFSIZ+1];
  char optstr[BUFSIZ+1];
  char flgstr[BUFSIZ+1];

  /** Initialize global parameters **/

  /* defaults */
  global->ncep_out = 1;   /* output to binary by default */

  global->startyear	= MISSING;
  global->startmonth	= MISSING;
  global->startday	= MISSING;
  global->starthour	= MISSING;

  global->endyear	= MISSING;
  global->endmonth	= MISSING;
  global->endday	= MISSING;
  global->endhour	= MISSING;

  global->nrecs		= MISSING;

  /* initialise character strings */
  global->ncep_met[0] = (char)NULL;
  global->forcing[0]  = (char)NULL;
  global->mask[0]     = (char)NULL;
  global->ncep_ext[0] = (char)NULL;


  /** Read through global control file to find parameters **/
  fgets(cmdstr,BUFSIZ,gp);

  while(!feof(gp)) {
    /* skip blank lines and those starting with '#' */
    if(cmdstr[0]!='#' && cmdstr[0]!='\n' && cmdstr[0]!='\0') {

      sscanf(cmdstr,"%s",optstr);

      /*******************************
        Get Model Global Parameters
	*****************************/
      if(strcasecmp("TIME_STEP",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->dt);
      }
      else if(strcasecmp("NRECS",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->nrecs);
      }
      else if(strcasecmp("STARTYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->startyear);
      }
      else if(strcasecmp("STARTMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->startmonth);
      }
      else if(strcasecmp("STARTDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->startday);
      }
      else if(strcasecmp("STARTHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->starthour);
      }
      else if(strcasecmp("ENDYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->endyear);
      }
      else if(strcasecmp("ENDMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->endmonth);
      }
      else if(strcasecmp("ENDDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->endday);
      }
      else if(strcasecmp("ENDHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global->endhour);
      }
      else if(strcasecmp("NCEP_MET",optstr)==0) {
        sscanf(cmdstr,"%*s %s",global->ncep_met);
      }
      else if(strcasecmp("FORCING1",optstr)==0) {
        sscanf(cmdstr,"%*s %s",global->forcing);
      }
      else if(strcasecmp("MASK",optstr)==0) {
        sscanf(cmdstr,"%*s %s",global->mask);
      }
      else if(strcasecmp("NCEP_QCLOG",optstr)==0) {
        sscanf(cmdstr,"%*s %s",global->qclog);
      }
      else if(strcasecmp("NCEP_EXT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",global->ncep_ext);
      }
      else if(strcasecmp("FORCE_FORMAT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
	if(strcmp(flgstr,"BINARY")==0)
	  global->ncep_out=1;
	else if(strcmp(flgstr,"ASCII")==0)
	  global->ncep_out=0;
	else{
	  fprintf(stderr,"Unrecognised option for \"FORCE_FORMAT\"\t%s\n", flgstr);
	  exit(EXIT_FAILURE);
	}
      }
      else {
	/*	fprintf(stderr,"WARNING: Unrecognized option in the global
		parameter file:\n\t%s is unknown - check your spelling\n",
		optstr); */
      }
    }
    fgets(cmdstr,BUFSIZ,gp);
  }

  /******************************************
    Check for undefined required parameters
  ******************************************/

  if(global->dt == MISSING ) {
     fprintf(stderr,"Must define time step <DT> in control file.\n");
     exit(EXIT_FAILURE);
  }
  if( global->startyear == MISSING || global->startmonth == MISSING ||
      global->startday  == MISSING || global->startyear  == MISSING ){
     fprintf(stderr,"Must define startdate as <year> <month> <day> <hour>\n");
     exit(EXIT_FAILURE);
  }

  /* specify either the nrecs or the end date */
  if( global->nrecs == MISSING ){
      if( global->endyear == MISSING || global->endmonth == MISSING ||
          global->endday  == MISSING || global->endyear  == MISSING ){
        fprintf(stderr,"Must specify EITHER <NRECS> or an end date\n");
      exit(EXIT_FAILURE);
      }
   }
   else if ( global->nrecs != MISSING ){
      if( global->endyear != MISSING || global->endmonth != MISSING ||
          global->endday  != MISSING || global->endyear  != MISSING ){
        fprintf(stderr,"Must specify EITHER <NRECS> or an end date\n");
      exit(EXIT_FAILURE);
      }
   }

   /* check the file name paths */
   if ( global->ncep_met[0] == (char)NULL ){
     fprintf(stderr,"NCEP met path <NCEP_MET> must be specified\n");
     exit(EXIT_FAILURE);
   }
   if ( global->forcing[0] == (char)NULL ){
     fprintf(stderr,"Processed met file path <FORCING> must be specified\n");
     exit(EXIT_FAILURE);
   }
   if ( global->mask[0] == (char)NULL ){
     fprintf(stderr,"Computational <MASK> file must be specified\n");
     exit(EXIT_FAILURE);
   }
   if ( global->ncep_ext[0] == (char)NULL ){
     fprintf(stderr,"NCEP filename extension must be specified\n");
     exit(EXIT_FAILURE);
   }

}
