#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "vicNl.h"
#include "vicNl_ldas.h"

/*
  Author: Greg O'Donnell [tempgd@hydro.washington.edu], June 2000

   Re-read the global file and get the LDAS specific options.
   The code is set up this was to avoid changing get_global_param.c

*/

static char Routine[] = "get_ldas_global_param";

void get_ldas_global_param(FILE *gp, ldas_option_struct *ldas_options )
{

  char cmdstr[BUFSIZ+1];
  char optstr[BUFSIZ+1];
  char flgstr[BUFSIZ+1];

  /* reset after the get_global_param read statements */
  rewind(gp);

  /* make sure the stream is still open and readable */
  fgets(cmdstr,BUFSIZ,gp);
  if(errno){
    fprintf(stderr,"Error in\t'%s'\n",Routine);
    perror(strerror(errno));
    exit(EXIT_FAILURE);
  }

  /* initialise variables */
  ldas_options->mask[0]           = '\0';
  ldas_options->processor_mask[0] = '\0';
  ldas_options->pweight_mask[0]   = '\0';
  ldas_options->pp_output[0]      = '\0';
  ldas_options->lats_table[0]     = '\0';
  ldas_options->processor         = 0;
  ldas_options->endhour           = -99;
  ldas_options->statehour         = 0;
 
  ldas_options->ncep_met[0] = (char)NULL;
  ldas_options->forcing[0]  = (char)NULL;
  ldas_options->mask[0]     = (char)NULL;
  ldas_options->ncep_ext[0] = (char)NULL;

  /* process necessary options */
  while(!feof(gp)) {
    if(cmdstr[0]!='#' && cmdstr[0]!='\n' && cmdstr[0]!='\0') {   

      sscanf(cmdstr,"%s",optstr); 
      if(strcasecmp("MASK",optstr)==0) {
        sscanf(cmdstr,"%*s %s", ldas_options->mask);
      }   
      else if(strcasecmp("PROCESSOR_MASK",optstr)==0) {
        sscanf(cmdstr,"%*s %s", ldas_options->processor_mask);
      }   
      else if(strcasecmp("PROCESSOR",optstr)==0) {
        sscanf(cmdstr,"%*s %i", &ldas_options->processor);
      }
      else if(strcasecmp("PWEIGHT_MASK",optstr)==0) {
        sscanf(cmdstr,"%*s %s", ldas_options->pweight_mask);
      }
      else if(strcasecmp("PP_OUTPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s", ldas_options->pp_output);
      }
      else if(strcasecmp("LATS_TABLE",optstr)==0) {
        sscanf(cmdstr,"%*s %s", ldas_options->lats_table);
      }
      /* there isnt a better place for this */
      else if(strcasecmp("ENDHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&ldas_options->endhour);
      }
      /* there isnt a better place for this */
      else if(strcasecmp("STATEHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&ldas_options->statehour);
      }
      else if(strcasecmp("NCEP_MET",optstr)==0) {
        sscanf(cmdstr,"%*s %s",&ldas_options->ncep_met[0]);
      }
      else if(strcasecmp("FORCING1",optstr)==0) {
        sscanf(cmdstr,"%*s %s",&ldas_options->forcing[0]);
      }
      else if(strcasecmp("MASK",optstr)==0) {
        sscanf(cmdstr,"%*s %s",&ldas_options->mask[0]);
      }
      else if(strcasecmp("NCEP_QCLOG",optstr)==0) {
        sscanf(cmdstr,"%*s %s",&ldas_options->qclog[0]);
      }
      else if(strcasecmp("NCEP_EXT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",&ldas_options->ncep_ext[0]);
      }
      else if(strcasecmp("FORCE_FORMAT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
	if(strcmp(flgstr,"BINARY")==0)
	  ldas_options->ncep_out=1;
	else if(strcmp(flgstr,"ASCII")==0)
	  ldas_options->ncep_out=0;
	else{
	  fprintf(stderr,"Unrecognised option for \"FORCE_FORMAT\"\t%s\n", flgstr);
	  exit(EXIT_FAILURE);
	}
      }
      else{
	/* holder */
	;
      }
    }
    fgets(cmdstr,BUFSIZ,gp);   
  }

  /* check the file name paths */
  if ( ldas_options->ncep_met[0] == (char)NULL ){
    fprintf(stderr,"NCEP met path <NCEP_MET> must be specified\n");
    exit(EXIT_FAILURE);
  }
  if ( ldas_options->forcing[0] == (char)NULL ){
    fprintf(stderr,"Processed met file path <FORCING> must be specified\n");
    exit(EXIT_FAILURE);
  }
  if ( ldas_options->mask[0] == (char)NULL ){
    fprintf(stderr,"Computational <MASK> file must be specified\n");
    exit(EXIT_FAILURE);
  }
  if ( ldas_options->ncep_ext[0] == (char)NULL ){
    fprintf(stderr,"NCEP filename extension must be specified\n");
    exit(EXIT_FAILURE);
  }

  /* check all all options are present */
  if(ldas_options->processor == 0){
    fprintf(stderr,"Must specify 'PROCESSOR' in global file.\n");
    exit(EXIT_FAILURE);
  }
  else if ( ldas_options->mask[0]=='\0'){
    fprintf(stderr,"Must specify 'MASK' in global file.\n");
    exit(EXIT_FAILURE);
  }
  else if ( ldas_options->lats_table[0]=='\0'){
    fprintf(stderr,"Must specify 'LATS_TABLE' in global file.\n");
    exit(EXIT_FAILURE);
  }
  else if ( ldas_options->processor_mask[0]=='\0') {
    fprintf(stderr,"Must specify 'PROCESSOR_MASK' in global file.\n");
    exit(EXIT_FAILURE);
  }
  else if ( ldas_options->pweight_mask[0]=='\0') {
    fprintf(stderr,"Must specify 'PWEIGHT_MASK' in global file.\n");
    exit(EXIT_FAILURE);
  }
  else if ( ldas_options->pp_output[0]=='\0') {
    fprintf(stderr,"Must specify 'PP_OUTPUT' in global file.\n");
    exit(EXIT_FAILURE);
  }
  else if(ldas_options->endhour < 0){
    fprintf(stderr,"Must specify 'ENDHOUR' in global file.\n");
    exit(EXIT_FAILURE);
  }  
}
