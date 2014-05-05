#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include "vicNl_ldas.h"

/*
  Author: Greg O'Donnell [tempgd@hydro.washington.edu], June 2000

  Process the input masks [set by MASK and PROCESSOR_MASK in
  the global file].

  the 'mask' file contains all cell numbers in the met forcing files.
  The 'processor_mask' file contains those cells to be done, using 
  the 'processor' flag.
  
  Modifications:
  	092000	Store the cell locations (USA=1, Canada=2 or Mexico=3)
  		in the ldas_index structure
  													JS
  	030501	Read in precip weighting map and store these in 
  		weighting array ot be used later. This replaces cell
  		location (USA, Canada, Mexico) array.
  													JS

*/

/* Prototypes for local funcunctions */
void alloc_2Darray( int ***array, int ncols, int nrows);
void free_2Darray ( int ***array, int nrows);


void process_masks(ldas_option_struct *ldas_options, ldas_index_struct *ldas_index)
{

  FILE *mfp;
  FILE *pmfp;
/** Start of changes -- JS **/
  FILE *pwfp;
/** End of changes -- JS **/
  
  int  pmask_nr;   /* processor mask file rows */
  int  pmask_nc;   /* processor mask file rows */
/** Start of changes -- JS **/
  int  pweight_nr; /* processor mask file rows */
  int  pweight_nc; /* processor mask file rows */
/** End of changes -- JS **/
  int  missing;    /* mask file missing value */

  int **mask_grid;  /* grid containing mask */
  int **pmask_grid; /* grid containing processor mask */
/** Start of changes -- JS **/
  int **pweight_grid; /* grid containing precip weights */
/** End of changes -- JS **/

  int i;
  int j;
  int k;
  int m;

  /* open mask files */
  mfp  = open_file(ldas_options->mask,          "r");
  pmfp = open_file(ldas_options->processor_mask,"r");
/** Start of changes -- JS **/
  pwfp = open_file(ldas_options->pweight_mask,"r");
/** End of changes -- JS **/

  /* read and compare headers */
  /*read mask header */
  fscanf(mfp,"%*s %d", &ldas_index->nc);
  fscanf(mfp,"%*s %d", &ldas_index->nr);
  fscanf(mfp,"%*s %f", &ldas_index->xll);
  fscanf(mfp,"%*s %f", &ldas_index->yll);
  fscanf(mfp,"%*s %f", &ldas_index->size);
  fscanf(mfp,"%*s %d", &missing); 

  /*read processor mask header */
  fscanf(pmfp,"%*s %d", &pmask_nc);
  fscanf(pmfp,"%*s %d", &pmask_nr);
  fscanf(pmfp,"%*s %*f");
  fscanf(pmfp,"%*s %*f");
  fscanf(pmfp,"%*s %*f");
  fscanf(pmfp,"%*s %*d"); 

/** Start of changes -- JS **/
  /*read precip weights header */
  fscanf(pwfp,"%*s %d", &pweight_nc);
  fscanf(pwfp,"%*s %d", &pweight_nr);
  fscanf(pwfp,"%*s %*f");
  fscanf(pwfp,"%*s %*f");
  fscanf(pwfp,"%*s %*f");
  fscanf(pwfp,"%*s %*d"); 
/** End of changes -- JS **/

  if(ldas_index->nc != pmask_nc || ldas_index->nr != pmask_nr ){
    fprintf(stderr,"Incompatible mask dimensions in %s and %s\n", 
	    ldas_options->mask, ldas_options->processor_mask);
    exit(EXIT_FAILURE);
  }

/** Start of changes -- JS **/
  if(ldas_index->nc != pweight_nc || ldas_index->nr != pweight_nr ){
    fprintf(stderr,"Incompatible mask dimensions in %s (%d x %d) and %s (%d x %d)\n", 
	    ldas_options->mask, ldas_index->nc, ldas_index->nr, ldas_options->pweight_mask, pweight_nc, pweight_nr);
    exit(EXIT_FAILURE);
  }  
/** End of changes -- JS **/
  
  /* read the masks into memory */
  alloc_2Darray(&mask_grid,  ldas_index->nc, ldas_index->nr); 
  alloc_2Darray(&pmask_grid, pmask_nc, pmask_nr); 
  alloc_2Darray(&pweight_grid, pweight_nc, pweight_nr); 

  /* read grids into memory - upside down - as NCEP grib data is stored*/
  /* mask grid */
  ldas_index->ncell_met=0;
  for(i=ldas_index->nr-1;i>=0;i--){
    for(j=0;j<ldas_index->nc;j++){
      fscanf(mfp,"%d",&mask_grid[i][j]);
      if(mask_grid[i][j]!=missing){
/** Start of changes -- JS **/
          /* remove this line so that mask_grid retains the actual value from the file */
/*        mask_grid[i][j]=1;*/
/** End of changes -- JS **/
        ldas_index->ncell_met++;
      }
      else
        mask_grid[i][j]=0;
    }
  }       
  /* processor mask grid */
  ldas_index->ncell_comp=0;
  for(i=pmask_nr-1;i>=0;i--){
    for(j=0;j<pmask_nc;j++){
      fscanf(pmfp,"%d",&pmask_grid[i][j]);
      if(pmask_grid[i][j]==ldas_options->processor){
        pmask_grid[i][j]=1;
	if(mask_grid[i][j]==0){
	  fprintf(stderr,
		  "Processor mask [%d][%d] is computational and mask is not\n",
		  i, j);
	  exit(EXIT_FAILURE);
	}
        ldas_index->ncell_comp++;
      }
      else
        pmask_grid[i][j]=0;
    }
  }
  /* precip weights grid */
  for(i=pweight_nr-1;i>=0;i--){
    for(j=0;j<pweight_nc;j++){
      fscanf(pwfp,"%d",&pweight_grid[i][j]);
    }
  }

#if 1
  fprintf(stderr,"# of cells in mask, and # of cells to run\t%d %d\n",
	  ldas_index->ncell_met, ldas_index->ncell_comp);
#endif

  /* allocate memory for list of computational cells */
  /* ldas_index.cell_num are the cell numbers in the soil file to run */
  if( !((*ldas_index).cell_num = (int *)calloc(ldas_index->ncell_comp,sizeof(int))) ){
    fprintf(stderr,"Cannot allocate memory to ldas_index.cell_num\n");
    exit(EXIT_FAILURE);
  }
  /* ldas_index.cell_index is used to extract the met data */
  if( !((*ldas_index).cell_index = (int *)calloc(ldas_index->ncell_comp,sizeof(int))) ){
    fprintf(stderr,"Cannot allocate memory to ldas_index.cell_index\n");
    exit(EXIT_FAILURE);
  }
/** Start of changes -- JS **/
  /* ldas_index.cell_pweight is for precip weights */
  if( !((*ldas_index).cell_pweight = (int *)calloc(ldas_index->ncell_comp,sizeof(int))) ){
    fprintf(stderr,"Cannot allocate memory to ldas_index.cell_pweight\n");
    exit(EXIT_FAILURE);
  }
/** End of changes -- JS **/

  /* create a list of cells to be computed */
  k = m = 0;
  for(i=0;i<pmask_nr;i++){
    for(j=0;j<pmask_nc;j++){
      if(pmask_grid[i][j]){
        (*ldas_index).cell_num[k]=1+i*pmask_nc+j;
        (*ldas_index).cell_index[k]=m;
/** Start of changes -- JS **/
        /* store the precip weights */
        (*ldas_index).cell_pweight[k]=pweight_grid[i][j];
/** End of changes -- JS **/
        k++;
      }
/** Start of changes -- JS **/
/*      if(mask_grid[i][j]==1)m++;*/
      /* mask may have location information in it so test for all locations >= 1 */
      if(mask_grid[i][j]>=1)m++;
/** End of changes -- JS **/
    }
  }

  
  /* free memory allocated to grids */
  free_2Darray(&mask_grid,  ldas_index->nr);     
  free_2Darray(&pmask_grid, pmask_nr);
  free_2Darray(&pweight_grid, pweight_nr);
     
}

/* 
   Greg O'Donnell [tempgd@hydro.washington.edu], May, 2000
   
   allocate memory for a 2D array
   array[nrows][ncols]
*/
void alloc_2Darray( int ***array, int ncols, int nrows)
{

  const char *Err = "\nCannot allocate memory in alloc_2Darray.\n";

  int i;

  if( !(*array = (int **) calloc(nrows,sizeof(int *))) ){
    fprintf(stderr,Err);
    fprintf(stderr,"%d %d\n", ncols, nrows);
    exit(EXIT_FAILURE);
  }

  for(i=0;i<nrows;i++){
    if( !((*array)[i] = (int *) calloc(ncols, sizeof(int))) ){
      fprintf(stderr,Err);
      exit(EXIT_FAILURE);
    }
  }

}


/* 
   Greg O'Donnell [tempgd@hydro.washington.edu], May, 2000
   
   free memory for a 2D array
   array[nrows][ncols]
*/
void free_2Darray( int ***array,  int nrows)
{

  int i;

  for(i=0;i<nrows;i++)
    free((*array)[i]);

  free(*array);

}
