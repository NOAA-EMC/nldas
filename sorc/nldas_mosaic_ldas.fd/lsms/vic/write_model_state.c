//-------------------------------------------------------------------------
// NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
// Released October 2005
//
// See SOFTWARE DISTRIBUTION POLICY for software distribution policies
//
// The LIS source code and documentation are in the public domain,
// available without fee for educational, research, non-commercial and
// commercial purposes.  Users may distribute the binary or source
// code to third parties provided this statement appears on all copies and
// that no charge is made for such copies.
//
// NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
// SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
// IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
// LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
//
// See COPYRIGHT.TXT for copyright details.
//
//-------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"

//#if SAVE_STATE

void FTN(write_model_state)(char *filename, int *len, int *ntiles, int *Nbands, int *Nlayer, int *Nnode)
/*********************************************************************
  write_model_state      Keith Cherkauer           April 14, 2000

  This subroutine saves the model state at hour 0 of the date 
  defined in the global control file using STATEDAY, STATEMONTH,
  and STATEYEAR.  The saved files can then be used to initialize 
  the model to the same state as when the files were created.

  Soil moisture, soil thermal, and snowpack variables  are stored 
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

*********************************************************************/
{
  extern model_state_struct *model_state;
  extern soil_con_struct *soil_con;

  char *realfile;
  FILE *fp;

  float tmpval;
  float Nsum;
  int   Nveg;
  int    veg;
  int    lidx;
  int    nidx;
  int    Ndist;
  int i, band;

//<kluge>
  model_state_struct * ptr_model_state;
#if ( defined OPENDAP ) && ( defined SPMD )
  extern par_struct par; 
  extern ctile_spmd tspmd;
  int nch;
  if ( par.npes > 1 )
  {
    MPI_Reduce(ntiles,&nch,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

    if ( par.rank == 0 )
    {
       ptr_model_state = (model_state_struct *) lis_calloc(nch, 
                                                    sizeof(model_state_struct),
                                                    "write_model_state"); 
    }
    MPI_Gatherv(model_state,tspmd.cdi_array[par.rank],
                par.MPI_MS_STRUCT,ptr_model_state,tspmd.cdi_array,tspmd.cdispls,
                par.MPI_MS_STRUCT,0,MPI_COMM_WORLD);
  }
  else
  {
    ptr_model_state = model_state;
  }
  if ( par.rank == 0 )
  {
#else
  ptr_model_state = model_state;
#endif
//</kluge>

  printf("write_initial_state\n");

  /* single tile version of the model */
  Nveg = 1;
  
  realfile = (char *)lis_malloc(*len+1, "write_model_state");
  i = 0;
  while(i < *len) {
    realfile[i] = filename[i];
    i++;
  }
  realfile[*len] = '\0';

  printf("FILENAME>>> %s \n", realfile);
  printf("Number of vegs, bands, layers, nodes>>> %d %d %d %d\n", Nveg, *Nbands, *Nlayer, *Nnode);
  fp = open_file(realfile, "wb");
  printf("Successfully opened the file..\n");
  free(realfile);
	    
#if 0
  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;
#else
  Ndist = 1;
#endif

  for ( i = 0; i < *ntiles; i++ ) 
  {
  
//    printf("write_initial_state: tile = %d\n", i);

    /* write cell information */
    fwrite(&Nveg, sizeof(int), 1, fp);
    fwrite(Nlayer, sizeof(int), 1, fp);
    fwrite(Nnode, sizeof(int), 1, fp);

#if 0
    /* Write soil thermal node depths */
    Nsum = 0;
    for ( nidx = 0; nidx < *Nnode; nidx++ )
    {
      fwrite(&Nsum, sizeof(float), 1, fp);
      if ( nidx < *Nnode - 1 )
      {
        Nsum += (soil_con->dz_node[nidx] + soil_con->dz_node[nidx+1]) / 2.;
      }
    }    
#else
    /* Write soil thermal node depths */
    Nsum = 0;
    for ( nidx = 0; nidx < *Nnode; nidx++ )
    {
      fwrite(&soil_con->dz_node[nidx], sizeof(float), 1, fp);
    }    
#endif

    /* Output for all vegetation types */
    for ( veg = 0; veg <= Nveg; veg++ ) {
//    for ( veg = 0; veg < Nveg; veg++ ) {

      /* Output for all snow bands */
      for ( band = 0; band < *Nbands; band++ ) {
      
        /* Write cell identification information */
        fwrite(&veg, sizeof(int), 1, fp);
        fwrite(&band, sizeof(int), 1, fp);
      
        /* Write average total soil moisture */
        for ( lidx = 0; lidx < *Nlayer; lidx++ ) {
#if 0
	  tmpval = prcp->mu[veg] * (cell[WET][veg][band].layer[lidx].moist)
	    + (1. - prcp->mu[veg]) * (cell[DRY][veg][band].layer[lidx].moist);
#else
	  tmpval = ptr_model_state[i].moist[WET][veg][band][lidx];
#endif
          fwrite(&tmpval, sizeof(float), 1, fp);
        }
      
        /* Write average ice content */
        for ( lidx = 0; lidx < *Nlayer; lidx++ ) {
#if 0
	  tmpval = prcp->mu[veg] * (cell[WET][veg][band].layer[lidx].ice)
	    + (1. - prcp->mu[veg]) * (cell[DRY][veg][band].layer[lidx].ice);
#else
          tmpval = ptr_model_state[i].ice[WET][veg][band][lidx];
#endif
          fwrite(&tmpval, sizeof(float), 1, fp);
        }
      
        /* Write average dew storage */
        if ( veg < Nveg ) {
#if 0
	  tmpval = prcp->mu[veg] * (veg_var[WET][veg][band].Wdew)
	    + (1. - prcp->mu[veg]) * (veg_var[DRY][veg][band].Wdew);
#else
          tmpval = ptr_model_state[i].Wdew[WET][veg][band];
#endif
          fwrite(&tmpval, sizeof(float), 1, fp);
        }
      
        /* Write snow data */
        fwrite(&ptr_model_state[i].last_snow[veg][band], sizeof(int), 1, fp);
        fwrite(&ptr_model_state[i].MELTING[veg][band], sizeof(char), 1, fp);
        fwrite(&ptr_model_state[i].coverage[veg][band], sizeof(float), 1, fp);
        fwrite(&ptr_model_state[i].swq[veg][band], sizeof(float), 1, fp);
        fwrite(&ptr_model_state[i].surf_temp[veg][band], sizeof(float), 1, fp);
        fwrite(&ptr_model_state[i].surf_water[veg][band], sizeof(float), 1, fp);
        fwrite(&ptr_model_state[i].pack_temp[veg][band], sizeof(float), 1, fp);
        fwrite(&ptr_model_state[i].pack_water[veg][band], sizeof(float), 1, fp);
        fwrite(&ptr_model_state[i].density[veg][band], sizeof(float), 1, fp);
        fwrite(&ptr_model_state[i].coldcontent[veg][band], sizeof(float),1, fp);
        fwrite(&ptr_model_state[i].snow_canopy[veg][band], sizeof(float),1, fp);
      
        /* Write soil thermal node temperatures */
        for ( nidx = 0; nidx < *Nnode; nidx++ ) { 
          fwrite(&ptr_model_state[i].T[veg][band][nidx], sizeof(float), 1, fp);
        }
      
      }
    }
  }

  /* Force file to be written */
  fflush(fp);
  fclose(fp);

//<kluge>
#if ( defined OPENDAP ) && ( defined SPMD )
   } // close if ( par.rank == 0 )
   if ( par.npes > 1 )
   {
      if ( par.rank == 0 )
      {
         free(ptr_model_state);
      }
   }
#endif
//</kluge>

}

