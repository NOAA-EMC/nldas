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

void FTN(read_initial_model_state)(char *filename, int *len, int *Nbands, int *Nlayer, int *Nnode)
/*********************************************************************
  read_initial_model_state   Keith Cherkauer         April 14, 2000

  This subroutine initializes the model state at hour 0 of the date 
  defined in the given state file.  

  Soil moisture, soil thermal, and snowpack variables  are stored 
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

*********************************************************************/
{

  extern dist_prcp_struct *prcp;
  extern soil_con_struct *soil_con;
  extern par_struct par;
  extern ctile_spmd tspmd;
    
  FILE *fp;
  char *realfile;
  int i;
 
  char   ErrStr[MAXSTRING];
  int Nveg;
  float depth_node[MAX_NODES];
  float sum;
  int    veg, iveg;
  int    band, iband;
  int    lidx;
  int    nidx;
  int    tmp_Nveg;
  int    tmp_Nlayer;
  int    tmp_Nnode;
  int ntiles;
  
  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  printf("read_initial_state\n");
  
  /* single tile version of the model */
  Nveg = 1;

  realfile = (char *) lis_malloc(*len+1,"read_initial_model_state");
  i = 0;
  while(i < *len) {
    realfile[i] = filename[i];
    i++;
  }
  realfile[*len] = '\0';

  printf("FILENAME>>>%s \n", realfile);
  fp = open_file(realfile, "rb");
  printf("Successfully opened the file..\n");
  free(realfile);
  printf("Number of vegs, layers, nodes>>> %d %d %d \n", Nveg, *Nlayer, *Nnode);

  /* loop over tiles */
  ntiles = tspmd.cdi_array[par.rank];
  for(i=0; i<ntiles; i++) {
  
    //printf("read_initial_state: tile i = %d\n", i);

    cell    = prcp[i].cell;
    veg_var = prcp[i].veg_var;
    snow    = prcp[i].snow;
    energy  = prcp[i].energy;

    fread(&tmp_Nveg, sizeof(int), 1, fp);
    fread(&tmp_Nlayer, sizeof(int), 1, fp);
    fread(&tmp_Nnode, sizeof(int), 1, fp);

    //printf("read nveg = %d, nlayer = %d, nnode = %d\n", tmp_Nveg, tmp_Nlayer, tmp_Nnode);
    if ( tmp_Nveg != Nveg ) {
      sprintf(ErrStr,"The number of vegetation types in cell %i (%i) does not equal that defined in vegetation parameter file (%i).  Check your input files.", i, tmp_Nveg, Nveg);
      nrerror(ErrStr);
    }
    if ( tmp_Nlayer != *Nlayer ) {
      sprintf(ErrStr,"The number of soil layers in cell %i (%i) does not equal that defined in the model (%i).  Check your input files.", i, tmp_Nlayer, 3);
      nrerror(ErrStr);
    }
  
    if ( tmp_Nnode != *Nnode ) {
      sprintf(ErrStr,"The number of nodes in cell %i (%i) does not equal that defined in the model definition (%i).  Check your input files.", i, tmp_Nnode, *Nnode);
      nrerror(ErrStr);
    }

    /* Read soil thermal node depths */
    sum = 0.0;
    for ( nidx = 0; nidx < *Nnode; nidx++ ) { 
      fread(&depth_node[nidx], sizeof(float), 1, fp);
      /* JS: save node depths in soil_con struct */
      /* soil_con->depth_node[nidx] = depth_node[nidx]; */
      soil_con->dz_node[nidx] = depth_node[nidx];
      sum += soil_con->dz_node[nidx];
    }
//    if (*Nnode > 1)
//      compute_dz(soil_con->dz_node, depth_node, *Nnode, soil_con->dp);
//    else soil_con->dz_node[0] = 0;
    if ( *Nnode == 1) soil_con->dz_node[0] = 0;
    sum -= 0.5 * soil_con->dz_node[0];
    sum -= 0.5 * soil_con->dz_node[*Nnode-1];
    /* if ( abs( sum - soil_con->dp ) > SMALL ) */
    if ( sum > soil_con->dp )
    {
      fprintf( stderr, "WARNING: Sum of soil nodes (%f) exceeds defined damping depth (%f).  Resetting damping depth.\n", sum, soil_con->dp );
      soil_con->dp = sum;
    }


    /* Input for all vegetation types */
    for ( veg = 0; veg <= Nveg; veg++ ) {
//    for ( veg = 0; veg < Nveg; veg++ ) {
    
      /* Input for all snow bands */
      for ( band = 0; band < *Nbands; band++ ) {
      
        /* Read cell identification information */
        if(fread(&iveg, sizeof(int), 1, fp) != 1)
          nrerror("End of model state file found unexpectedly 1");
        if(fread(&iband, sizeof(int), 1, fp) != 1)
          nrerror("End of model state file found unexpectedly 2");
        if ( iveg != veg || iband != band ) {
          fprintf(stderr,"The vegetation and snow band indexs in the model state file (veg = %i, band = %i) do not match those currently requested (veg = %i , band = %i).  Model state file must be stored with variables for all vegetation indexed by variables for all snow bands.", iveg, iband, veg, band);
          nrerror(ErrStr);
        }
      
        /* Read average total soil moisture */
        for ( lidx = 0; lidx < *Nlayer; lidx++ ) {
          if(fread(&cell[WET][veg][band].layer[lidx].moist, sizeof(float), 1, fp) == 0)
            nrerror("End of model state file found unexpectedly 3");
        }
      
        /* Read average ice content */
        for ( lidx = 0; lidx < *Nlayer; lidx++ ) {
          if(fread(&cell[WET][veg][band].layer[lidx].ice, sizeof(float), 1, fp) == 0)
            nrerror("End of model state file found unexpectedly 4");
        }
      
        /* Read average dew storage */
        if ( veg < Nveg ) {
          if(fread(&veg_var[WET][veg][band].Wdew, sizeof(float), 1, fp) == 0)
            nrerror("End of model state file found unexpectedly 5");
        }
      
        /* Read snow data */
        if(fread(&snow[veg][band].last_snow, sizeof(int), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 6");
        if(fread(&snow[veg][band].MELTING, sizeof(char), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 7");
        if(fread(&snow[veg][band].coverage, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 8");
        if(fread(&snow[veg][band].swq, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 9");
        if(fread(&snow[veg][band].surf_temp, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 10");
        if(fread(&snow[veg][band].surf_water, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 11");
        if(fread(&snow[veg][band].pack_temp, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 12");
        if(fread(&snow[veg][band].pack_water, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 13");
        if(fread(&snow[veg][band].density, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 14");
        if(fread(&snow[veg][band].coldcontent, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 15");
        if(fread(&snow[veg][band].snow_canopy, sizeof(float), 1, fp) == 0)
          nrerror("End of model state file found unexpectedly 16");
        if(snow[veg][band].density > 0.) 
          snow[veg][band].depth = 1000. * snow[veg][band].swq 
          / snow[veg][band].density;
      
        /* Read soil thermal node temperatures */
        for ( nidx = 0; nidx < *Nnode; nidx++ ) 
          if(fread(&energy[veg][band].T[nidx], sizeof(float), 1, fp) == 0)
            nrerror("End of model state file found unexpectedly 17");
      
      }
    }
  }
  fclose(fp);
}
