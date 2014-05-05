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
#include "vicNl.h"
#include "ftn.h"

void FTN(read_veglib)(char *filename, int *len, int *Ntype, int *Nroot)
/**********************************************************************
  read_veglib.c               Keith Cherkauer                 1997

  This routine reads in a library of vegetation parameters for all
  vegetation classes used in the model.  The veg class number is used
  to reference the information in this library.

  Modifications:
  09-24-98 Modified to remove root fractions from the library file.
           See read_vegparam.c and calc_root_fraction.c for new
           root fraction distribution information.               KAC
	   
**********************************************************************/
{
  extern FILE *open_file(char string[], char type[]);
  extern veg_lib_struct *veg_lib;
  char *realfile; 
  FILE *veglib; 
  
  int    i, j;
  int    tmpflag;
  int    Nveg_type;
  char   str[MAXSTRING];
  char   ErrStr[MAXSTRING];
  float maxd;

  float sum, depth_sum; 
  realfile = (char *) lis_malloc(*len+1,"read_veglib"); 
  i = 0; 
  while(i<*len) {
    realfile[i] = filename[i]; 
    i++; 
  }
  realfile[*len] = '\0'; 
  printf("veg lib filename >> %s \n",realfile); 
  veglib = open_file(realfile, "r");
  printf("Successfully opened veg lib file \n");
  free(realfile);

  rewind(veglib);
  fgets(str,MAXSTRING,veglib);
  Nveg_type = 0;
  while(!feof(veglib)) {
    if(str[0]<=57 && str[0]>=48) Nveg_type++;
    fgets(str,MAXSTRING,veglib);
  }
  rewind(veglib);
  veg_lib = (veg_lib_struct *)lis_calloc(Nveg_type,
                                         sizeof(veg_lib_struct),
                                         "read_veglib");

  fscanf(veglib, "%s", str);
  i=0;
  while (!feof(veglib)) {
    if(str[0]<=57 && str[0]>=48) {
      veg_lib[i].veg_class = atoi(str);
      sprintf(ErrStr, "veg lib class %d %d \n",i,veg_lib[i].veg_class);
      lis_log_msgC(ErrStr);
      fscanf(veglib, "%i",  &tmpflag);
      if(tmpflag==0) veg_lib[i].overstory = FALSE;
      else veg_lib[i].overstory = TRUE;
      fscanf(veglib, "%f", &veg_lib[i].rarc);
      fscanf(veglib, "%f", &veg_lib[i].rmin);
      //printf("RARC ..%d %f %f\n",i,veg_lib[i].rarc, veg_lib[i].rmin);
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%f", &veg_lib[i].LAI[j]);
        veg_lib[i].Wdmax[j] = LAI_WATER_FACTOR * veg_lib[i].LAI[j];
      }
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%f", &veg_lib[i].albedo[j]);
      }
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%f", &veg_lib[i].roughness[j]);
      }
      veg_lib[i].wind_h = 0.;
      maxd = 0;
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%f", &veg_lib[i].displacement[j]);
        if(veg_lib[i].displacement[j] > maxd) maxd = veg_lib[i].displacement[j];
        if(veg_lib[i].LAI[j] > 0 && veg_lib[i].displacement[j] <= 0) {
          sprintf(str,"Vegetation has leaves (LAI = %f), but no displacement (%f)",
	          veg_lib[i].LAI[j], veg_lib[i].displacement[j]);
          nrerror(str);
        }
        if(veg_lib[i].albedo[j] < 0 || veg_lib[i].albedo[j] > 1) {
          sprintf(str,"Albedo must be between 0 and 1 (%f)",
	          veg_lib[i].albedo[j]);
          nrerror(str);
        }
      }
      fscanf(veglib, "%f", &veg_lib[i].wind_h);
      if(veg_lib[i].wind_h < maxd && veg_lib[i].overstory) {
        sprintf(str,"Vegetation reference height (%f) for vegetation class %i, must be greater than the maximum displacement height (%f) when OVERSTORY has been set TRUE.",
                veg_lib[i].wind_h,veg_lib[i].veg_class,maxd);
        nrerror(str);
      }
      fscanf(veglib, "%f",  &veg_lib[i].RGL);         /* minimum value of incoming
						    solar radiation at which there
						   will still be transpiration */
      if(veg_lib[i].RGL < 0) {
        sprintf(str,"Minimum value of incoming solar radiation at which there is transpiration (RGL) must be greater than 0 for vegetation class %i.  Check that the vegetation library has the correct number of columns.",
                veg_lib[i].veg_class);
        nrerror(str);
      }
      fscanf(veglib, "%f", &veg_lib[i].rad_atten);   /* vegetation radiation 
						      attenuation factor */
      if(veg_lib[i].rad_atten < 0 || veg_lib[i].rad_atten > 1) {
        sprintf(str,"The vegetation radiation attenuation factor must be greater than 0, and less than 1 for vegetation class %i.  Check that the vegetation library has the correct number of columns.",
                veg_lib[i].veg_class);
        nrerror(str);
      }
      fscanf(veglib, "%f", &veg_lib[i].wind_atten);  /* canopy wind speed
						      attenuation factor */
      fscanf(veglib, "%f", &veg_lib[i].trunk_ratio); /* ratio of tree height that
						      is trunk */
      /* start of changes -- JS */
      /* read in the rooting depths and fractions */
      depth_sum = 0.0;
      sum = 0.0;
      for(j=0; j<(*Nroot); j++) {
        fscanf(veglib, "%f %f", &veg_lib[i].root_depth[j], &veg_lib[i].root_frac[j]);
       	//printf("ROOT depths.. %d %d %f %f \n",i,j,veg_lib[i].root_depth[j], veg_lib[i].root_frac[j]);
        depth_sum += veg_lib[i].root_depth[j];
        sum += veg_lib[i].root_frac[j];	    
      }
      if(depth_sum <= 0) {
        sprintf(str,"Root zone depths must sum to a value greater than 0.");
        nrerror(str);
      }
      if(fabs((sum - 1.0)>1e-5)) {
        fprintf(stderr,"WARNING: Root zone fractions sum to more than 1 ( = %f), normalizing fractions.  If the sum is large, check that your vegetation parameter file is in the form - <zone 1 depth> <zone 1 fract> <zone 2 depth> <zone 2 fract> ...\n", sum);
        for(j=0; j<*Nroot; j++) {
          veg_lib[i].root_frac[j] /= sum;
        }
      }
      /* end of changes -- JS */
      fgets(str, MAXSTRING, veglib);	/* skip over end of line comments */
      i++;
    }
    else fgets(str, MAXSTRING, veglib);
    fscanf(veglib, "%s", str);
  }
  if(i!=Nveg_type) {
    sprintf(ErrStr,"ERROR: Problem reading vegetation library file - make sure the file has the right number of columns.\n");
    nrerror(ErrStr);
  }
  *Ntype = Nveg_type;
  for(i=0; i<Nveg_type; i++)
  {
    sprintf(ErrStr, "veg lib ..%d %d \n",i,veg_lib[i].veg_class);
    lis_log_msgC(ErrStr);
  }
  fclose(veglib);
} 





