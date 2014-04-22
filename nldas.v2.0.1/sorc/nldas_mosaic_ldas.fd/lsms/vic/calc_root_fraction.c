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
#include <math.h>
#include "vicNl.h"
#include "ftn.h"
void FTN(calc_root_fractions)(int *Nlayer,int *Root_zones)
/**********************************************************************
  calc_root_fraction.c    Keith Cherkauer      September 24, 1998

  This routine computes the fraction of roots in each soil layer based
  on the root zone distribution defined in the vegetation parameter
  file.  Roots are assumed to be linearly distributed within each
  root zone.

**********************************************************************/
{
  extern par_struct par; 
  extern ctile_spmd tspmd; 
  extern veg_con_struct **veg_con; 
  extern soil_con_struct *soil_con; 
  /* char   ErrStr[MAXSTRING]; */
  int    Nveg;
  int    veg;
  int    layer;
  int    zone;
  int    i,k;
  float  sum_depth;
  float  sum_fract;
  float  dum;
  float Zstep;
  float Zsum;
  float Lstep;
  float Lsum;
  float Zmin_fract;
  float Zmin_depth;
  float Zmax;

  //  Nveg      = veg_con[0].vegetat_type_num;
  Nveg      = 1;
  for(k=0; k<tspmd.cdi_array[par.rank]; k++){
  for(veg=0;veg<Nveg;veg++) {
    sum_depth  = 0;
    sum_fract  = 0;
    layer      = 0;
    Lstep      = soil_con[k].depth[layer];
    Lsum       = Lstep;
    Zsum       = 0;
    zone       = 0;
    /** mostly for 11 and 12..*/
    for(layer=0; layer<*Nlayer; layer++)
      veg_con[k][veg].root[layer] = 0;
    layer = 0;
    
    while(zone<*Root_zones) {
      //printf("here1 %d %d %d %f %f\n",k,veg,veg_con[k][veg].veg_class, veg_con[k][veg].zone_depth[zone],veg_con[k][veg].zone_fract[zone]);
      Zstep = (float) veg_con[k][veg].zone_depth[zone];
      if((Zsum + Zstep) <= Lsum && Zsum >= Lsum - Lstep) {
	/** CASE 1: Root Zone Completely in Soil Layer **/
	sum_fract += veg_con[k][veg].zone_fract[zone];
      }      
      else {
	/** CASE 2: Root Zone Partially in Soil Layer **/
	if(Zsum < Lsum - Lstep) {
	  /** Root zone starts in previous soil layer **/
	  Zmin_depth = Lsum - Lstep;
	  Zmin_fract = linear_interp(Zmin_depth,Zsum,Zsum+Zstep,0,
				     veg_con[k][veg].zone_fract[zone]);
	}
	else {
	  /** Root zone starts in current soil layer **/
	  Zmin_depth = Zsum;
	  Zmin_fract = 0.;
	}
	if(Zsum + Zstep <= Lsum) {
	  /** Root zone ends in current layer **/
	  Zmax = Zsum + Zstep;
	}
	else {
	  /** Root zone extends beyond bottom of current layer **/
	  Zmax = Lsum;
	}
	sum_fract += linear_interp(Zmax,Zsum,Zsum+Zstep,0,
				   veg_con[k][veg].zone_fract[zone]) - Zmin_fract;
      }

      /** Update Root Zone and Soil Layer **/
      if(Zsum + Zstep < Lsum) {
	Zsum += Zstep;
	zone ++;
      }
      else if(Zsum + Zstep == Lsum) {
	Zsum += Zstep;
	zone ++;
	if(layer<*Nlayer) {
	  veg_con[k][veg].root[layer] = sum_fract;
	  sum_fract = 0.;
	}
	layer++;
	if(layer<*Nlayer) {
	  Lstep  = soil_con[k].depth[layer];
	  Lsum  += Lstep;
	}
	else if(layer==*Nlayer) {
	  Lstep  = Zsum + Zstep - Lsum;
	  if(zone<*Root_zones-1) {
	    for(i=zone+1;i<*Root_zones;i++) {
	      Lstep += veg_con[k][veg].zone_depth[i];
	    }
	  }
	  Lsum  += Lstep;
	}
      }
      else if(Zsum + Zstep > Lsum) {
	if(layer<*Nlayer) {
	  veg_con[k][veg].root[layer] = sum_fract;
	  sum_fract = 0.;

	}
	layer++;
	if(layer<*Nlayer) {
	  Lstep  = soil_con[k].depth[layer];
	  Lsum  += Lstep;
	}
	else if(layer==*Nlayer) {
	  Lstep  = Zsum + Zstep - Lsum;
	  if(zone<*Root_zones-1) {
	    for(i=zone+1;i<*Root_zones;i++) {
	      Lstep += veg_con[k][veg].zone_depth[i];
	    }
	  }
	  Lsum  += Lstep;
	}
      }
	
    }

    if(sum_fract > 0 && layer >= *Nlayer) {
      veg_con[k][veg].root[*Nlayer-1] += sum_fract;
    }
    else if(sum_fract > 0) {
      veg_con[k][veg].root[layer] += sum_fract;
    }

    dum=0.;
    for (layer=0;layer<*Nlayer;layer++) {
      if(veg_con[k][veg].root[layer] < 1.e-4) veg_con[k][veg].root[layer] = 0.;
      dum += veg_con[k][veg].root[layer];
    }
    /**    if(dum == 0.0){
      printf("Crashing on tile %d %f %f\n",k,veg_con[k][veg].zone_fract[0],
	     veg_con[k][veg].zone_fract[1]);
      sprintf(ErrStr,"Root fractions sum equals zero: %f , Vege Class: %d\n",
	      dum, veg_con[k][veg].veg_class);
      nrerror(ErrStr);
      }*/
    if(dum!=0.0){
      for (layer=0;layer<*Nlayer;layer++) {
	veg_con[k][veg].root[layer] /= dum;
      }
    }
    

  }
  }
  /**
  for(k=0; k<tspmd.cdi_array[par.rank]; k++){
      for(veg=0;veg<Nveg;veg++) {
	for (layer=0;layer<*Nlayer;layer++) {
	  printf("ROOT %d %d %d %f \n",par.rank,k,veg,veg_con[k][veg].root[layer]);
	}
      }
  }
  */
}


