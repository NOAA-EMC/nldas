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
#include "ftn.h"

void  FTN(read_vegparam)(int *veg)
{
  extern par_struct par; 
  extern ctile_spmd tspmd; 
  extern veg_con_struct **veg_con; 
  extern veg_lib_struct *veg_lib;
  int i,j,k;
  int type; 
  int vegetat_type_num = 1; 
 
  for(i=0; i<tspmd.cdi_array[par.rank]; i++){
    veg_con[i] = (veg_con_struct *)lis_calloc(vegetat_type_num,
                                              sizeof(veg_con_struct),
                                              "read_vegparam");
  }
  for(i=0; i<tspmd.cdi_array[par.rank]; i++){
    type = veg[i];
    //printf("types ..%d %d \n",i,type);
    for(j=0; j<vegetat_type_num; j++){
      //veg_con[i][j].Cv = 1.0; 
      //veg_con[i][j].Cv_sum = 1.0;
      for(k=0; k<MAX_LAYERS; k++)
	veg_con[i][j].root[k] = 0.0; 
#if 1
      if(type==12 || type==13) {
	veg_con[i][j].zone_depth[0] = 0;
	veg_con[i][j].zone_fract[0] = 0;
	veg_con[i][j].zone_depth[1] = 0;
	veg_con[i][j].zone_fract[1] = 0;
      }
      else{
	veg_con[i][j].zone_depth[0] = veg_lib[type-1].root_depth[0];
	veg_con[i][j].zone_fract[0] = veg_lib[type-1].root_frac[0];
	veg_con[i][j].zone_depth[1] = veg_lib[type-1].root_depth[1];
	veg_con[i][j].zone_fract[1] = veg_lib[type-1].root_frac[1];
      }
#else
      veg_con[i][j].zone_depth[0] = 0.5; 
      veg_con[i][j].zone_fract[0] = 0.5;      
      veg_con[i][j].zone_depth[1] = 0.5;
      veg_con[i][j].zone_fract[1] = 0.5;
#endif
      /**
      if(veg_con[i][j].zone_depth[0]==0.0 ||
	 veg_con[i][j].zone_depth[1]==0.0 ) printf("ERROR>> %d %f %f\n",type,veg_lib[type].root_depth[0],
						   veg_lib[type].root_depth[1]);
      */
      if(type==12 || type==13) 
	veg_con[i][j].veg_class = 0 ;
      else 
	veg_con[i][j].veg_class = type - 1;
      // veg_con[i][j].vegetat_type_num = 1; 
    }
  }
}




