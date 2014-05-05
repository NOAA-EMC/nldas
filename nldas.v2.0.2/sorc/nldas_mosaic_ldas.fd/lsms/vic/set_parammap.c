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
#include "vicNl.h"
#include "ftn.h"

#define SMALL_NUMBER ( 1e-06 )

void FTN(set_parammap)(float *Ds,
                       float *Dsmax,
                       float *Ws,
                       float *b_infilt,
                       float *depth1,
                       float *depth2,
                       float *depth3,
                       int   *nch,
                       int   *nlayer)
{

  extern soil_con_struct *soil_con; 
  int i, j; 

  const float udef = MISSING;
  printf("DBG: set_parammap -- nch = %d\n",*nch);
  printf("DBG: set_parammap -- nlayer = %d\n",*nlayer);

  for ( i = 0; i < *nch; i++)
  {
    if ( Ds[i] > udef + SMALL_NUMBER || 
         Ds[i] < udef - SMALL_NUMBER )
    {
      soil_con[i].Ds = Ds[i];
    }
    if ( Dsmax[i] > udef + SMALL_NUMBER || 
         Dsmax[i] < udef - SMALL_NUMBER )
    {
      soil_con[i].Dsmax = Dsmax[i]; 
    }
    if ( Ws[i] > udef + SMALL_NUMBER || 
         Ws[i] < udef - SMALL_NUMBER )
    {
      soil_con[i].Ws = Ws[i]; 
    }
    if ( b_infilt[i] > udef + SMALL_NUMBER || 
         b_infilt[i] < udef - SMALL_NUMBER )
    {
      soil_con[i].b_infilt = b_infilt[i];
    }
    if ( depth1[i] > udef + SMALL_NUMBER || 
         depth1[i] < udef - SMALL_NUMBER )
    {
      soil_con[i].depth[0] = depth1[i]; 
    }
    if ( depth2[i] > udef + SMALL_NUMBER || 
         depth2[i] < udef - SMALL_NUMBER )
    {
      soil_con[i].depth[1] = depth2[i];
    } 
    if ( depth3[i] > udef + SMALL_NUMBER || 
         depth3[i] < udef - SMALL_NUMBER )
    {
      soil_con[i].depth[2] = depth3[i];       
      //printf("ds %d %f \n", i, soil_con[i].Ds); 
    }
  }


  /* Set initial soil moisture content at critical value. */
  /* This can be overwritten by a parameter map if necessary. */
  for ( i = 0; i < *nch; i++)
  {
    for ( j = 0; j < *nlayer; j++)
    {
      soil_con[i].init_moist[j] = 1000.0 * soil_con[i].Wcr[j]
                                         * soil_con[i].porosity[j]
                                         * soil_con[i].depth[j];
      //      printf("init moist %d %f \n",i,soil_con[i].depth[j]); 

      soil_con[i].max_moist[j]  = 1000.0 * soil_con[i].depth[j] 
                                         * soil_con[i].porosity[j];
      if ( soil_con[i].init_moist[j] > soil_con[i].max_moist[j] )
      {
        soil_con[i].init_moist[j] = soil_con[i].max_moist[j];
      }                               

      //printf("%d Wcr = %f, poro = %f, depth = %f, max = %f, init = %f\n", 
      //       j, soil_con[i].Wcr[j], soil_con[i].porosity[j], 
      //       soil_con[i].depth[j], soil_con[i].max_moist[j], 
      //       soil_con[i].init_moist[j]);

    }
  }

}
 

