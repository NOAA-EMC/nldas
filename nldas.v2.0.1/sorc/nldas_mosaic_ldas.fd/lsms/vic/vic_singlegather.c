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
#include <vicNl.h>
#include <math.h>
#include "misc.h"
#if (defined SPMD)
#include <mpi.h>
#endif
#include "ftn.h"
extern void FTN(get_vicvar)(int *, int *, float *);
//BOP
//
// !ROUTINE: vic_singlegather.c
//
// !DESCRIPTION: 
//  Gathers each VIC variable to write output
//
// !INTERFACE:
void FTN(vic_singlegather)(int *index, float *tmp){
  //EOP
  extern par_struct par; 
  extern ctile_spmd tspmd; 
  extern out_data_struct *outdata1;
  extern atmos_data_struct *atmos;
  extern dmy_struct dmy; 
  int i; 
  /* float var_temp[tspmd.cdi_array[par.rank]];  */
  float * var_temp;

  var_temp = (float *) lis_calloc(tspmd.cdi_array[par.rank], sizeof(float),
                                  "vic_singlegather");

  FTN(get_vicvar)(index,&(tspmd.cdi_array[par.rank]),var_temp);

  //BOC
#if 0
  if ( *index == 1 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].swnet/outdata1[i].count;
    }
  }
  else if ( *index == 2 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].lwnet/outdata1[i].count;
    }
  }
  else if ( *index == 3 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].qle/outdata1[i].count;
    }
  }
  else if ( *index == 4 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = -1.*outdata1[i].qh/outdata1[i].count;
    }
  }
  else if ( *index == 5 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = -1.*outdata1[i].qg/outdata1[i].count;
    }
  }
  else if ( *index == 6 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = atmos[i].out_prec/(dmy.dt);
    } 
  } 
  else if ( *index == 7 ) 
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].snowf/(dmy.dt);
    }
  }
  else if ( *index == 8 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].evap/(dmy.dt);
    }
  }
  else if ( *index == 9 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].qs;
    }
  }
  else if ( *index == 10 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].qsb;
    }
  }
  else if ( *index == 11 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].qfz;
    }
  }
  else if ( *index == 12 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].snowt+KELVIN;
    }
  }
  else if ( *index == 13 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].avgsurft;
    }
  }
  else if ( *index == 14 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].radt;
    }
  }
  else if ( *index == 15 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].albedo;
    }
  }
  else if ( *index == 16 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].moist[0];
    }
  }
  else if ( *index == 17 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].moist[1];
    }
  }
  else if ( *index == 18 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].moist[2];
    }
  }
  else if ( *index == 19 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].tveg/outdata1[i].count;
    }
  }
  else if ( *index == 20 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].esoil/outdata1[i].count;
    }
  }
  else if ( *index == 21 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].soilwet/outdata1[i].count;
    }
  }
  else if ( *index == 22 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].rootmoist;
    }
  }
  else if ( *index == 23 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].swe*1000.0;
    }
  }
  else if ( *index == 24 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      var_temp[i] = outdata1[i].qsm;
    }
  }
  else if ( *index == 27 )
  {
    for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
    {
      if ( outdata1[i].acond != 0.0 )
      {
        var_temp[i] = 1.0/outdata1[i].acond;
      }
    else
      {
        var_temp[i] = 0.0;
      }
    }
  }
#endif
#if  ( defined SPMD )
  MPI_Gatherv(var_temp,tspmd.cdi_array[par.rank],
              MPI_FLOAT,tmp,tspmd.cdi_array,tspmd.cdispls,
              MPI_FLOAT,0,MPI_COMM_WORLD);
#else
   for ( i = 0; i < tspmd.cdi_array[par.rank]; ++i )
   {
      tmp[i] = var_temp[i];
   }
#endif
  free((void *)var_temp);
  //EOC
}
