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
#include "ftn.h"
//BOP
//
// !ROUTINE: vic_totinit
//
// !DESCRIPTION: 
//  Reinitializes VIC variables before land model run
//
// !INTERFACE:
void FTN(vic_totinit)()
//EOP
{
  extern out_data_struct *outdata1; 
  extern par_struct par; 
  extern ctile_spmd tspmd; 
  int i; 
  //BOC
  for(i=0; i<tspmd.cdi_array[par.rank]; i++){
    outdata1[i].moist_prev = outdata1[i].moist[0] +
                             outdata1[i].moist[1] +
                             outdata1[i].moist[2];
    outdata1[i].swe_prev = outdata1[i].swe; 
    outdata1[i].swnet = 0; 
    outdata1[i].lwnet = 0;
    outdata1[i].qle = 0;
    outdata1[i].qh = 0;
    outdata1[i].qg = 0;
    outdata1[i].rainf = 0;
    outdata1[i].snowf = 0;
    outdata1[i].evap = 0;
    outdata1[i].qs = 0;
    outdata1[i].qsb = 0;
    outdata1[i].qfz = 0;
    outdata1[i].tveg = 0;
    outdata1[i].esoil = 0;
    outdata1[i].qsm = 0;
    outdata1[i].count = 0; 
  }
  //EOC
}
