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
void FTN(vic_delcalcs)(){
  extern par_struct par; 
  extern ctile_spmd tspmd; 
  extern out_data_struct *outdata1; 
  
  int i ; 
  for(i=0; i<tspmd.cdi_array[par.rank]; i++){
    outdata1[i].moist_prev = outdata1[i].moist[0]+
      outdata1[i].moist[1]+outdata1[i].moist[2]; 
    outdata1[i].swe_prev = outdata1[i].swe; 
  }
}
