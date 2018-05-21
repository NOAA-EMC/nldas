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

void compute_dz(float *dz, float *thermdepths, int Nnodes, float dp) {
/***********************************************************************
  compute_dz              Keith Cherkauer	 	May 16, 2000

  This routines computes the soil thermal node thicknesses using an
  array of node depths.

***********************************************************************/
 
  char ErrStr[MAXSTRING];
  int  j;
  float sum;
    for(j=Nnodes-1;j>0;j--) {
      printf("node %d, depth = %lf\n", j, thermdepths[j]);
    }
	    
  if((int)(thermdepths[Nnodes-1]*1000 + 0.5) != (int)(dp*1000 + 0.5)) {
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %f, but is equal to %f",
	    Nnodes-1,dp,thermdepths[Nnodes-1]);
    nrerror(ErrStr);
  }

  for(j=Nnodes-1;j>0;j--) {
    thermdepths[j] -= thermdepths[j-1];
    thermdepths[j] = (float)(int)(thermdepths[j] * 10000. + 0.5) / 10000.;
  }

  sum = 0;
  dz[0] = (float)(int)(thermdepths[1] * 10000. + 0.5) / 10000.;
  for(j=1;j<Nnodes;j++) {
    dz[j] = 2. * (float)(int)((thermdepths[j] - dz[j-1] / 2.) * 10000. + 0.5) / 10000.;
    if(dz[j] < 0) {
      sprintf(ErrStr,"Check spacing between thermal layers %i and %i\n",
	      j-1,j);
      nrerror(ErrStr);
    }
    sum += (dz[j-1] + dz[j]) / 2.;
  }

  if((int)(sum*1000 + 0.5) != (int)(dp*1000 + 0.5)) {
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %f, but is equal to %f",
	    Nnodes-1,dp,sum);
    nrerror(ErrStr);
  }

}
