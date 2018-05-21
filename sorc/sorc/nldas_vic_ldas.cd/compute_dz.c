#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id: compute_dz.c,v 4.2 2000/05/16 21:57:54 vicadmin Exp vicadmin $";

void compute_dz(double *dz, double *thermdepths, int Nnodes, double dp) {
/***********************************************************************
  compute_dz              Keith Cherkauer	 	May 16, 2000

  This routines computes the soil thermal node thicknesses using an
  array of node depths.
  
  Modifications:
  11/13/00	Added some lee-way to the thermal depth check to handle
  		floating point representation errors in the data read in
  		from the STATE file.								JS
  12/19/00	Changed all uses of rint function to (int) casts

***********************************************************************/
 
  char ErrStr[MAXSTRING];
  int  j;
  double sum;

/** Start of changes -- JS **/
/* avoid this kind of comparison because of floating point errors in the state file */
/*  if(thermdepths[Nnodes-1] != dp) {*/
/*  if(rint(thermdepths[Nnodes-1]*1000) != rint(dp*1000)) {*/
  if((int)(thermdepths[Nnodes-1]*1000 + 0.5) != (int)(dp*1000 + 0.5)) {
/** End of changes -- JS **/
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %f, but is equal to %f",
	    Nnodes-1,dp,thermdepths[Nnodes-1]);
    nrerror(ErrStr);
  }

  for(j=Nnodes-1;j>0;j--) {
    thermdepths[j] -= thermdepths[j-1];
/** Start of changes -- JS **/
/*    thermdepths[j] = rint(thermdepths[j] * 10000.) / 10000.;*/
    thermdepths[j] = (int)(thermdepths[j] * 10000. + 0.5) / 10000.;
/** End of changes -- JS **/
  }

  sum = 0;
/** Start of changes -- JS **/
/*  dz[0] = rint(thermdepths[1] * 10000.) / 10000.;*/
  dz[0] = (int)(thermdepths[1] * 10000. + 0.5) / 10000.;
/** End of changes -- JS **/
  for(j=1;j<Nnodes;j++) {
/** Start of changes -- JS **/
/*    dz[j] = 2. * rint((thermdepths[j] - dz[j-1] / 2.) * 10000.) / 10000.;*/
    dz[j] = 2. * (int)((thermdepths[j] - dz[j-1] / 2.) * 10000. + 0.5) / 10000.;
/** End of changes -- JS **/
    if(dz[j] < 0) {
      sprintf(ErrStr,"Check spacing between thermal layers %i and %i\n",
	      j-1,j);
      nrerror(ErrStr);
    }
    sum += (dz[j-1] + dz[j]) / 2.;
  }

/** Start of changes -- JS **/
/*  if(rint(sum*1000) != rint(dp*1000)) {*/
  if((int)(sum*1000 + 0.5) != (int)(dp*1000 + 0.5)) {
/** End of changes -- JS **/
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %f, but is equal to %f",
	    Nnodes-1,dp,sum);
    nrerror(ErrStr);
  }

}
