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

float estimate_T1(float Ts, 
		   float T1_old,
		   float T2,
		   float D1, 
		   float D2, 
		   float kappa1, 
		   float kappa2, 
		   float Cs1, 
		   float Cs2, 
		   float dp,
		   float delta_t) {
/**********************************************************************
  estimate_T1                Keith Cherkauer          July 15, 1998

  uses Xu Liangs 3-layer energy balance formulation to estimate the 
  temperature between the first and second layers.  Formerly calculated
  independently in each of the surface energy balance equation routines.

  Modifications:
  01-20-00 removed from end of func_surf_energy_bal.c and put into a
           separate file                                           KAC

**********************************************************************/

  float C1;
  float C2;
  float C3;
  float T1;

  C1 = Cs2 * dp / D2 * ( 1. - exp(-D2/dp));
  C2 = - ( 1. - exp(D1/dp) ) * exp(-D2/dp);
  C3 = kappa1/D1 - kappa2/D1 + kappa2/D1*exp(-D1/dp);

  T1 = (kappa1/2./D1/D2*(Ts) + C1/delta_t*T1_old
     + (2.*C2-1.+exp(-D1/dp))*kappa2/2./D1/D2*T2)
     / (C1/delta_t + kappa2/D1/D2*C2 + C3/2./D2);

  return(T1);

}
