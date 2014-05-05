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


float linear_interp(float x,float lx,float ux,float ly,float uy) {
  return((x-lx)/(ux-lx)*(uy-ly)+ly);
}

float exp_interp(float x,float lx,float ux,float ly,float uy) {
  return(uy+(ly-uy)*exp(-(x-lx)));
}

float modify_Ksat(float Temp, int frozen_soil) {
/**********************************************************************
	modify_Ksat	Keith Cherkauer		February 12, 1997

  This subroutine returns a parameter to multiply with Ksat to modify
  it for the effects of temperature on the viscosity and density of 
  water.  It is assumed that the given Ksat value was measured at 
  20C (68F).

  Viscosity and density taken from Linsley, "Hydrology for Engineers", 
      A-10

	Temp	Rho	Mu	Factor
	C	kg/m^3	mPa-s	
	0	999.84	1.79	0.560
	5	999.96	1.52	0.659
	10	999.70	1.31	0.770
	15	999.10	1.14	0.878
	20	998.21	1.00	1.00
	25	997.05	0.890	1.12
	30	995.65	0.798	1.25
	35	994.04	0.719	1.39
	40	992.22	0.653	1.52

**********************************************************************/


  float Factor;

  /** formula generated by multiple regression against kinematic
      viscosity data from the Handbook of Chemistry and Physics **/
  if(frozen_soil) {
    Factor = 0.003557 / (0.006534 - 0.0002282 * Temp + 4.794e-6 * (Temp) 
			 * (Temp) - 4.143e-8 * (Temp) * (Temp) * (Temp));
  }
  else Factor = 1.;

  if(Factor>2.) Factor=2.;

  /*return (Factor);*/
  return (1.0);

}