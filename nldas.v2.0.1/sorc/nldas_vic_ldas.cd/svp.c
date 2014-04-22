#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: svp.c,v 3.2 1999/07/07 16:32:11 vicadmin Exp $";

double svp(double temp)
/**********************************************************************
  This routine computes the saturated vapor pressure using Handbook
  of Hydrology eqn 4.2.2

  Pressure in kPa

  Modifications:
  8/30/00	Function returns svp of zero if the input temperature
  		is below -C_SVP. Otherwise the function fails and
  		returns a non-sensical value. Note: this should have 
  		been fixed in the stabdard version with the use of the
  		constant MIN_TDEW, but I don't think this was implemented.
  														JS

**********************************************************************/
{
  double SVP;
  
/** Start of changes - JS **/
  if(temp < -C_SVP) return (0.0);
/** End of changes - JS **/

  SVP = A_SVP * exp((B_SVP * temp)/(C_SVP+temp));

  if(temp<0) SVP *= 1.0 + .00972 * temp + .000042 * temp * temp;

  return (SVP);
}

double svp_slope(double temp)
/**********************************************************************
  This routine computes the gradient of d(svp)/dT using Handbook
  of Hydrology eqn 4.2.3
**********************************************************************/
{
  return (B_SVP * C_SVP) / ((C_SVP + temp) * (C_SVP + temp)) * svp(temp);
}


 
