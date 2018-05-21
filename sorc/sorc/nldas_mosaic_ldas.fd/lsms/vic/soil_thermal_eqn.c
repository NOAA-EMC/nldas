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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

float soil_thermal_eqn(float T, va_list ap) {

  float value;

  float TL;
  float TU;
  float T0;
  float moist;
  float max_moist;
  float bubble;
  float expt;
  float ice0;
  float gamma;
  float fprime;
  float A;
  float B;
  float C;
  float D;
  float E;
  float ice;

  TL         = (float) va_arg(ap, double);
  TU         = (float) va_arg(ap, double);
  T0         = (float) va_arg(ap, double);
  moist      = (float) va_arg(ap, double);
  max_moist  = (float) va_arg(ap, double);
  bubble     = (float) va_arg(ap, double);
  expt       = (float) va_arg(ap, double);
  ice0       = (float) va_arg(ap, double);
  gamma      = (float) va_arg(ap, double);
  fprime     = (float) va_arg(ap, double);
  A          = (float) va_arg(ap, double);
  B          = (float) va_arg(ap, double);
  C          = (float) va_arg(ap, double);
  D          = (float) va_arg(ap, double);
  E          = (float) va_arg(ap, double);

  if(T<0.) {

    ice = moist - maximum_unfrozen_water(T,max_moist,bubble,expt);
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;
  value = T*E - A*(TL-TU) - B*(TL+TU-gamma*fprime) - C - D*(ice-ice0);

  return(value);

}
