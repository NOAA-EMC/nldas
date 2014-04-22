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

#define MIN_PREC     1.e-5      /* smallest amount of precipitation that
				   is allowed to fall as snow or rain in
				   a mixed precipitation event */


float calc_rainonly(float air_temp,
		     float prec,
		     float mu) {
/**********************************************************************
  calc_rainonly.c	Keith Cherkauer		March 7, 1998

  Determines from the air temperature what fraction of incoming
  precipitation is frozen and unfrozen (snow and rain).

  Modifications:
  09-22-98 Modified to filter out very small fractions of snow
           or rain in mixed precipitation.  Minimum value MIN_PREC
	   is adjusted to account for the size of mu (minimum
	   is based of fractional precipitation with mu=1, since
	   snow cannot be solved for when mu<1).                  KAC

**********************************************************************/

  float rainonly;

  rainonly = 0.;
  if(MAX_SNOW_TEMP<=MIN_RAIN_TEMP)
    printf("ERROR: MAX_SNOW_TEMP must be greater than MIN_RAIN_TEMP\n");
    //    vicerror("ERROR: MAX_SNOW_TEMP must be greater then MIN_RAIN_TEMP");
  //printf("temps..%f %f %f \n",air_temp, MAX_SNOW_TEMP, mu);
  if(air_temp < MAX_SNOW_TEMP && air_temp > MIN_RAIN_TEMP) {
    rainonly = (air_temp - MIN_RAIN_TEMP)
        / (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
  }
  else if(air_temp >= MAX_SNOW_TEMP) {
    rainonly = prec;
  }

  if(mu < 1) rainonly = prec;

  return(rainonly);

}

#undef MIN_RAIN
