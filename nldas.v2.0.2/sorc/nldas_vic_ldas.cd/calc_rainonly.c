#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

#define MIN_PREC     1.e-5      /* smallest amount of precipitation that
				   is allowed to fall as snow or rain in
				   a mixed precipitation event */

static char vcid[] = "$Id: calc_rainonly.c,v 3.1 1999/02/16 18:02:07 vicadmin Exp $";

double calc_rainonly(double air_temp,
		     double prec,
		     double MAX_SNOW_TEMP,
		     double MIN_RAIN_TEMP,
		     double mu) {
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
	  
  10-12-00 Logical test of air_temp with MAX_SNOW_TEMP changed from 
           > to >= to handle all cases.				JS

  07-03-01 CHanged test for MAX_SNOW_TEMP<=MIN_RAIN_TEMP to < so that
           we can have snow temp = rain temp (e.g. for 0deg C)
								JS
**********************************************************************/

  double rainonly;

  rainonly = 0.;
  if(MAX_SNOW_TEMP<MIN_RAIN_TEMP)
    vicerror("ERROR: MAX_SNOW_TEMP must be greater then MIN_RAIN_TEMP");
  if(air_temp < MAX_SNOW_TEMP && air_temp > MIN_RAIN_TEMP) {
    rainonly = (air_temp - MIN_RAIN_TEMP)
        / (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
  }
/** Start of changes -- JS **/
/* must be >= and not just >, otherwise it fails if air_temp = MAX_SNOW_TEMP because rainfall = 0 and snowfall = 100% of precip */
  else if(air_temp >= MAX_SNOW_TEMP) {
/** End of changes -- JS **/
    rainonly = prec;
  }

  if(mu < 1) rainonly = prec;

  return(rainonly);

}

#undef MIN_RAIN
