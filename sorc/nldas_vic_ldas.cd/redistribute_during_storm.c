#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id: redistribute_during_storm.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void redistribute_during_storm(cell_data_struct ***cell,
			       veg_var_struct   ***veg_var,
			       int                 veg,
			       int                 Nveg,
			       int                 rec,
			       double              Wdmax,
			       double              old_mu,
			       double              new_mu,
			       double             *max_moist) {
/**********************************************************************
  redistribute_during_storm.c     Keith Cherkauer     January 13, 1998

  This subroutine redistributes soil moisture when the current storm 
  changes intensity.

  Modified:
  06-24-98 Changed to run on only one vegetation type at a time    KAC
  07-13-98 modified to redistribute Wdew within all defined
           elevation bands                                         KAC
  08-19-99 simplified logic, and added check to make sure soil
           moisture does not exceed maximum soil moisture content  Bart

**********************************************************************/
 
  extern option_struct   options;

  unsigned char error;
  char          ErrorString[MAXSTRING];
  int           layer;
  int           band;
  double        temp_wet;
  double        temp_dry;

  /** Redistribute Soil Moisture **/
  for(layer = 0; layer < options.Nlayer; layer++) {

    for(band = 0; band < options.SNOW_BAND; band++) {

      temp_wet = cell[WET][veg][band].layer[layer].moist;
      temp_dry = cell[DRY][veg][band].layer[layer].moist;
      error = redistribute_moisture_for_storm(&temp_wet, &temp_dry, 
					      max_moist[layer], old_mu, 
					      new_mu);
      if(error) {
	sprintf(ErrorString,"%s: Error in moist accounting %f -> %f record %i\n",
		__FILE__,cell[WET][veg][band].layer[layer].moist*old_mu
		+ cell[DRY][veg][band].layer[layer].moist*(1.-old_mu),
		temp_wet*new_mu+temp_dry*(1.-new_mu),rec);
	vicerror(ErrorString);
      }
      cell[WET][veg][band].layer[layer].moist = temp_wet;
      cell[DRY][veg][band].layer[layer].moist = temp_dry;
      
      temp_wet = cell[WET][veg][band].layer[layer].ice;
      temp_dry = cell[DRY][veg][band].layer[layer].ice;
      error = redistribute_moisture_for_storm(&temp_wet, &temp_dry, 
					      max_moist[layer], old_mu, 
					      new_mu);
      if(error) {
	sprintf(ErrorString,"%s: Error in ice accounting %f -> %f record %i\n",
		__FILE__,cell[WET][veg][band].layer[layer].ice*old_mu
		+ cell[DRY][veg][band].layer[layer].ice*(1.-old_mu),
		temp_wet*new_mu+temp_dry*(1.-new_mu),rec);
	vicerror(ErrorString);
      }
      cell[WET][veg][band].layer[layer].ice = temp_wet;
      cell[DRY][veg][band].layer[layer].ice = temp_dry; 
    }
  }

  /****************************************
    Redistribute Stored Water in Vegetation
  ****************************************/
  if(veg < Nveg) {
    for(band = 0; band < options.SNOW_BAND; band++) {
      temp_wet = veg_var[WET][veg][band].Wdew;
      temp_dry = veg_var[DRY][veg][band].Wdew;
      error = redistribute_moisture_for_storm(&temp_wet, &temp_dry, Wdmax, 
					      old_mu, new_mu);
      if(error) {
	sprintf(ErrorString,"%s: Error in Wdew accounting %f -> %f record %i\n",
		__FILE__, veg_var[WET][veg][band].Wdew * old_mu
		+ veg_var[DRY][veg][band].Wdew * (1. - old_mu),
		temp_wet * new_mu + temp_dry * (1. - new_mu), rec);
	vicerror(ErrorString);
      }
      veg_var[WET][veg][band].Wdew = temp_wet;
      veg_var[DRY][veg][band].Wdew = temp_dry;
    }
  }
}

unsigned char redistribute_moisture_for_storm(double *wet_value,
					      double *dry_value,
					      double  max_moist,
					      double  old_mu,
					      double  new_mu) {
/**********************************************************************
  This subroutine redistributes the given parameter between wet and
  dry cell fractions when the precipitation changes in intensity.

  Modified 08-19-99 to add check that maximum soil moisture is not
    exceeded.                                              Bart
**********************************************************************/

  unsigned char error;
  double old_wet;
  double old_dry;
  double diff1 = 0.;
  double diff2 = 0.;
  double diff3 = 0.;

  if (fabs(*wet_value - *dry_value) < SMALL)
    return FALSE;

  old_wet = *wet_value;
  old_dry = *dry_value;

  if (old_mu > new_mu && (1.-new_mu) > SMALL && new_mu > SMALL) {
    *dry_value = (old_mu-new_mu) * old_wet + (1.-old_mu) * old_dry;
    *dry_value /= (1.-new_mu);
  }
  else if ((1.-new_mu) > SMALL && new_mu > SMALL) {
    *wet_value = (new_mu-old_mu) * old_dry + old_mu * old_wet;
    *wet_value /= new_mu;
  }
  else {
    *wet_value = (1.-old_mu) * old_dry + old_mu * old_wet;
    *dry_value = *wet_value;
  }
  diff1 = fabs((old_mu * old_wet + (1.-old_mu) * old_dry) -
    (new_mu * *wet_value + (1.-new_mu) * *dry_value));
  if (*dry_value > max_moist) {
    diff2 = fabs(*dry_value - max_moist);
    *dry_value = max_moist;
  }
  if (*wet_value > max_moist) {
    diff3 = fabs(*wet_value - max_moist);
    *wet_value = max_moist;
  }
/** Start of changes -- JS **/
/* value of SMALL is just too small in this instance, so temporarily change it */
#define SMALL        1.e-5
/** End of changes -- JS **/
  if(diff1 > SMALL || diff2 > SMALL || diff3 > SMALL) {
    fprintf(stderr, "diffs: %f %f %f, SMALL = %f\n", diff1, diff2, diff3, SMALL);
    error = TRUE;
    }
  else 
    error = FALSE;
#define SMALL        1.e-12
  
  return (error);
}
