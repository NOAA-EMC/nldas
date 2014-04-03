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

/********************************************************************************
  filename  : penman.c
  purpose   : Calculate daily evapotranspiration using the combination equation
  interface : - input : - net radiation (W/m2)
                        - vapor pressure deficit (Pa)
			- aerodynamical resistance (s/m)
			- minimum stomatal resistance (s/m)
			- architectural resistance (s/m)
			- LAI
			- soil moisture stress factor
			- air temperature (to calculate slope of saturated 
			  vapor pressure curve) (C)
			- elevation (m)
              - output: - daily evapotranspiration (mm/day)
  programmer: Bart Nijssen
  date      : August 7, 1995
  changes   :
  references: 
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>


#define CLOSURE 4000		/** Pa **/
#define RSMAX 5000
#define VPDMINFACTOR 0.1

float vic_penman(float rad, 
	      float vpd, 
	      float ra, 
	      float rs, 
	      float rarc, 
              float lai, 
	      float gsm_inv, 
	      float tair, 
              float net_short, 
	      float  elevation, 
	      float  RGL)
{
  float evap;                  /* Penman-Monteith evapotranspiration */
  float slope;                 /* slope of saturated vapor pressure curve */
  float rc;                    /* canopy resistance */
  float r_air;                 /* density of air in kg/m3 */
  float h;                     /* scale height in the atmosphere (m) */
  float lv;                    /* latent heat of vaporization (J/kg) */
  float pz;                    /* surface air pressure */
  float gamma;                 /* psychrometric constant (Pa/C) */
  float Tfactor;               /* factor for canopy resistance based on temperature */
  float vpdfactor;             /* factor for canopy resistance based on vpd */
  float DAYfactor;             /* factor for canopy resistance based on photosy
nthesis */
  float f;
  //printf("penman %f %f %f %f %f %f %f %f %f %f %f\n",rad,vpd,ra,rs,rarc,lai,gsm_inv,tair,net_short,elevation,RGL);
  /* calculate the slope of the saturated vapor pressure curve in Pa/K */
  slope = svp_slope(tair)*1000;

  /* calculate resistance factors (Wigmosta et al., 1994) */

  if(rs>0.) {
    f = net_short / RGL;
    DAYfactor = (1. + f)/(f + rs/RSMAX);
  }
  else DAYfactor = 1.;

  Tfactor = .08 * tair - 0.0016 * tair * tair;
  Tfactor = (Tfactor <= 0.0) ? 1e-10 : Tfactor;

  vpdfactor = 1 - vpd/CLOSURE;
  vpdfactor = (vpdfactor < VPDMINFACTOR) ? VPDMINFACTOR : vpdfactor;

  /* calculate canopy resistance in s/m */
  rc = rs/(lai * gsm_inv * Tfactor * vpdfactor) * DAYfactor;
  rc = (rc > RSMAX) ? RSMAX : rc;
  //printf("rc %f %f %f %f %f %f \n",rc,rs,lai,gsm_inv, Tfactor, vpdfactor); 
  /* calculate scale height based on average temperature in the column */
  h  = 287/9.81 * ((tair + 273.15) + 0.5 * (float)elevation * LAPSE_PM);

  /* use hypsometric equation to calculate p_z, assume that virtual
     temperature is equal air_temp */
  pz = PS_PM * exp(-(float)elevation/h);
  
  /* calculate latent heat of vaporization. Eq. 4.2.1 in Handbook of 
     Hydrology, assume Ts is Tair */
  lv = 2501000 - 2361 * tair;
  
  /* calculate gamma. Eq. 4.2.28. Handbook of Hydrology */
  gamma = 1628.6 * pz/lv;
  
  /* calculate factor to be applied to rc/ra */
  
  /* calculate the air density, using eq. 4.2.4 Handbook of Hydrology */
  r_air = 0.003486 * pz/(275 + tair);
 
  /* calculate the evaporation in mm/day (by not dividing by the density 
     of water (~1000 kg/m3)), the result ends up being in mm instead of m */ 
    
  evap = (slope * rad + r_air * CP_PM * vpd/ra)/
         (lv * (slope + gamma * (1 + (rc + rarc)/ra))) * SEC_PER_DAY;
  //printf("evap %f %f %f %f \n",evap,ra,rc,rarc);
  if (vpd > 0.0 && evap < 0.0) 
    evap = 0.0;
  //printf("returning from pen..%f \n",evap);
  return evap;
}

float quick_penman(float rad, 
		    float vpd, 
		    float gsm_inv,
		    float CONST_1,
		    float CONST_2,
		    float CONST_3,
		    float CONST_4,
		    float CONST_5)
{
  float evap;			/* Penman-Monteith evapotranspiration */
  float rc;			/* canopy resistance */

  /* calculate canopy resistance in s/m */
  rc = CONST_5 / gsm_inv;
  rc = (rc > RSMAX) ? RSMAX : rc;

  /* calculate the evaporation in mm/day (by not dividing by the density 
     of water (~1000 kg/m3)), the result ends up being in mm instead of m */ 
    
  evap = (CONST_1 * rad + CONST_2) / (CONST_3 + CONST_4 * rc);

/*   if (vpd > 0.0 && evap < 0.0)  */
/*     evap = 0.0; */

  return evap;
}

void compute_penman_constants(float  vpd, 
			      float  ra, 
			      float  rs, 
			      float  rarc, 
			      float  lai, 
			      float  tair, 
			      float  net_short, 
			      float   elevation, 
			      float   RGL,
			      float *CONST_1,
			      float *CONST_2,
			      float *CONST_3,
			      float *CONST_4,
			      float *CONST_5)
{
  float slope;			/* slope of saturated vapor pressure curve */
  float r_air;			/* density of air in kg/m3 */
  float h;			/* scale height in the atmosphere (m) */
  float lv;			/* latent heat of vaporization (J/kg) */
  float pz;			/* surface air pressure */
  float gamma;			/* psychrometric constant (Pa/C) */
  float Tfactor;		/* factor for canopy resistance based on 
				   temperature */
  float vpdfactor;		/* factor for canopy resistance based on vpd */
  float DAYfactor;		/* factor for canopy resistance based on 
				   photosynthesis */
  float f;

  /* calculate the slope of the saturated vapor pressure curve in Pa/K */
  slope = svp_slope(tair)*1000;

  /* calculate resistance factors (Wigmosta et al., 1994) */

  if(rs>0.) {
    f = net_short / RGL;
    DAYfactor = (1. + f)/(f + rs/RSMAX);
  }
  else DAYfactor = 1.;
  
  Tfactor = .08 * tair - 0.0016 * tair * tair;
  Tfactor = (Tfactor <= 0.0) ? 1e-10 : Tfactor;
  
  vpdfactor = 1 - vpd/CLOSURE;
  vpdfactor = (vpdfactor < VPDMINFACTOR) ? VPDMINFACTOR : vpdfactor;
  
  /* calculate scale height based on average temperature in the column */
  h  = 287/9.81 * ((tair + 273.15) + 0.5 * (float)elevation * LAPSE_PM);
  
  /* use hypsometric equation to calculate p_z, assume that virtual
     temperature is equal air_temp */
  pz = PS_PM * exp(-(float)elevation/h);
  
  /* calculate latent heat of vaporization. Eq. 4.2.1 in Handbook of 
     Hydrology, assume Ts is Tair */
  lv = 2501000 - 2361 * tair;
  
  /* calculate gamma. Eq. 4.2.28. Handbook of Hydrology */
  gamma = 1628.6 * pz/lv;
  
  /* calculate factor to be applied to rc/ra */
  
  /* calculate the air density, using eq. 4.2.4 Handbook of Hydrology */
  r_air = 0.003486 * pz/(275 + tair);

  /* Set constant values for later use */
  *CONST_1 = slope;
  *CONST_2 = r_air * CP_PM * vpd / ra;
  *CONST_3 = lv * (slope + gamma * (1 + rarc / ra));
  *CONST_4 = lv * gamma / ra;
  *CONST_5 = rs / (lai * Tfactor * vpdfactor) * DAYfactor;
  
}
