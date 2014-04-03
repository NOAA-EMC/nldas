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

/*
 * SUMMARY:      SnowPackEnergyBalance.c - Calculate snow pack energy balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 09:09:29
 * LAST-MOD: Thu Apr  6 12:43:30 2000 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate snow pack energy balance
 * DESCRIP-END.
 * FUNCTIONS:    SnowPackEnergyBalance()
 * COMMENTS:     
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

#define GRAMSPKG 1000.
#define CH_WATER 4186.8e3
#define JOULESPCAL     4.1868   /* Joules per calorie */
#define EPS            0.622    /* ratio of molecular weight of water vapor to
                                   that for dry air */
#define CP          1013.0      /* Specific heat of moist air at constant 
                                   pressure (J/(kg*C)) */
/*****************************************************************************
  Function name: SnowPackEnergyBalance()

  Purpose      : Calculate the surface energy balance for the snow pack

  Required     :
    float TSurf           - new estimate of effective surface temperature
    va_list ap            - Argument list initialized by va_start().  For
                            elements of list and order, see beginning of
                            routine

  Returns      :
    float RestTerm        - Rest term in the energy balance

  Modifies     : 
    float *RefreezeEnergy - Refreeze energy (W/m2) 
    float *VaporMassFlux  - Mass flux of water vapor to or from the
                            intercepted snow (m/s)

  Comments     :
    Reference:  Bras, R. A., Hydrology, an introduction to hydrologic
                science, Addisson Wesley, Inc., Reading, etc., 1990.
*****************************************************************************/
float SnowPackEnergyBalance(float TSurf, va_list ap)
{


  //  const char *Routine = "SnowPackEnergyBalance";

  /* start of list of arguments in variable argument list */

  float Dt;                     /* Model time step (hours) */
  float Ra;                     /* Aerodynamic resistance (s/m) */
  float Z;                      /* Reference height (m) */
  float Displacement;           /* Displacement height (m) */
  float Z0;                     /* surface roughness height (m) */
  float Wind;                   /* Wind speed (m/s) */
  float ShortRad;               /* Net incident shortwave radiation (W/m2) */
  float LongRadIn;              /* Incoming longwave radiation (W/m2) */
  float AirDens;                /* Density of air (kg/m3) */
  float Lv;                     /* Latent heat of vaporization (J/kg3) */
  float Tair;                   /* Air temperature (C) */
  float Press;                  /* Air pressure (Pa) */
  float Vpd;			/* Vapor pressure deficit (Pa) */
  float EactAir;                /* Actual vapor pressure of air (Pa) */
  float Rain;                   /* Rain fall (m/timestep) */
  float SweSurfaceLayer;        /* Snow water equivalent in surface layer (m)
                                 */ 
  float SurfaceLiquidWater;     /* Liquid water in the surface layer (m) */
  float OldTSurf;               /* Surface temperature during previous time
                                   step */ 
  float *GroundFlux;		/* Ground Heat Flux (W/m2) */
  float *RefreezeEnergy;        /* Refreeze energy (W/m2) */
  float *VaporMassFlux;          /* Mass flux of water vapor to or from the
                                   intercepted snow */
  float *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
  float *DeltaColdContent;       /* Change in cold content (W/m2) */
  float TGrnd;     
  float SnowDepth;  
  float SnowDensity; 
  float SurfAttenuation; 
  
  /* end of list of arguments in variable argument list */

  float Density;                /* Density of water/ice at TMean (kg/m3) */
  float EsSnow;                 /* saturated vapor pressure in the snow pack
                                   (Pa)  */
  float *LatentHeat;		/* Latent heat exchange at surface (W/m2) */
  float LongRadOut;		/* long wave radiation emitted by surface
				   (W/m2) */
  float Ls;                     /* Latent heat of sublimation (J/kg) */
  float NetRad;			/* Net radiation exchange at surface (W/m2) */
  float RestTerm;		/* Rest term in surface energy balance
				   (W/m2) */
  float *SensibleHeat;		/* Sensible heat exchange at surface (W/m2) */
  float  TMean;                /* Average temperature for time step (C) */

  /* Assign the elements of the array to the appropriate variables.  The list
     is traversed as if the elements are floats, because:

     In the variable-length part of variable-length argument lists, the old
     ``default argument promotions'' apply: arguments of type float are
     always promoted (widened) to type double, and types char and short int
     are promoted to int. Therefore, it is never correct to invoke
     va_arg(argp, float); instead you should always use va_arg(argp,
     double). 

     (quoted from the comp.lang.c FAQ list)
     */
  Dt                 = (float) va_arg(ap, double);
  Ra                 = (float) va_arg(ap, double);
  Z                  = (float) va_arg(ap, double);
  Displacement       = (float) va_arg(ap, double);
  Z0                 = (float) va_arg(ap, double);
  Wind               = (float) va_arg(ap, double);
  ShortRad           = (float) va_arg(ap, double);
  LongRadIn          = (float) va_arg(ap, double);
  AirDens            = (float) va_arg(ap, double);
  Lv                 = (float) va_arg(ap, double);
  Tair               = (float) va_arg(ap, double);
  Press              = (float) va_arg(ap, double);
  Vpd                = (float) va_arg(ap, double);
  EactAir            = (float) va_arg(ap, double);
  Rain               = (float) va_arg(ap, double);
  SweSurfaceLayer    = (float) va_arg(ap, double);
  SurfaceLiquidWater = (float) va_arg(ap, double);
  OldTSurf           = (float) va_arg(ap, double);
  RefreezeEnergy     = (float *) va_arg(ap, double *);
  VaporMassFlux      = (float *) va_arg(ap, double *);
  AdvectedEnergy     = (float *) va_arg(ap, double *);
  DeltaColdContent   = (float *) va_arg(ap, double *);
  TGrnd              = (float) va_arg(ap, double);
  SnowDepth          = (float) va_arg(ap, double);
  SnowDensity        = (float) va_arg(ap, double);
  SurfAttenuation    = (float) va_arg(ap, double);
  GroundFlux         = (float *) va_arg(ap, double *);
  LatentHeat         = (float *) va_arg(ap, double *);
  SensibleHeat       = (float *) va_arg(ap, double *);
  
  /* Calculate active temp for energy balance as average of old and new  */
  
/*   TMean = 0.5 * (OldTSurf + TSurf); */
  TMean = TSurf;
  Density = RHO_W;
  
  /* Correct aerodynamic conductance for stable conditions
     Note: If air temp >> snow temp then aero_cond -> 0 (i.e. very stable)
     velocity (vel_2m) is expected to be in m/sec */                        
  
  /* Apply the stability correction to the aerodynamic resistance 
     NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
     think that it is more correct to calculate ALL fluxes at the same
     reference level */


  if (Wind > 0.0)
    Ra /= StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0);
    /*Ra /= StabilityCorrection(2.f, 0.f, TMean, Tair, Wind, Z0);*/
  else
    Ra = HUGE_RESIST;

  /* Calculate longwave exchange and net radiation */
 
  LongRadOut = STEFAN * (TMean+273.15) * (TMean+273.15) 
    * (TMean+273.15) * (TMean+273.15);
  NetRad = SurfAttenuation * ShortRad + LongRadIn - LongRadOut;
  
  /* Calculate the sensible heat flux */

  *SensibleHeat = AirDens * CP * (Tair - TMean)/Ra;

  /* Calculate the mass flux of ice to or from the surface layer */
 
  /* Calculate the saturated vapor pressure in the snow pack, 
     (Equation 3.32, Bras 1990) */

  EsSnow = svp(TMean) * 1000.;

/*   EsSnow = 610.78 * exp((float)((17.269 * TMean) / (237.3 + TMean))); */

/*   if (TMean < 0.0) */
/*     EsSnow *= 1.0 + .00972 * TMean + .000042 * pow((float)TMean,(float)2.0); */
  
  *VaporMassFlux = AirDens * (EPS/Press) * (EactAir - EsSnow)/Ra;
  *VaporMassFlux /= Density;
  if (Vpd == 0.0 && *VaporMassFlux < 0.0)
    *VaporMassFlux = 0.0;
  
  /* Calculate latent heat flux */
 
  if (TMean >= 0.0) {
    /* Melt conditions: use latent heat of vaporization */
    *LatentHeat = Lv * *VaporMassFlux * Density;
  }
  else {
    /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
    Ls = (677. - 0.07 * TMean) * JOULESPCAL * GRAMSPKG;
    *LatentHeat = Ls * *VaporMassFlux * Density;
  }
  
  /* Calculate advected heat flux from rain 
     WORK IN PROGRESS:  Should the following read (Tair - Tsurf) ?? */
  
#if TIMESTEPSECS 
  *AdvectedEnergy = (CH_WATER * Tair * Rain) / (Dt);
#else
  *AdvectedEnergy = (CH_WATER * Tair * Rain) / (Dt*SECPHOUR);
#endif
  
  /* Calculate change in cold content */
  
#if TIMESTEPSECS 
  *DeltaColdContent = CH_ICE * SweSurfaceLayer * (TSurf - OldTSurf)/
    (Dt);
#else
  *DeltaColdContent = CH_ICE * SweSurfaceLayer * (TSurf - OldTSurf)/
    (Dt * SECPHOUR);
#endif

  /* Calculate Ground Heat Flux */
  if(SnowDepth>0.) {
    *GroundFlux = 2.9302e-6 * SnowDensity * SnowDensity
        * (TGrnd - TMean) / SnowDepth;
  }
  else *GroundFlux=0;
  *DeltaColdContent -= *GroundFlux;
  
  /* Calculate net energy exchange at the snow surface */
  
  RestTerm = NetRad + *SensibleHeat + *LatentHeat + *AdvectedEnergy 
           - *DeltaColdContent ;

#if TIMESTEPSECS 
  *RefreezeEnergy = (SurfaceLiquidWater * Lf * Density)/(Dt);
#else
  *RefreezeEnergy = (SurfaceLiquidWater * Lf * Density)/(Dt * SECPHOUR);
#endif

  if (TSurf == 0.0 && RestTerm > -(*RefreezeEnergy)) {
    *RefreezeEnergy = -RestTerm;  /* available energy input over cold content
                                    used to melt, i.e. Qrf is negative value
                                    (energy out of pack)*/ 
    RestTerm = 0.0;
  }
  else {
    RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
  }
  
  return RestTerm;
}

#undef GRAMSPKG
#undef CH_WATER
#undef JOULESPCAL
#undef EPS
#undef CP
