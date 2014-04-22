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
 * SUMMARY:      SnowMelt.c - Calculate snow accumulation and melt
 * USAGE:        
 *
 * AUTHOR:       Mark Wigmosta and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 08:50:06
 * LAST-MOD: Tue Mar  7 17:02:40 2000 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate snow accumulation and melt using an energy balance
 *               approach for a two layer snow model
 * DESCRIP-END.
 * FUNCTIONS:    SnowMelt()
 * COMMENTS:     
 */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"


/*****************************************************************************
  Function name: SnowMelt()

  Purpose      : Calculate snow accumulation and melt using an energy balance
                 approach for a two layer snow model

  Required     :
    float delta_t               - Model timestep (hours)
    float z2           - Reference height (m) 
    float displacement          - Displacement height (m)
    float aero_resist  - Aerodynamic resistance (uncorrected for
                                   stability) (s/m)
    float atmos->density        - Density of air (kg/m3)
    float atmos->vp             - Actual vapor pressure of air (Pa) 
    float Le           - Latent heat of vaporization (J/kg3)
    float atmos->net_short      - Net exchange of shortwave radiation (W/m2)
    float atmos->longwave       - Incoming long wave radiation (W/m2)
    float atmos->pressure       - Air pressure (Pa)
    float RainFall              - Amount of rain (m)
    float Snowfall              - Amount of snow (m)
    float atmos->air_temp       - Air temperature (C)
    float atmos->vpd            - Vapor pressure deficit (Pa)
    float wind                  - Wind speed (m/s)
    float snow->pack_water      - Liquid water content of snow pack 
    float snow->surf_water	 - Liquid water content of surface layer 
    float snow->swq             - Snow water equivalent at current pixel (m)
    float snow->vapor_flux;     - Mass flux of water vapor to or from the
                                   intercepted snow (m/time step)
    float snow->pack_temp       - Temperature of snow pack (C)
    float snow->surf_temp       - Temperature of snow pack surface layer (C)
    float snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)

  Modifies     :
    float atmos->melt           - Amount of snowpack outflow (m)
    float snow->pack_water      - Liquid water content of snow pack 
    float snow->surf_water	 - Liquid water content of surface layer 
    float snow->swq             - Snow water equivalent at current pixel (m)
    float snow->vapor_flux;     - Mass flux of water vapor to or from the
                                   intercepted snow (m)
    float snow->pack_temp       - Temperature of snow pack (C)
    float snow->surf_temp       - Temperature of snow pack surface layer (C)
    float snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)

  Comments     :
*****************************************************************************/
void snow_melt(soil_con_struct  *soil_con, 
	       int               iveg,
               float            z2,
               float            aero_resist,
               float            Le,
               snow_data_struct *snow, 
               float            delta_t,
               float            displacement,
               float            Z0,
               float            surf_atten,
               float            rainfall,
               float            snowfall,
               float            wind,
               float            grnd_surf_temp,
	       float            air_temp,
	       float            net_short,
	       float            longwave,
	       float            density,
	       float            pressure,
	       float            vpd,
	       float            vp,
	       float           *melt,
               float           *save_advection,
               float           *save_deltaCC,
               float           *save_grnd_flux,
               float           *save_latent,
               float           *save_sensible,
               float           *save_Qnet,
	       float           *save_refreeze_energy)
{

  float DeltaPackCC;            /* Change in cold content of the pack */
  float DeltaPackSwq;           /* Change in snow water equivalent of the
                                   pack (m) */
  float Ice;                    /* Ice content of snow pack (m)*/
  float InitialSwq;             /* Initial snow water equivalent (m) */
  float MassBalanceError;       /* Mass balance error (m) */
  float MaxLiquidWater;         /* Maximum liquid water content of pack (m) */
  float OldTSurf;               /* Old snow surface temperature (C) */
  float PackCC;                 /* Cold content of snow pack (J) */
  float PackSwq;                /* Snow pack snow water equivalent (m) */
  float Qnet;                   /* Net energy exchange at the surface (W/m2) */
  float RefreezeEnergy;         /* refreeze energy (W/m2) */
  float RefrozenWater;          /* Amount of refrozen water (m) */
  float SnowFallCC;             /* Cold content of new snowfall (J) */
  float SnowMelt;               /* Amount of snow melt during time interval
                                   (m water equivalent) */
  float SurfaceCC;              /* Cold content of snow pack (J) */
  float SurfaceSwq;             /* Surface layer snow water equivalent (m) */
  float SnowFall;
  float RainFall;
  float vapor_flux;
  float advection;
  float deltaCC;
  float grnd_flux;		/* thermal flux through snowpack from ground */
  float latent_heat;		/* thermal flux through snowpack from ground */
  float sensible_heat;		/* thermal flux through snowpack from ground */
  float melt_energy = 0.;

  SnowFall = snowfall / 1000.; /* convet to m */
  RainFall = rainfall / 1000.; /* convet to m */

  InitialSwq = snow->swq;
  OldTSurf = snow->surf_temp;

  /* Initialize snowpack variables */
  
  Ice  = snow->swq - snow->pack_water - snow->surf_water;
  
  /* Reconstruct snow pack */
  if (Ice > MAX_SURFACE_SWE)
    SurfaceSwq = MAX_SURFACE_SWE;
  else
    SurfaceSwq = Ice;
  PackSwq = Ice - SurfaceSwq;
  
  /* Calculate cold contents */
  SurfaceCC = CH_ICE * SurfaceSwq * snow->surf_temp;
  PackCC = CH_ICE * PackSwq * snow->pack_temp;
  if (air_temp > 0.0)
    SnowFallCC = 0.0;
  else
    SnowFallCC = CH_ICE * SnowFall * air_temp;
  
  /* Distribute fresh snowfall */
  if (SnowFall > (MAX_SURFACE_SWE - SurfaceSwq)) {
    DeltaPackSwq = SurfaceSwq + SnowFall - MAX_SURFACE_SWE;
    if (DeltaPackSwq > SurfaceSwq)
      DeltaPackCC = SurfaceCC + (SnowFall - MAX_SURFACE_SWE)/SnowFall *
        SnowFallCC;
    else
      DeltaPackCC = DeltaPackSwq/SurfaceSwq * SurfaceCC;
    SurfaceSwq = MAX_SURFACE_SWE;
    SurfaceCC += SnowFallCC - DeltaPackCC;
    PackSwq += DeltaPackSwq;
    PackCC += DeltaPackCC;
  }
  else {
    SurfaceSwq += SnowFall;
    SurfaceCC += SnowFallCC;
  }
  if (SurfaceSwq > 0.0)
    snow->surf_temp = SurfaceCC/(CH_ICE * SurfaceSwq);
  else 
    snow->surf_temp = 0.0;
  if (PackSwq > 0.0)    
    snow->pack_temp = PackCC/(CH_ICE * PackSwq);
  else
    snow->pack_temp = 0.0;

  /* Adjust ice and snow->surf_water */
  Ice += SnowFall;
  snow->surf_water += RainFall;
  
  /* Calculate the surface energy balance for snow_temp = 0.0 */
  
  vapor_flux = snow->vapor_flux;

  Qnet = CalcSnowPackEnergyBalance((float)0.0, delta_t, aero_resist,
				   z2, displacement, Z0, wind, net_short, 
				   longwave, density, Le, air_temp,
				   pressure * 1000., vpd * 1000., vp * 1000.,
				   RainFall, SurfaceSwq, snow->surf_water, 
				   OldTSurf, &RefreezeEnergy, &vapor_flux, 
				   &advection, &deltaCC, grnd_surf_temp,
				   snow->depth, snow->density, surf_atten, 
				   &grnd_flux, &latent_heat, &sensible_heat);

  snow->vapor_flux = vapor_flux;
  save_refreeze_energy[0] = RefreezeEnergy;

  /* If Qnet == 0.0, then set the surface temperature to 0.0 */
  if (Qnet == 0.0) {
    snow->surf_temp = 0.0;
    if (RefreezeEnergy >= 0.0) {
      RefrozenWater = RefreezeEnergy/(Lf * RHO_W) 
#if TIMESTEPSECS
          * delta_t;
#else
          * delta_t * SECPHOUR; 
#endif
      if (RefrozenWater > snow->surf_water) {
        RefrozenWater = snow->surf_water;
        RefreezeEnergy = RefrozenWater * Lf * RHO_W/
#if TIMESTEPSECS
          (delta_t);
#else
          (delta_t * SECPHOUR);
#endif
      } 
      melt_energy  += RefreezeEnergy;
      SurfaceSwq   += RefrozenWater;
      Ice          += RefrozenWater;
      snow->surf_water   -= RefrozenWater;
      assert(snow->surf_water >= 0.0);
      SnowMelt      = 0.0;
    }
    else {
      
      /* Calculate snow melt */      
      SnowMelt = fabs(RefreezeEnergy)/(Lf * RHO_W) * 
#if TIMESTEPSECS 
        delta_t;
#else
        delta_t * SECPHOUR;
#endif
      melt_energy += RefreezeEnergy;
    }

    /* Convert vapor mass flux to a depth per timestep and adjust snow->surf_water */
#if TIMESTEPSECS
    snow->vapor_flux *= delta_t;
#else
    snow->vapor_flux *= delta_t * SECPHOUR;
#endif
    
    if (snow->surf_water < -(snow->vapor_flux)) {
      snow->vapor_flux = -(snow->surf_water);
      snow->surf_water    = 0.0;
    }
    else
      snow->surf_water += snow->vapor_flux;
    
    /* If SnowMelt < Ice, there was incomplete melting of the pack */
    
    if (SnowMelt < Ice) {
      if (SnowMelt <= PackSwq) {
        snow->surf_water += SnowMelt;
        PackSwq    -= SnowMelt;
        Ice        -= SnowMelt;
      }
      else {
        snow->surf_water += SnowMelt + snow->pack_water;
        snow->pack_water  = 0.0;
        PackSwq     = 0.0;
        Ice        -= SnowMelt;
        SurfaceSwq  = Ice;
      }
    }
    
    /* Else, SnowMelt > Ice and there was complete melting of the pack */
    else {
      SnowMelt    = Ice;
      snow->surf_water += Ice;
      SurfaceSwq  = 0.0;
      snow->surf_temp      = 0.0;
      PackSwq     = 0.0;
      snow->pack_temp      = 0.0;
      Ice         = 0.0;
    }
  }
  
  /* Else, SnowPackEnergyBalance(T=0.0) <= 0.0 */
  else  {
    /* Calculate surface layer temperature using "Brent method" */

    vapor_flux = snow->vapor_flux;

    snow->surf_temp = root_brent((float)(snow->surf_temp-SNOW_DT), 
				 (float)0.0,
				 SnowPackEnergyBalance, delta_t, 
				 aero_resist, z2, 
				 displacement, Z0, wind, net_short, longwave,
				 density, Le, air_temp, pressure * 1000.,
				 vpd * 1000., vp * 1000., RainFall, 
				 SurfaceSwq,
				 snow->surf_water, OldTSurf, &RefreezeEnergy, 
				 &vapor_flux, &advection, &deltaCC, 
				 grnd_surf_temp, snow->depth,
				 snow->density,surf_atten,&grnd_flux,
				 &latent_heat, &sensible_heat);

    if(snow->surf_temp <= MISSING + 1 )
      ErrorSnowPackEnergyBalance(snow->surf_temp, delta_t, aero_resist,
				 z2, displacement, Z0, wind, net_short,
				 longwave, density, Le, air_temp,
				 pressure * 1000., vpd * 1000., vp * 1000.,
				 RainFall, SurfaceSwq, snow->surf_water, 
				 OldTSurf,
				 &RefreezeEnergy, &vapor_flux, &advection, 
				 &deltaCC, grnd_surf_temp,
				 snow->depth, snow->density,surf_atten,
				 &grnd_flux,&latent_heat,
				 &sensible_heat);

    Qnet = CalcSnowPackEnergyBalance(snow->surf_temp, delta_t, aero_resist,
				     z2, displacement, Z0, wind, net_short,
				     longwave, density, Le, air_temp,
				     pressure * 1000., vpd * 1000., 
				     vp * 1000.,
				     RainFall, SurfaceSwq, snow->surf_water, 
				     OldTSurf,
				     &RefreezeEnergy, &vapor_flux, 
				     &advection, &deltaCC, grnd_surf_temp,
				     snow->depth, snow->density,surf_atten,
				     &grnd_flux,&latent_heat,
				     &sensible_heat);

    snow->vapor_flux = vapor_flux;
    save_refreeze_energy[0] = RefreezeEnergy;

    /* since we iterated, the surface layer is below freezing and no snowmelt
     */ 
    
    SnowMelt = 0.0;
    
    /* Since updated snow_temp < 0.0, all of the liquid water in the surface
       layer has been frozen */ 
    
    SurfaceSwq += snow->surf_water;
    Ice        += snow->surf_water;
    snow->surf_water  = 0.0;
    melt_energy += snow->surf_water * Lf 
#if TIMESTEPSECS
                       * RHO_W/(delta_t);
#else 
                       * RHO_W/(delta_t * SECPHOUR);
#endif
    
    /* Convert mass flux to a depth per timestep and adjust SurfaceSwq */
    
#if TIMESTEPSECS
    snow->vapor_flux *= delta_t;
#else
    snow->vapor_flux *= delta_t * SECPHOUR;
#endif

    if (SurfaceSwq < -(snow->vapor_flux)) {
      snow->vapor_flux = -SurfaceSwq;
      SurfaceSwq    = 0.0;
      Ice           = PackSwq;
    }
    else {
      SurfaceSwq += snow->vapor_flux;
      Ice += snow->vapor_flux;
    }

  }
  
  /* Done with iteration etc, now Update the liquid water content of the
     surface layer */ 
  
  MaxLiquidWater = LIQUID_WATER_CAPACITY * SurfaceSwq;
  if  (snow->surf_water > MaxLiquidWater) {
    melt[0]    = snow->surf_water - MaxLiquidWater;
    snow->surf_water = MaxLiquidWater;
  }
  else
    melt[0] = 0.0;
  
  /* Refreeze liquid water in the pack.                                   
     variable 'RefreezeEnergy' is the heat released to the snow pack
     if all liquid water were refrozen.                                   
     if RefreezeEnergy < PackCC then all water IS refrozen           
     PackCC always <=0.0 

     WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this does 
     not involve energy transported to the pixel.  Instead heat from the snow 
     pack is used to refreeze water */
  
  snow->pack_water += melt[0]; /* add surface layer outflow to pack 
                                      liquid water*/
  RefreezeEnergy = snow->pack_water * Lf * RHO_W;

  /* calculate energy released to freeze*/
  
  if (PackCC < -RefreezeEnergy) { /* cold content not fully depleted*/
    PackSwq   += snow->pack_water;    /* refreeze all water and update*/
    Ice       += snow->pack_water;
    snow->pack_water = 0.0;
    if (PackSwq > 0.0) {
      PackCC = PackSwq * CH_ICE * snow->pack_temp + RefreezeEnergy;
      snow->pack_temp = PackCC / (CH_ICE * PackSwq);
      if(snow->pack_temp > 0.) snow->pack_temp = 0.;
    }
    else 
      snow->pack_temp = 0.0;
  }
  else { 
    /* cold content has been either exactly satisfied or exceeded. If
       PackCC = refreeze then pack is ripe and all pack water is
       refrozen, else if energy released in refreezing exceeds PackCC 
       then exactly the right amount of water is refrozen to satify PackCC.
       The refrozen water is added to PackSwq and Ice */

    snow->pack_temp      = 0.0;    
    DeltaPackSwq = -PackCC/(Lf * RHO_W); 
    snow->pack_water -= DeltaPackSwq;
    PackSwq += DeltaPackSwq;
    Ice += DeltaPackSwq;
  }
  
  /* Update the liquid water content of the pack */
  
  MaxLiquidWater = LIQUID_WATER_CAPACITY * PackSwq;
  if (snow->pack_water > MaxLiquidWater) {
    melt[0]    = snow->pack_water - MaxLiquidWater;
    snow->pack_water = MaxLiquidWater;
  }
  else
    melt[0] = 0.0;

  /* Update snow properties */
  
  Ice  = PackSwq + SurfaceSwq;

  if (Ice > MAX_SURFACE_SWE) {
    SurfaceCC   = CH_ICE * snow->surf_temp * SurfaceSwq;
    PackCC      = CH_ICE * snow->pack_temp * PackSwq;
    if (SurfaceSwq > MAX_SURFACE_SWE) {
      PackCC     += SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq; 
      SurfaceCC  -= SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq; 
      PackSwq    += SurfaceSwq - MAX_SURFACE_SWE;
      SurfaceSwq -= SurfaceSwq - MAX_SURFACE_SWE;
    }
    else if ( SurfaceSwq < MAX_SURFACE_SWE) {
      PackCC     -= PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / PackSwq; 
      SurfaceCC  += PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / PackSwq; 
      PackSwq    -= MAX_SURFACE_SWE - SurfaceSwq;
      SurfaceSwq += MAX_SURFACE_SWE - SurfaceSwq;
    }
    snow->pack_temp      = PackCC / (CH_ICE * PackSwq);
    snow->surf_temp      = SurfaceCC / (CH_ICE * SurfaceSwq);
  }
  else {
    PackSwq = 0.0;
    PackCC  = 0.0;
    snow->pack_temp  = 0.0;
  }
  snow->swq = Ice + snow->pack_water + snow->surf_water;

  if (snow->swq == 0.0) {
    snow->surf_temp = 0.0;
    snow->pack_temp = 0.0;
  }
    
  /* Mass balance test */
  
  MassBalanceError = (InitialSwq - snow->swq) + (RainFall + SnowFall) 
    - melt[0] + snow->vapor_flux; 
  
/*  printf("%d %d %g\n", y, x, MassBalanceError);*/

  melt[0] *= 1000.; /* converts back to mm */
  //snow->mass_error = MassBalanceError;
  snow->coldcontent = SurfaceCC;
  snow->vapor_flux *= -1.;
  *save_advection = advection;
  *save_deltaCC = deltaCC;
  *save_grnd_flux = grnd_flux;
  *save_latent = latent_heat;
  *save_sensible = sensible_heat;
  *save_Qnet = Qnet;

}

/*****************************************************************************
  Function name: CalcSnowPackEnergyBalance()

  Purpose      : Dummy function to make a direct call to
                 SnowEnergyBalance() possible.

  Required     : 
    float TSurf - SnowPack surface temperature (C)
    other arguments required by SnowPackEnergyBalance()

  Returns      :
    float Qnet - Net energy exchange at the SnowPack snow surface (W/m^2)

  Modifies     : none

  Comments     : function is local to this module
*****************************************************************************/
float CalcSnowPackEnergyBalance(float Tsurf, ...)
{

  va_list ap;                   /* Used in traversing variable argument list
                                 */ 
  float Qnet;                   /* Net energy exchange at the SnowPack snow
                                   surface (W/m^2) */

  va_start(ap, Tsurf);
  Qnet = SnowPackEnergyBalance(Tsurf, ap);
  va_end(ap);
  
  return Qnet;
}

float ErrorSnowPackEnergyBalance(float Tsurf, ...)
{

  va_list ap;                   /* Used in traversing variable argument list
                                 */ 
  float Qnet;                   /* Net energy exchange at the SnowPack snow
                                   surface (W/m^2) */

  va_start(ap, Tsurf);
  Qnet = ErrorPrintSnowPackEnergyBalance(Tsurf, ap);
  va_end(ap);
  
  return Qnet;
}

float ErrorPrintSnowPackEnergyBalance(float TSurf, va_list ap)
{


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

  float *LatentHeat;		/* Latent heat exchange at surface (W/m2) */
  float *SensibleHeat;		/* Sensible heat exchange at surface (W/m2) */

  char lismsg[MAXSTRING];
  int  CELL_NUM = -1;

  /* initialize variables */
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
  
  /* print variables */
  sprintf(lismsg,"DBG: snow_melt -- Dt = %f [ %d ]\n",Dt,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Ra = %f [ %d ]\n",Ra,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Z = %f [ %d ]\n",Z,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Displacement = %f [ %d ]\n",Displacement,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Z0 = %f [ %d ]\n",Z0,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Wind = %f [ %d ]\n",Wind,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- ShortRad = %f [ %d ]\n",ShortRad,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- LongRadIn = %f [ %d ]\n",LongRadIn,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- AirDens = %f [ %d ]\n",AirDens,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Lv = %f [ %d ]\n",Lv,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Tair = %f [ %d ]\n",Tair,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Press = %f [ %d ]\n",Press,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Vpd = %f [ %d ]\n",Vpd,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- EactAir = %f [ %d ]\n",EactAir,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Rain = %f [ %d ]\n",Rain,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- SweSurfaceLayer = %f [ %d ]\n",SweSurfaceLayer,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- SurfaceLiquidWater = %f [ %d ]\n",SurfaceLiquidWater,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- OldTSurf = %f [ %d ]\n",OldTSurf,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- RefreezeEnergy = %f [ %d ]\n",RefreezeEnergy[0],CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- VaporMassFlux = %f [ %d ]\n",VaporMassFlux[0],CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- AdvectedEnergy = %f [ %d ]\n",AdvectedEnergy[0],CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- DeltaColdContent = %f [ %d ]\n",DeltaColdContent[0],CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- TGrnd = %f [ %d ]\n",TGrnd,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- SnowDepth = %f [ %d ]\n",SnowDepth,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- SnowDensity = %f [ %d ]\n",SnowDensity,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- SurfAttenuation = %f [ %d ]\n",SurfAttenuation,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- GroundFlux = %f [ %d ]\n",GroundFlux[0],CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- LatentHeat = %f [ %d ]\n",LatentHeat[0],CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- SensibleHeat = %f [ %d ]\n",SensibleHeat[0],CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: snow_melt -- Finished dumping snow_melt variables.\nTry increasing SNOW_DT to get model to complete cell.\nThen check output for instabilities. [ %d ]\n",CELL_NUM);
  lis_log_msgC(lismsg);
  //  vicerror("Finished dumping snow_melt variables.\nTry increasing SNOW_DT to get model to complete cell.\nThen check output for instabilities.");

  return(0.0);

}

