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
 * SUMMARY:      SnowInterception.c - simulates snow interception and release
 * USAGE:        
 *
 * AUTHOR:       Brian Connelly and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    29-Aug-1996 at 13:42:17
 * LAST-MOD: Tue Feb 29 14:36:36 2000 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculates the interception and subsequent release of
 *               by the forest canopy using an energy balance approach
 * DESCRIP-END.
 * FUNCTIONS:    SnowInterception()
 * COMMENTS:     Modified for use with VIC-NL code by Keith Cherkauer
 *               on 4-9-98   
 */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"


/*****************************************************************************
  Function name: SnowInterception()

  Purpose      : Calculate snow interception and release by the canopy

  Required     :
    int Dt                 - Model timestep (hours)
    float F                - Fractional coverage
    float LAI              - Leaf Area Index
    float MaxInt           - Maximum rainfall interception storage (m)
    float BaseRa           - Aerodynamic resistance (uncorrected for
                             stability) (s/m)
    float AirDens          - Density of air (kg/m3)
    float EactAir          - Actual vapor pressure of air (Pa) 
    float Lv               - Latent heat of vaporization (J/kg3)
    PIXRAD *LocalRad       - Components of radiation balance for current pixel
                             (W/m2) 
    float Press            - Air pressure (Pa)
    float Tair             - Air temperature (C) 
    float Vpd	           - Vapor pressure deficit (Pa) 
    float Wind             - Wind speed (m/s)
    float *RainFall        - Amount of rain (m)
    float *Snowfall        - Amount of snow (m)
    float *IntRain         - Intercepted rain (m) 
    float *IntSnow         - Snow water equivalent of intercepted snow (m)
    float *TempIntStorage  - Temporary storage for snowmelt and rainfall
                             involved in mass release calculations (m)
    float *VaporMassFlux   - Vapor mass flux to/from intercepted snow
                             (m/timestep)
    float *Tcanopy         - Canopy temperature (C)
    float *MeltEnergy      - Energy used in heating and melting of the snow 
                             (W/m2)

  Returns      : none

  Modifies     :
    float *RainFall        - Amount of rain (m)
    float *Snowfall        - Amount of snow (m)
    float *IntRain         - Intercepted rain (m) 
    float *IntSnow         - Snow water equivalent of intercepted snow (m)
    float *TempIntStorage  - Temporary storage for snowmelt and rainfall
                             involved in mass release calculations (m)
    float *VaporMassFlux   - Vapor mass flux to/from intercepted snow
                             (m/timestep)  
    float *Tcanopy         - Canopy temperature (C)

  Comments     : Only the top canopy layer is taken into account for snow
                 interception.  Snow interception by lower canopy is
                 disregarded.  Rain water CAN be intercepted by lower canopy
                 layers (similar to InterceptionStorage()).
                 Of course:  NO vegetation -> NO interception

  Modifications:
  06-98 included maximum structural loading to prevent the model
        from loading the canopy with more snow than it can handle
        structurally.                                             PXS
  09-98 aerodynamic resistances in the canopy when snow has been  
        intercepted is increased by a factor of 10:  include REF
        Journal of Hydrology, 1998                            KAC, GO'D

*****************************************************************************/
void snow_intercept(float Dt, 
		    float F,  
		    float LAI, 
		    float MaxInt, 
		    float Ra, 
		    float AirDens,
		    float EactAir, 
		    float Lv, 
		    float Shortwave,
		    float Longwave, 
		    float Press, 
		    float Tair, 
		    float Vpd, 
		    float Wind,  
		    float *RainFall,
		    float *SnowFall, 
		    float *IntRain, 
		    float *IntSnow,
		    float *TempIntStorage, 
		    float *VaporMassFlux,
		    float *Tcanopy, 
		    float *MeltEnergy, 
		    int month, 
		    int hour)
{
  float AdvectedEnergy;         /* Energy advected by the rain (W/m2) */
  float BlownSnow;              /* Depth of snow blown of the canopy (m) */
  float DeltaSnowInt;           /* Change in the physical swe of snow
				    interceped on the branches. (m) */
  float Drip;                   /* Amount of drip from intercepted snow as a
				    result of snowmelt (m) */
  float ExcessSnowMelt;         /* Snowmelt in excess of the water holding
				    capacity of the tree (m) */
  float EsSnow;                 /* saturated vapor pressure in the snow pack
                                   (Pa)  */
  float InitialSnowInt;         /* Initial intercepted snow (m) */ 
  float InitialWaterInt;        /* Initial intercepted water (snow and rain)
                                    (m) */ 
  float LatentHeat;             /* Latent heat flux (W/m2) */
  float LongOut;                /* Longwave radiation emitted by canopy 
                                   (W/m2) */
  float Ls;                     /* Latent heat of sublimation (J/(kg K) */
  float MassBalanceError;       /* Mass blalnce to make sure no water is
				    being destroyed/created (m) */
  float MaxWaterInt;            /* Water interception capacity (m) */  
  float MaxSnowInt;             /* Snow interception capacity (m) */
  float NetRadiation;
  float PotSnowMelt;            /* Potential snow melt (m) */
  float RainThroughFall;        /* Amount of rain reaching to the ground (m)
				  */ 
  float RefreezeEnergy;         /* Energy available for refreezing or melt */
  float ReleasedMass;           /* Amount of mass release of intercepted snow
				    (m) */ 
  float SensibleHeat;           /* Sensible heat flux (W/m2) */
  float SnowThroughFall;        /* Amount of snow reaching to the ground (m)
				  */ 
  float Tmp;                    /* Temporary variable */

  float Imax1;                  /* maxium water intecept regardless of temp */
  float IntRainFract;           /* Fraction of intercpeted water which is 
				    liquid */
  float IntSnowFract;           /* Fraction of intercepted water which is 
				    solid */
  float Overload;               /* temp variable to calculated structural 
				    overloading */

  /* Convert Units from VIC (mm -> m) */
  *RainFall /= 1000.;
  *SnowFall /= 1000.;
  *IntRain  /= 1000.;
  MaxInt    /= 1000.;

  /* Initialize Drip, H2O balance, and mass release variables. */
  
  InitialWaterInt = *IntSnow + *IntRain;
  
  *IntSnow /= F;
  *IntRain /= F;
  
  InitialSnowInt = *IntSnow;
  
  Drip = 0.0;
  ReleasedMass = 0.0;
  
  /* Determine the maximum snow interception water equivalent.           
     Kobayashi, D., 1986, Snow Accumulation on a Narrow Board,           
     Cold Regions Science and Technology, (13), pp. 239-245.           
     Figure 4. */  
  
  Imax1 = 4.0* LAI_SNOW_MULTIPLIER * LAI;

  if (Tair < -1.0 && Tair > -3.0)
    MaxSnowInt = (Tair*3.0/2.0) + (11.0/2.0);
  else if (Tair > -1.0) 
    MaxSnowInt = 4.0;  
  else
    MaxSnowInt = 1.0;
  
  /* therefore LAI_ratio decreases as temp decreases */
  
  MaxSnowInt *= LAI_SNOW_MULTIPLIER * LAI;
  
  /* Calculate snow interception. */  
  
  DeltaSnowInt = (1-*IntSnow/MaxSnowInt) * *SnowFall; 
  if (DeltaSnowInt + *IntSnow > MaxSnowInt) 
    DeltaSnowInt = MaxSnowInt - *IntSnow;
  if (DeltaSnowInt < 0.0)  
    DeltaSnowInt = 0.0;
  
  /* Reduce the amount of intercepted snow if windy and cold.         
     Ringyo Shikenjo Tokyo, #54, 1952.                                
     Bulletin of the Govt. Forest Exp. Station,                       
     Govt. Forest Exp. Station, Meguro, Tokyo, Japan.                 
     FORSTX 634.9072 R475r #54.                                       
     Page 146, Figure 10.                                               
     
     Reduce the amount of intercepted snow if snowing, windy, and     
     cold (< -3 to -5 C).                                             
     Schmidt and Troendle 1992 western snow conference paper. */  
  
  if (Tair < -3.0 && DeltaSnowInt > 0.0 && Wind > 1.0) {
    BlownSnow = (0.2 * Wind - 0.2) * DeltaSnowInt;
    if (BlownSnow >= DeltaSnowInt) 
      BlownSnow = DeltaSnowInt;
    DeltaSnowInt -= BlownSnow;
  }
  
  /* now update snowfall and total accumulated intercepted snow amounts */

  if (*IntSnow +  DeltaSnowInt > Imax1) DeltaSnowInt =0.0; 
  
  /* pixel depth    */ 
  SnowThroughFall = (*SnowFall - DeltaSnowInt) * F + (*SnowFall) * (1 - F);
  
  /* physical depth */
  *IntSnow += DeltaSnowInt;
  
  /* Calculate amount of rain intercepted on branches and stored in
     intercepted snow. */  
  
  /* physical depth */
  MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;
  
  if ((*IntRain + *RainFall) <= MaxWaterInt) {
    /* physical depth */
    *IntRain += *RainFall;
    /* pixel depth */
    //printf("h1 f=%f rain=%f\n",F,*RainFall);
    RainThroughFall = *RainFall * (1 - F);      
  }
  else {
    /* pixel depth */
    RainThroughFall = (*IntRain + *RainFall - MaxWaterInt) * F + 
      (*RainFall * (1 - F));
    /* physical depth */
    //printf("h2 through=%f\n",RainThroughFall);
    *IntRain = MaxWaterInt;
  }

  /* at this point we have calculated the amount of snowfall intercepted and
     the amount of rainfall intercepted.  These values have been 
     appropriately subtracted from SnowFall and RainFall to determine 
     SnowThroughfall and RainThroughfall.  However, we can end up with the 
     condition that the total intercepted rain plus intercepted snow is 
     greater than the maximum bearing capacity of the tree regardless of air 
     temp (Imax1).  The following routine will adjust *IntRain and *IntSnow 
     by triggering mass release due to overloading.  Of course since *IntRain
     and *IntSnow are mixed, we need to slough them of as fixed fractions  */

  if (*IntRain + *IntSnow > Imax1) { /*then trigger structural unloading*/
    Overload = (*IntSnow + *IntRain) - Imax1;
    IntRainFract= *IntRain/(*IntRain + *IntSnow);
    IntSnowFract = *IntSnow/(*IntRain + *IntSnow);
    *IntRain = *IntRain - Overload*IntRainFract;
    *IntSnow = *IntSnow - Overload*IntSnowFract;
    RainThroughFall = RainThroughFall + (Overload*IntRainFract)*F;
    //printf("h3 %f \n",RainThroughFall);
    SnowThroughFall = SnowThroughFall + (Overload*IntSnowFract)*F;
  }
  
  /* The canopy temperature is assumed to be equal to the air temperature if 
     the air temperature is below 0C, otherwise the canopy temperature is 
     equal to 0C */
  
  if (Tair > 0.)
    *Tcanopy = 0.;
  else
    *Tcanopy = Tair;

  /* Calculate the net radiation at the canopy surface, using the canopy 
     temperature.  The outgoing longwave is subtracted twice, because the 
     canopy radiates in two directions */

  Tmp = *Tcanopy + 273.15;
  LongOut = STEFAN_B * (Tmp * Tmp * Tmp * Tmp);
  NetRadiation = (1.-NEW_SNOW_ALB)*Shortwave + Longwave - 2 * F * LongOut;
  NetRadiation /= F;

  /* Calculate the vapor mass flux between the canopy and the surrounding 
     air mass */
  
  EsSnow = svp(*Tcanopy); 

  /** Added division by 10 to incorporate change in canopy resistance due
      to smoothing by intercepted snow **/
  *VaporMassFlux = AirDens * (0.622/Press) * (EactAir - EsSnow) / Ra / 10.; 
  *VaporMassFlux /= RHO_W; 

  if (Vpd == 0.0 && *VaporMassFlux < 0.0)
    *VaporMassFlux = 0.0;
  
  /* Calculate the latent heat flux */

  Ls = (677. - 0.07 * *Tcanopy) * 4.1868 * 1000.;
  LatentHeat = Ls * *VaporMassFlux * RHO_W;

  /* Calculate the sensible heat flux */

  SensibleHeat = AirDens * Cp * (Tair - *Tcanopy)/Ra;
  
  /* Calculate the advected energy */

#if TIMESTEPSECS
  AdvectedEnergy = (4186.8 * Tair * *RainFall)/(Dt);
#else
  AdvectedEnergy = (4186.8 * Tair * *RainFall)/(Dt * SECPHOUR);
#endif
      
  /* Calculate the amount of energy available for refreezing */

  RefreezeEnergy = SensibleHeat + LatentHeat + NetRadiation + AdvectedEnergy;

#if TIMESTEPSECS
  RefreezeEnergy *= Dt;
#else
  RefreezeEnergy *= Dt * SECPHOUR;
#endif

  /* if RefreezeEnergy is positive it means energy is available to melt the
     intercepted snow in the canopy.  If it is negative, it means that 
     intercepted water will be refrozen */
  
  /* Update maximum water interception storage */
  
  MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;

  /* Convert the vapor mass flux from a flux to a depth per interval */
#if TIMESTEPSECS
  *VaporMassFlux *= Dt;
#else
  *VaporMassFlux *= Dt * SECPHOUR;
#endif
  
  if (RefreezeEnergy > 0.0) {

    if (-(*VaporMassFlux) > *IntRain) {
      *VaporMassFlux = -(*IntRain);
      *IntRain = 0.;
    }
    else
      *IntRain += *VaporMassFlux;

    PotSnowMelt = min((RefreezeEnergy/Lf/RHO_W), *IntSnow);

#if TIMESTEPSECS
    *MeltEnergy -= (Lf * PotSnowMelt * RHO_W) / (Dt);
#else
    *MeltEnergy -= (Lf * PotSnowMelt * RHO_W) / (Dt *SECPHOUR);
#endif

    if ((*IntRain + PotSnowMelt) <= MaxWaterInt) {

      *IntSnow -= PotSnowMelt;
      *IntRain += PotSnowMelt;
      PotSnowMelt = 0.0;

    }
    
    else {

      ExcessSnowMelt = PotSnowMelt + *IntRain - MaxWaterInt;
      
      *IntSnow -= MaxWaterInt - (*IntRain);
      *IntRain = MaxWaterInt;
      if (*IntSnow < 0.0) 
        *IntSnow = 0.0;
      
      if (SnowThroughFall > 0.0 && 
	  InitialSnowInt <= MIN_INTERCEPTION_STORAGE) {
        /* Water in excess of MaxWaterInt has been generated.  If it is 
           snowing and there was little intercepted snow at the beginning 
	   of the time step ( <= MIN_INTERCEPTION_STORAGE), then allow the 
	   snow to melt as it is intercepted */
        Drip += ExcessSnowMelt; 
        *IntSnow -= ExcessSnowMelt;
        if (*IntSnow < 0.0) 
          *IntSnow = 0.0;
      }
      else 
      /* Else, SnowThroughFall = 0.0 or SnowThroughFall > 0.0 and there is a 
         substantial amount of intercepted snow at the beginning of the time 
         step ( > MIN_INTERCEPTION_STORAGE).  Snow melt may generate mass 
         release. */
        *TempIntStorage += ExcessSnowMelt;
      
      MassRelease(IntSnow, TempIntStorage, &ReleasedMass, &Drip);
      //printf("drip2 = %f \n",Drip);
    }
    
    /* If intercepted snow has melted, add the water it held to drip */
    
    MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;
    if (*IntRain > MaxWaterInt) {
      Drip += *IntRain - MaxWaterInt;
      //printf("drip3 = %f \n",Drip);
      *IntRain = MaxWaterInt;
    }
  }  
  
  else /* else (RefreezeEnergy <= 0.0) */ {
    
    /* Reset *TempIntStorage to 0.0 when energy balance is negative */
    
    *TempIntStorage = 0.0;
    
    /* Refreeze as much surface water as you can */
    
    if (RefreezeEnergy > - (*IntRain) * Lf) {
      *IntSnow += fabs(RefreezeEnergy) / Lf;
      *IntRain -= fabs(RefreezeEnergy) / Lf;

#if TIMESTEPSECS
      *MeltEnergy += (fabs(RefreezeEnergy) * RHO_W) / (Dt);
#else
      *MeltEnergy += (fabs(RefreezeEnergy) * RHO_W) / (Dt *SECPHOUR);
#endif

      RefreezeEnergy = 0.0;
    }

    else {
      
      /* All of the water in the surface layer has been frozen. */
      
      *IntSnow += *IntRain;
     
      /* Added on April 8 as a test */
      /*       RefreezeEnergy += *IntRain*Lf; */
      /*       *VaporMassFlux = MAX(*VaporMassFlux,  */
      /*                            RefreezeEnergy/(Ls * RHO_W)); */
      
      /* Energy released by freezing of intercepted water is added to the 
         MeltEnergy */

#if TIMESTEPSECS
      *MeltEnergy += (Lf * *IntRain * RHO_W) / (Dt);
#else
      *MeltEnergy += (Lf * *IntRain * RHO_W) / (Dt *SECPHOUR);
#endif
      *IntRain = 0.0;
      
    } 
    
    if (-(*VaporMassFlux) > *IntSnow) {
      *VaporMassFlux = -(*IntSnow);
      *IntSnow = 0.0;
    }
    else
      *IntSnow += *VaporMassFlux;
  } 
  
  *IntSnow *= F;
  *IntRain *= F;
  *MeltEnergy *= F;
  *VaporMassFlux *= F;
  Drip           *= F;
  ReleasedMass   *= F;
  
  /* Calculate intercepted H2O balance. */  
  
  MassBalanceError = (InitialWaterInt - (*IntSnow + *IntRain))
                   + (*SnowFall + *RainFall) - (SnowThroughFall
                   + RainThroughFall + Drip + ReleasedMass)
                   + *VaporMassFlux;
  //printf("rain comps %f %f \n",RainThroughFall, Drip);
  //printf("snow comps %f %f %f\n",SnowThroughFall, ReleasedMass, *SnowFall);
  *RainFall = RainThroughFall + Drip;
  *SnowFall = SnowThroughFall + ReleasedMass;
  //printf("after %f \n", *SnowFall);
  /* Convert Units to VIC (m -> mm) */
  *VaporMassFlux *= -1.;
  *RainFall *= 1000.;
  *SnowFall *= 1000.;
  *IntRain  *= 1000.;

}
