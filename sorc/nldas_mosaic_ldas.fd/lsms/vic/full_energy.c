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
#include <math.h>
#include "vicNl.h"

void full_energy(int rec,
                 atmos_data_struct   *atmos, 
                 soil_con_struct     *soil_con,
                 veg_con_struct      *veg_con,
                 dist_prcp_struct    *prcp,
                 dmy_struct          dmy,
                 int                 Nlayer, 
                 int                 Nnode,
                 int                 snowband,
                 int                 full_energy_flag,
                 int                 frozen_soil_flag,
                 int                 grnd_flux_flag,
                 int                 quick_flux_flag)
/**********************************************************************
        full_energy     Keith Cherkauer	        January 8, 1997

  This subroutine controls the model core, it solves both the energy
  and water balance models, as well as frozen soils.  

  modifications:
  07-98 restructured to fix problems with distributed precipitation, 
        and to add the ability to solve the snow model at different 
        elevation bands within a single grid cell.                 KAC
  01-19-00 modified to work with the new atmosphere data structure 
           implemented when the radiation forcing routines were 
           updated.  Also modified to use the new simplified 
           soil moisture storage for the frozen soil algorithm.    KAC

**********************************************************************/
{
  /* extern par_struct par; */
  extern veg_lib_struct *veg_lib;
  char                   overstory;
  int                     j;
  int                    Ndist;
  int                    iveg;
  int                    Nveg;
  int                    veg_class;
  int                    band;
  int                    Nbands;
  float                 out_prec[2*MAX_BANDS];
  float                 dp;
  float                 ice0;
  float                 moist;
  float                 surf_atten;
  float                 wind_h;
  float                 height;
  float                 displacement;
  float                 roughness;
  float                 ref_height;
  //  float                 Cv;
  float                 Le;
  float                 Ls;
  float                 Melt[2*MAX_BANDS];
  float                 bare_albedo;
  float                 snow_inflow[MAX_BANDS];
  float                 tmp_wind[3];
  float                 gauge_correction[2];
  veg_var_struct      ***veg_var;
  cell_data_struct    ***cell;
  energy_bal_struct    **energy;
  snow_data_struct     **snow;
  veg_var_struct        *wet_veg_var;
  veg_var_struct        *dry_veg_var;
  veg_var_struct         empty_veg_var;
  /* set local pointers */
  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;

  /* set variables for distributed precipitation */
  //  if(options.DIST_PRCP) 
  //    Ndist = 2;
  //  else 
    Ndist = 1;
  Nbands = snowband;

  /* Set number of vegetation types */
  //Nveg      = veg_con[0].vegetat_type_num;
  Nveg    = 1; 

  /** Set Damping Depth **/
  dp        = soil_con->dp;

  /**    if(par.rank==3) {
      printf("SOIL PARAM..%d %f\n",rec,soil_con->dp);
      }
  if(par.rank==3){
      printf("Tfactor %d %f\n",rec,soil_con->Pfactor[0]);
      }*/
  /** No correction for now.. */
  gauge_correction[0] = 1; 
  gauge_correction[1] = 1;   
  atmos->out_prec = 0;
  
  /**************************************************
    Solve Energy and/or Water Balance for Each
    Vegetation Type
  **************************************************/
  for ( iveg = 0; iveg <= Nveg; iveg++)
  {
    /** Solve Veg Type only if Coverage Greater than 0% **/
    //    if ((iveg <  Nveg && veg_con[iveg].Cv  > 0.) || 
    //	(iveg == Nveg && veg_con[0].Cv_sum < 1.)) 
    if ( iveg < Nveg )
    {
      //printf("iveg %d %d \n",iveg,Nveg);
      //if ( iveg < Nveg ) Cv = 1; 
      //else Cv = 0;

      /**************************************************
        Initialize Model Parameters
      **************************************************/

      /* Initialize energy balance variables */
      for ( band = 0; band < Nbands; band++ )
      {
         //if(soil_con->AreaFract[band] > 0) {
         energy[iveg][band].shortwave = 0;
         energy[iveg][band].longwave  = 0.;
         //}
      }

      /* Initialize snow variables */
      for ( band = 0; band < Nbands; band++ )
      {
         //if(soil_con->AreaFract[band] > 0) {
         snow[iveg][band].vapor_flux        = 0.;
         snow[iveg][band].canopy_vapor_flux = 0.;
         snow_inflow[band]                  = 0.;
         Melt[band*2]                       = 0.;
         //}
      }

      /* Initialize precipitation storage */
      for ( j = 0; j < 2*MAX_BANDS; j++ )
      {
         out_prec[j] = 0;
      }
    
      /** Define vegetation class number **/
      if ( iveg < Nveg ) 
      {
         veg_class = veg_con[iveg].veg_class;
      }
      else 
      {
         veg_class = 0;
      }
      if ( iveg < Nveg ) 
      {
         wind_h = veg_lib[veg_class].wind_h;
      }
      else
      {
        /** HARD CODED from gp->wind_h*/
        wind_h = 10.0;   
      }
      /** Compute Surface Attenuation due to Vegetation Coverage **/
      if ( iveg < Nveg )
      {
         surf_atten = exp(-veg_lib[veg_class].rad_atten 
                          * veg_lib[veg_class].LAI[dmy.month-1]);
      }
      else 
      {
         surf_atten = 1.;
      }
        
      if ( (full_energy_flag == 1) || (frozen_soil_flag == 1) )
      {
         prepare_full_energy(iveg,Nveg,Nnode,prcp,soil_con,
                             &moist,&ice0,Nlayer,snowband,frozen_soil_flag);
      }

      /** Compute Bare Soil (free of snow) Albedo **/
      if ( iveg != Nveg ) 
      {
         bare_albedo = veg_lib[veg_class].albedo[dmy.month-1];
      }
      else 
      {
         bare_albedo = BARE_SOIL_ALBEDO;
      }

      /*************************************
       Compute the aerodynamic resistance 
      *************************************/

      /* Set surface descriptive variables */
      overstory = FALSE;
      if ( iveg < Nveg )
      {
        displacement = veg_lib[veg_class].displacement[dmy.month-1];
        roughness = veg_lib[veg_class].roughness[dmy.month-1];
        overstory = veg_lib[veg_class].overstory;
      }
      if ( iveg == Nveg || roughness == 0 )
      {
        displacement = 0.;
        roughness = soil_con->rough;
        overstory = FALSE;
      }

      /* Initialize wind speeds */
      tmp_wind[0] = atmos->wind;
      tmp_wind[1] = MISSING;
      tmp_wind[2] = MISSING;
 
      /* Estimate vegetation height */
      height = calc_veg_height(displacement);

      /* Estimate reference height */
      if ( displacement < wind_h ) 
      {
         ref_height = wind_h;
      }
      else 
      {
         ref_height = displacement + wind_h + roughness;
      }

      /* Compute aerodynamic resistance over various surface types */
      CalcAerodynamic(overstory, iveg, Nveg, veg_lib[veg_class].wind_atten,
                      height, soil_con->rough, soil_con->snow_rough,
                      &displacement, &roughness, &ref_height,
                      veg_lib[veg_class].trunk_ratio,
                      tmp_wind, cell[WET][iveg][0].aero_resist);

      /******************************
        Solve ground surface fluxes 
      ******************************/

      for ( band = 0; band < Nbands; band++ ) 
      {
         //if( soil_con->AreaFract[band] > 0 ) {

         if ( iveg < Nveg )
         {
            wet_veg_var = &(veg_var[WET][iveg][band]);
            dry_veg_var = &(veg_var[DRY][iveg][band]);
         }
         else
         {
            wet_veg_var = &(empty_veg_var);
            dry_veg_var = &(empty_veg_var);
         }

         surface_fluxes(rec,overstory, band, veg_class, iveg, Nveg, Ndist, 
                        Nbands, Nlayer, Nnode,dp, ice0, moist, 
                        surf_atten, height, displacement, roughness, 
                        ref_height, bare_albedo, 
                        cell[WET][iveg][0].aero_resist, 
                        &(cell[WET][iveg][band].baseflow), 
                        &(cell[DRY][iveg][band].baseflow), 
                        &(cell[WET][iveg][band].runoff), 
                        &(cell[DRY][iveg][band].runoff), &out_prec[band*2], 
                        tmp_wind, &Le, &Ls, &(Melt[band*2]), 
                        &(cell[WET][iveg][band].inflow), 
                        &(cell[DRY][iveg][band].inflow), 
                        &snow_inflow[band], gauge_correction,
                        veg_con[iveg].root, atmos, soil_con, dmy, 
                        &(energy[iveg][band]), &(snow[iveg][band]),
                        cell[WET][iveg][band].layer, 
                        cell[DRY][iveg][band].layer, 
                        wet_veg_var, dry_veg_var,full_energy_flag,
                        frozen_soil_flag,grnd_flux_flag,quick_flux_flag);
         atmos->out_prec += out_prec[band*2];
         //} /** End Loop Through Elevation Bands **/
      } /** End Full Energy Balance Model **/
 
    } /** end current vegetation type **/
  } /** end of vegetation loop **/
}

