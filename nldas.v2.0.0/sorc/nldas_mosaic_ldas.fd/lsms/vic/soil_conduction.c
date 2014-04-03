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
#include "vicNl.h"

float soil_conductivity(float moist, 
			 float Wu, 
			 float soil_density, 
			 float bulk_density,
			 float quartz) {
/**********************************************************************
  Soil thermal conductivity calculated using Johansen's method.

  Reference: Farouki, O.T., "Thermal Properties of Soils" 1986
	Chapter 7: Methods for Calculating the Thermal Conductivity 
		of Soils

  H.B.H. - refers to the handbook of hydrology.

  porosity = n = porosity
  ratio = Sr = fractionaldegree of saturation
  All K values are conductivity in W/mK
  Wu is the fractional volume of unfrozen water

  UNITS: input in m, kg, s

  Returns K in W/m/K

  float moist         total moisture content (mm/mm)
  float Wu            liquid water content (mm/mm)
  float soil_density  soil density (kg m-3)
  float bulk_density  soil bulk density (kg m-3)
  float quartz        soil quartz content (fraction)

**********************************************************************/
  float Ke;
  float Ki = 2.2;	/* thermal conductivity of ice (W/mK) */
  float Kw = 0.57;	/* thermal conductivity of water (W/mK) */
  float Ksat;
  float Ks;		/* thermal conductivity of solid (W/mK)
			   function of quartz content */
  float Kdry;
  float Sr;		/* fractional degree of saturation */
  float K;
  float porosity;

  Kdry = (0.135*bulk_density+64.7)/(soil_density-0.947*bulk_density);

  if(moist>0.) {

    porosity = 1.0 - bulk_density / soil_density;

    Sr = moist/porosity;

    Ks = pow(7.7,quartz) * pow(2.2,1.0-quartz);

    if(Wu==moist) {

      /** Soil unfrozen **/
      Ksat = pow(Ks,1.0-porosity) * pow(Kw,porosity);
      Ke = 0.7 * log10(Sr) + 1.0;

    }
    else {

      /** Soil frozen **/
      Ksat = pow(Ks,1.0-porosity) * pow(Ki,porosity-Wu) * pow(Kw,Wu);
      Ke = Sr;

    }

    K = (Ksat-Kdry)*Ke+Kdry;
    if(K<Kdry) K=Kdry;

  }
  else K=Kdry;

  return (K); 
}


#define organic_fract 0.00

float volumetric_heat_capacity(float soil_fract,
                                float water_fract,
                                float ice_fract) {
/**********************************************************************
  This subroutine calculates the soil volumetric heat capacity based 
  on the fractional volume of its component parts.

  Constant values are volumetric heat capacities in J/m^3/K
	Soil value is for clay or quartz - assumed for all other types

  float soil_fract   fraction of soil volume composed of actual soil (fract)
  float water_fract  fraction of soil volume composed of liquid water (fract)
  float ice_fract    fraction of soil volume composed of ice (fract)

**********************************************************************/

  float Cs;

  Cs = 2.0e6 * (soil_fract - organic_fract);
  Cs += 4.2e6 * water_fract;
  Cs += 1.9e6 * ice_fract;
  Cs += 2.7e6 * organic_fract;
  Cs += 1.3e3 * (1. - (soil_fract + water_fract + ice_fract + organic_fract));

  return (Cs);

}

#undef organic_fract

void set_node_parameters(float   *dz_node,
			 float   *max_moist_node,
			 float   *expt_node,
			 float   *bubble_node,
			 float   *alpha,
			 float   *beta,
			 float   *gamma,
			 float   *depth,
			 float   *max_moist,
			 float   *expt,
			 float   *bubble,
			 float   *quartz,
                         float   layer_node_fract[MAX_LAYERS+1][MAX_NODES],
			 int       Nnodes,
			 int       Nlayers,
			 char      FS_ACTIVE) {
/**********************************************************************
  This subroutine sets the thermal node soil parameters to constant 
  values based on those defined for the current grid cells soil type.
  Thermal node propertiers for the energy balance solution are also 
  set (these constants are used to reduce the solution time required
  within each iteration).

  float   *dz_node          thermal node thicknes (m)
  float   *max_moist_node   maximum moisture content at thermal node (mm/mm)
  float   *expt_node        exponential at thermal node ()
  float   *bubble_node      bubbling pressure at thermal node (cm)
  float   *alpha            first thermal eqn term ()
  float   *beta             second thermal eqn term ()
  float   *gamma            third thermal eqn term ()
  float   *depth            soil moisture layer thickness (m)
  float   *max_moist        soil moisture layer maximum moisture content (mm)
  float   *bubble           soil moisture layer bubbling pressure (cm)
  float    quartz           soil quartz content (fract)
  float ***ufwc_table_node  table of unfrozen water contents ()
  int       Nnodes           number of soil thermal nodes
  int       Nlayers          number of soil moisture layers
  char      FS_ACTIVE        TRUE if frozen soils are active in grid cell

  Modifications:

  02-11-00 Modified to remove node zone averages, node parameters
           are now set based on the parameter value of the layer 
           in which they fall.  Averaging of layer properties 
	   only occurs if the node falls directly on a layer
	   boundary.                                        KAC
  11-May-04 (Port from 4.1.0) Modified to correct differences between
            calculations to determine maximum node moisture and node
            moisture, so that nodes on the boundary between soil
            layers are computed the same way for both.          TJB

**********************************************************************/


  char   PAST_BOTTOM;
  int    nidx, lidx;
  float Lsum; /* cumulative depth of moisture layer */
  float Zsum; /* cumulative depth of thermal node */
  float deltaL[MAX_LAYERS+1];
  char NO_FLUX = FALSE; 

  PAST_BOTTOM = FALSE;
  lidx = 0;
  Lsum = 0.;
  Zsum = 0.;

  /* set node parameters */
  for(nidx=0;nidx<Nnodes;nidx++) {

    if(Zsum == Lsum + depth[lidx] && nidx != 0 && lidx != Nlayers-1) {
      /* node on layer boundary */
      max_moist_node[nidx] = (max_moist[lidx] / depth[lidx] 
			      + max_moist[lidx+1] / depth[lidx+1]) / 1000 / 2.;
      expt_node[nidx]      = (expt[lidx] + expt[lidx+1]) / 2.;
      bubble_node[nidx]    = (bubble[lidx] + bubble[lidx+1]) / 2.;
    }
    else { 
      /* node completely in layer */
      max_moist_node[nidx] = max_moist[lidx] / depth[lidx] / 1000;
      expt_node[nidx]      = expt[lidx];
      bubble_node[nidx]    = bubble[lidx];
    }      
    Zsum += (dz_node[nidx] + dz_node[nidx+1]) / 2.;
    if(Zsum > Lsum + depth[lidx] && !PAST_BOTTOM) {
      Lsum += depth[lidx];
      lidx++;
      if( lidx == Nlayers ) {
	PAST_BOTTOM = TRUE;
	lidx = Nlayers-1;
      }
    }

  }

  /* set constant variables for thermal calculations */
  for(nidx=0;nidx<Nnodes-2;nidx++) {
    alpha[nidx] = ((dz_node[nidx+2] + dz_node[nidx+1]) / 2.0 
		   + (dz_node[nidx+1] + dz_node[nidx]) / 2.0);
    beta[nidx] = ((dz_node[nidx+2] + dz_node[nidx+1]) 
		  * (dz_node[nidx+2] + dz_node[nidx+1])) / 4.0 
      + ((dz_node[nidx+1]+dz_node[nidx]) 
	 * (dz_node[nidx+1]+dz_node[nidx])) / 4.0;
    gamma[nidx] = ((dz_node[nidx+2] + dz_node[nidx+1]) / 2.0 
		   - (dz_node[nidx+1] + dz_node[nidx]) / 2.0);
  }
  if(NO_FLUX) {
    /* no flux bottom boundary activated */
    alpha[Nnodes-2] = ((dz_node[Nnodes-1] + dz_node[Nnodes-1]) / 2.0 
		       + (dz_node[Nnodes-1] + dz_node[Nnodes-2]) / 2.0);
    beta[Nnodes-2] = ((dz_node[Nnodes-1] + dz_node[Nnodes-1]) 
		      * (dz_node[Nnodes-1] + dz_node[Nnodes-1]))/4.0 
      + ((dz_node[Nnodes-1]+dz_node[Nnodes-2])
	 * (dz_node[Nnodes-1]+dz_node[Nnodes-2])) / 4.0;
    gamma[Nnodes-2] = ((dz_node[Nnodes-1] + dz_node[Nnodes-1]) / 2.0 
		       - (dz_node[Nnodes-1] + dz_node[Nnodes-2]) / 2.0);
  }

  /* set fraction of soil thermal node in each soil layer */
  Lsum = 0;
  deltaL[Nlayers] = 0;
  for(lidx=0;lidx<Nlayers;lidx++) {
    deltaL[lidx] = depth[lidx];
    deltaL[Nlayers] -= depth[lidx];
  }
  for(nidx=0;nidx<Nnodes-1;nidx++) 
    deltaL[Nlayers] += (dz_node[nidx] + dz_node[nidx+1]) / 2.;
  for(lidx=0;lidx<=Nlayers;lidx++) {
    Zsum = -dz_node[0] / 2.;
    for(nidx=0;nidx<Nnodes;nidx++) {
      if ( Zsum < Lsum && Zsum + dz_node[nidx] >= Lsum )
      {
	layer_node_fract[lidx][nidx] = 1. 
	  - (float)linear_interp(Lsum, Zsum, Zsum + dz_node[nidx], 0, 1);
	if ( Lsum + deltaL[lidx] < Zsum + dz_node[nidx] )
        {
	  layer_node_fract[lidx][nidx] -= 
	    (float)linear_interp(Lsum + deltaL[lidx], Zsum, Zsum + 
                                 dz_node[nidx], 1, 0);
        }
      }
      else if ( Zsum < Lsum + deltaL[lidx] && 
	        Zsum + dz_node[nidx] >= Lsum + deltaL[lidx] )
      {
	layer_node_fract[lidx][nidx] = 
	  (float)linear_interp(Lsum + deltaL[lidx], Zsum, 
			       Zsum + dz_node[nidx], 0, 1);
      }
      else if ( Zsum >= Lsum && Zsum + dz_node[nidx] <= Lsum + deltaL[lidx] )
      {
	layer_node_fract[lidx][nidx] = 1;
      }
      else 
      {
         layer_node_fract[lidx][nidx] = 0;
      }
      Zsum += dz_node[nidx];
    }
    Lsum += deltaL[lidx];
  }
    


}

#define N_INTS 5

void distribute_node_moisture_properties(float *moist_node,
					 float *ice_node,
					 float *kappa_node,
					 float *Cs_node,
					 float *dz_node,
					 float *T_node,
					 float *max_moist_node,
					 float *expt_node,
					 float *bubble_node,
					 float *moist,
					 float *depth,
					 float *soil_density,
					 float *bulk_density,
					 float *quartz,
					 int     Nnodes,
					 int     Nlayers,
					 char    FS_ACTIVE,
					 int     frozen_soil) {
/*********************************************************************
  This subroutine determines the moisture and ice contents of each 
  soil thermal node based on the current node temperature and layer
  moisture content.  Thermal conductivity and volumetric heat capacity
  are then estimated for each node based on the division of moisture 
  contents..

  float *moist_node      thermal node moisture content (mm/mm)
  float *ice_node        thermal node ice content (mm/mm)
  float *kappa_node      thermal node thermal conductivity (W m-1 K-1)
  float *Cs_node         thermal node heat capacity (J m-3 K-1)
  float *dz_node         thermal node thickness (m)
  float *T_node          thermal node temperature (C)
  float *max_moist_node  thermal node maximum moisture content (mm/mm)
  float *expt_node       thermal node exponential
  float *bubble_node     thermal node bubbling pressure (cm)
  float *moist           soil layer moisture (mm)
  float *depth           soil layer thickness (m)
  float  soil_density    soil density (kg m-3)
  float *bulk_density    soil layer bulk density (kg m-3)
  float  quartz          soil quartz content (fract)
  int     Nnodes          number of soil thermal nodes
  int     Nlayers         number of soil moisture layers

  Modifications:

  02-11-00 Modified to remove node zone averages, node parameters
           are now set based on the parameter value of the layer 
           in which they fall.  Averaging of layer properties 
	   only occurs if the node falls directly on a layer
	   boundary.                                        KAC
  11-May-04 (Port from 4.1.0) Modified to check that node soil
            moisture is less than or equal to maximum node soil
            moisture, otherwise an error is printed to the screen
            and the model exits.                                TJB

*********************************************************************/

  char PAST_BOTTOM;
  int nidx, lidx;
  float Lsum; /* cumulative depth of moisture layer */
  float Zsum; /* cumulative depth of thermal node */
  float soil_fract;

  lidx = 0;
  Lsum = 0.;
  Zsum = 0.;
  PAST_BOTTOM = FALSE;

  /* node estimates */
  for(nidx=0;nidx<Nnodes;nidx++) {
 
    if(Zsum == Lsum + depth[lidx] && nidx != 0 && lidx != Nlayers-1) {
      /* node on layer boundary */
      moist_node[nidx] = (moist[lidx] / depth[lidx] 
			      + moist[lidx+1] / depth[lidx+1]) / 1000 / 2.;
      soil_fract = (bulk_density[lidx] / soil_density[lidx] 
		    + bulk_density[lidx+1] / soil_density[lidx+1]) / 2.;
    }
    else { 
      /* node completely in layer */
      moist_node[nidx] = moist[lidx] / depth[lidx] / 1000;
      soil_fract = (bulk_density[lidx] / soil_density[lidx]);
    }      

    // Check that node moisture does not exceed maximum node moisture
    if ( moist_node[nidx] > max_moist_node[nidx] ) {
      fprintf(stderr, "Node soil moisture, %f, exceeds maximum node soil moisuttre, %f.", moist_node[nidx], max_moist_node[nidx] );
    }

    if(T_node[nidx] < 0 && (FS_ACTIVE && (frozen_soil==1))) {
      /* compute moisture and ice contents */
      ice_node[nidx] 
	= moist_node[nidx] - maximum_unfrozen_water(T_node[nidx],
					     max_moist_node[nidx], 
					     bubble_node[nidx],
					     expt_node[nidx]);
      if(ice_node[nidx]<0) ice_node[nidx]=0;

      /* compute thermal conductivity */
      kappa_node[nidx] 
	= soil_conductivity(moist_node[nidx], moist_node[nidx] 
			    - ice_node[nidx], soil_density[lidx],
			    bulk_density[lidx], quartz[lidx]);

    }
    else {
      /* compute moisture and ice contents */
      ice_node[nidx]   = 0;
      /* compute thermal conductivity */
      kappa_node[nidx] 
	= soil_conductivity(moist_node[nidx], moist_node[nidx], 
			    soil_density[lidx], bulk_density[lidx], 
			    quartz[lidx]);
    }
    /* compute volumetric heat capacity */
    Cs_node[nidx] = volumetric_heat_capacity(bulk_density[lidx] 
					     / soil_density[lidx],
					     moist_node[nidx] - ice_node[nidx],
					     ice_node[nidx]);

    Zsum += (dz_node[nidx] + dz_node[nidx+1]) / 2.;
    if(Zsum > Lsum + depth[lidx] && !PAST_BOTTOM) {
      Lsum += depth[lidx];
      lidx++;
      if( lidx == Nlayers ) {
	PAST_BOTTOM = TRUE;
	lidx = Nlayers-1;
      }
    }
  }

}

#undef N_INTS

void estimate_layer_ice_content(layer_data_struct *layer,
				float            *dz,
				float            *T,
				float            *max_moist_node,
				float            *expt_node,
				float            *bubble_node,
				float            *depth,
				float            *max_moist,
				float            *expt,
				float            *bubble,
				float            *bulk_density,
				float            *soil_density,
				float            *quartz,
                                float layer_node_fract[MAX_LAYERS+1][MAX_NODES],
				int                Nnodes, 
				int                Nlayers,
				char               FS_ACTIVE,
				int                frozen_soil) {
/**************************************************************
  This subroutine estimates the ice content of all soil 
  moisture layers based on the distribution of soil thermal
  node temperatures.

  layer_struct *layer           structure with all soil moisture layer info
  float       *dz              soil thermal node thicknesses (m)
  float       *T               soil thermal node temperatures (C)
  float       *max_moist_node  soil thermal node max moisture content (mm/mm)
  float       *expt_node       soil thermal node exponential ()
  float       *bubble_node     soil thermal node bubbling pressure (cm)
  float       *depth           soil moisture layer thickness (m)
  float       *max_moist       soil layer maximum soil moisture (mm)
  float       *expt            soil layer exponential ()
  float       *bubble          soil layer bubling pressure (cm)
  float       *bulk_density    soil layer bulk density (kg m-3)
  float        soil_density    soil layer soil density (kg m-3)
  float        quartz          soil layer quartz content (fract)
  int           Nnodes          number of soil thermal nodes
  int           Nlayer          number of soil moisture layers

**************************************************************/


  int    nidx;
  int    lidx;
  float tmp_ice;

  for(lidx=0;lidx<Nlayers;lidx++) {
    layer[lidx].T = 0.;
    layer[lidx].ice = 0.;
    for(nidx=0;nidx<Nnodes;nidx++) {
      if(layer_node_fract[lidx][nidx] > 0) {
	layer[lidx].T += T[nidx] * layer_node_fract[lidx][nidx] * dz[nidx];
	if(T[nidx] < 0 && (frozen_soil==1) && FS_ACTIVE) {
	  tmp_ice = layer[lidx].moist 
	- maximum_unfrozen_water(T[nidx], max_moist[lidx], bubble[lidx], 
				 expt[lidx]);
	  if(tmp_ice < 0) tmp_ice = 0;
	}
	else tmp_ice = 0;
	layer[lidx].ice += tmp_ice * layer_node_fract[lidx][nidx] * dz[nidx];
      }
    }
    layer[lidx].T /= depth[lidx];
    layer[lidx].ice /= depth[lidx];
  }

}

void compute_soil_layer_thermal_properties(layer_data_struct *layer,
					   float            *depth,
					   float            *bulk_density,
					   float            *soil_density,
					   float            *quartz,
					   int                Nlayers) {
/********************************************************************
  This subroutine computes the thermal conductivity and volumetric
  heat capacity of each soil layer based on its current moisture
  and ice contents.  Ice is only present if the frozen soil
  algorithm is activated.

  layer_data_struct *layer          structure with all soil layer variables
  float            *depth          soil layer depths (m)
  float            *bulk_density   soil layer bulk density (kg/m^3)
  float             soil_density   soil layer soil density (kg/m^3)
  float             quartz         soil layer quartz content (fract)
  int                Nlayers        number of soil layers

********************************************************************/

  int lidx;
  float moist, ice;

  /* compute layer thermal properties */
  for(lidx=0;lidx<Nlayers;lidx++) {
    moist = layer[lidx].moist / depth[lidx] / 1000;
    ice = layer[lidx].ice / depth[lidx] / 1000;
    layer[lidx].kappa 
      = soil_conductivity(moist, moist - ice, soil_density[lidx], 
			  bulk_density[lidx], quartz[lidx]);
    layer[lidx].Cs 
      = volumetric_heat_capacity(bulk_density[lidx] / soil_density[lidx], 
				 moist - ice, ice);
  }
}

void find_0_degree_fronts(energy_bal_struct *energy,
			  float            *dz,
			  float            *T,
			  int                Nnodes) {
/***********************************************************************
  This subroutine reads through the soil thermal nodes and determines 
  the depths of all thawing and freezing fronts that are present.

  energy_bal_struct *energy  energy balance variable structure
  float            *dz      thermal node thicknesses (m)
  float            *T       thermal node temperatures (C)
  int                Nnodes  number of defined thermal nodes

***********************************************************************/

  int    nidx, fidx; 
  int    Nthaw; /* number of thawing fronts found */
  int    Nfrost; /* number of frost fronts found */
  float tdepth[MAX_FRONTS]; /* thawing frost depths */
  float fdepth[MAX_FRONTS]; /* freezing front depths */
  float Zsum;
  float deltaz;
  
  /* Initialize parameters */
  Zsum = 0;
  for(nidx=0;nidx<Nnodes-1;nidx++)
    Zsum += (dz[nidx] + dz[nidx+1]) / 2.;
  Nthaw = Nfrost = 0;
  for(fidx=0;fidx<MAX_FRONTS;fidx++) {
    fdepth[fidx] = MISSING;
    tdepth[fidx] = MISSING;
  }

  /* find 0 degree fronts */
  for(nidx=Nnodes-2;nidx>=0;nidx--) {
    deltaz = (dz[nidx] + dz[nidx+1]) / 2.;
    if(T[nidx] > 0 && T[nidx+1] <= 0 && Nthaw<MAX_FRONTS) {
      tdepth[Nthaw] = linear_interp(0,T[nidx],T[nidx+1],Zsum-deltaz,Zsum);
      Nthaw++;
    }
    else if(T[nidx] < 0 && T[nidx+1] >= 0 && Nfrost<MAX_FRONTS) {
      fdepth[Nfrost] = linear_interp(0,T[nidx],T[nidx+1],Zsum-deltaz,Zsum);
      Nfrost++;
    }
    Zsum -= deltaz;
  }

  /* store thaw depths */
  for(fidx=0;fidx<MAX_FRONTS;fidx++) energy->tdepth[fidx] = tdepth[fidx];
  /* store frost depths */
  for(fidx=0;fidx<MAX_FRONTS;fidx++) energy->fdepth[fidx] = fdepth[fidx];
  energy->Nthaw = Nthaw;
  energy->Nfrost = Nfrost;

}

float maximum_unfrozen_water(float T,
                              float max_moist,
                              float bubble,
                              float expt) {
/**********************************************************************
  This subroutine computes the maximum amount of unfrozen water that
  can exist at the current temperature.
**********************************************************************/

  float unfrozen;

  unfrozen = max_moist * pow((-Lf * T) / (T
      + 273.16) / (9.18 * bubble / 100.), -(2.0 / (expt - 3.0)));
  if(unfrozen > max_moist) unfrozen = max_moist;
  if(unfrozen < 0) unfrozen = 0;

  return (unfrozen);
}


layer_data_struct find_average_layer(layer_data_struct *wet,
				     layer_data_struct *dry,
				     float             depth,
				     float             mu) {
/*************************************************************
  This subroutine computes the average soil layer moistures
  between the wet and dry fraction for use in computing 
  energy balance parameters.  Other layer variables are copied 
  from the wet fraction structure since they are they same for 
  wet and dry fractions.
**************************************************************/

  layer_data_struct layer;

  layer = *wet;

  layer.ice = ((wet->ice * mu) + (dry->ice * (1. - mu)));
  layer.moist = ((wet->moist * mu) + (dry->moist * (1. - mu)));

  return(layer);

}

