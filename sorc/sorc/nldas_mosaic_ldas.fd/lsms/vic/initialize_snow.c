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

void initialize_snow (snow_data_struct **snow, 
		      int                veg_num,
		      int                cellnum,
		      int                snowband)
/**********************************************************************
	initialize_snow		Keith Cherkauer		January 22, 1997

  This routine initializes the snow variable arrays for each new
  grid cell.

  VARIABLES INITIALIZED:
    snow[i][j].snow;	          TRUE = snow, FALSE = no snow 
    snow[i][j].last_snow;         time steps since last snowfall 
    snow[i][j].snow_canopy;       amount of snow on canopy (m) 
    snow[i][j].swq;               snow water equivalent of the entire pack (m) 
    snow[i][j].surf_water;        liquid water content of the surface 
                                  layer (m) 
    snow[i][j].pack_water;        liquid water content of the snow pack (m) 
    snow[i][j].surf_temp;         depth averaged temperature of the snow pack
                                  surface layer (C) 
    snow[i][j].pack_temp;         depth averaged temperature of the snow pack
                                  (C) 
    snow[i][j].vapor_flux;        depth of water evaporation, sublimation, or 
                                  condensation from snow pack (m) 
    snow[i][j].canopy_vapor_flux; depth of water evaporation, sublimation, or 
                                  condensation from intercepted snow (m) 
    snow[i][j].albedo;            snow surface albedo (fraction) 
    snow[i][j].coldcontent;       cold content of snow pack 
    snow[i][j].mass_error;        snow mass balance error 
    snow[i][j].density;	          snow density (kg/m^3) 
    snow[i][j].depth;	          snow depth (m) 
    snow[i][j].tmp_int_storage;   temporary canopy storage, used in 
                                  snow_canopy 
    snow[i][j].Qnet;              Net energy error in snow model 
    snow[i][j].band_elev;         median elevation of the current snow band 
    snow[i][j].prec_frac;         fracton of precipitation that falls in the 
	  		          current snow band 

  modifications:
  07-09-98 modified to initialize snow variables for each defined
           snow elevation band.                                   KAC
  01-11-99 modified to read new initial snow conditions file format KAC
  04-17-00 removed call for read_initial_snow properties file, the
           file read is now incorporated into a single model state
           file.                                                  KAC

**********************************************************************/
{

  int i, j;
  //  int startlayer;

  //  startlayer=0;

  for ( i = 0 ; i <= veg_num ; i++ ) {
    for ( j = 0 ; j < snowband ; j++ ) {
      snow[i][j].snow            = 0;
      snow[i][j].last_snow       = 0;
      snow[i][j].swq             = 0.0;
      snow[i][j].surf_temp       = 0.0;
      snow[i][j].density         = 0.0;
      snow[i][j].coverage        = 0.0;
      snow[i][j].pack_water      = 0.0;
      snow[i][j].surf_water      = 0.0;
      snow[i][j].vapor_flux      = 0.0;
      snow[i][j].pack_temp       = 0.0;
      snow[i][j].snow_canopy     = 0.0;
      snow[i][j].tmp_int_storage = 0.0;
      if(snow[i][j].density>0.) 
	snow[i][j].depth = 1000. * snow[i][j].swq / snow[i][j].density;
      else snow[i][j].depth = 0.;
      snow[i][j].melt = 0.0; 
    }
  }
}
