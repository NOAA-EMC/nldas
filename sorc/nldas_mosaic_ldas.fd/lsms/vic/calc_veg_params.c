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
#include <math.h>
#include "vicNl.h"


float calc_veg_displacement(float height) {
/**********************************************************************
  calc_veg_displacement		Keith Cherkauer		January 27, 1997

  This subroutine estimates the displacement height of vegetation
  with a given average height based on equations on page 4.12 of the
  Handbook of Hydrology.
**********************************************************************/

  float value;

  value = 0.67 * height;

  return (value);

}

float calc_veg_height(float displacement) {
/**********************************************************************
  calc_veg_height		Keith Cherkauer		March 3, 1997

  This subroutine backs the vegetation height out of the given
  displacement using the reverse procedure from cal_veg_displacement.
**********************************************************************/

  float value;

  value = displacement / 0.67;

  return (value);

}

float calc_veg_roughness(float height) {
/**********************************************************************
  calc_veg_roughness		Keith Cherkauer		January 27, 1997

  This subroutine estimates the roughness height of vegetation
  with a given average height based on equations on page 4.12 of the
  Handbook of Hydrology.
**********************************************************************/

  float value;

  value = 0.123 * height;

  return (value);

}
