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
#include "vicNl.h"
 

snow_data_struct **make_snow_data(int nveg, int snowband)
/**********************************************************************
	make_snow_data	Keith Cherkauer		January 22, 1997

  This routine makes an array of snow cover data structures, one 
  for each vegetation type plus bare soil.

  modifications:
  07-09-98 modified to make te make a two dimensional array which 
           also accounts for a variable number of snow elevation
           bands                                               KAC

**********************************************************************/
{

  int                i;
  snow_data_struct **temp;

  temp = (snow_data_struct **) lis_calloc(nveg, 
				          sizeof(snow_data_struct *),
                                          "make_snow_data");

  for(i=0;i<nveg;i++) {
    temp[i] = (snow_data_struct *) lis_calloc(snowband, 
					  sizeof(snow_data_struct),
                                          "make_snow_data");
  }
    
  return temp;
}
