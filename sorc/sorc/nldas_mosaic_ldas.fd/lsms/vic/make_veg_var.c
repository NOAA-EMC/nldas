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
 

veg_var_struct **make_veg_var(int veg_type_num, int snowband)
/**********************************************************************
	make_veg_var	Dag Lohman		January 1996

  This routine makes an array of vegitation variables for each vegitation
  type.

  Modifications:
  07-13-98 modified to add structure definitions for all defined 
           elevation bands                                       KAC

**********************************************************************/
{
  
  int              i;
  veg_var_struct **temp;

  temp = (veg_var_struct **) lis_calloc(veg_type_num,
                                        sizeof(veg_var_struct *),
                                        "make_veg_var");

  for(i=0;i<veg_type_num;i++)
  {
     temp[i] = (veg_var_struct *) lis_calloc(snowband,
                                             sizeof(veg_var_struct),
                                             "make_veg_var");
  }

  return temp;
}
