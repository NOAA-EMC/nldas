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

void initialize_veg(veg_var_struct      **veg_var,
		    int                  vegetat_type_num,
		    int                  snowband)
/**********************************************************************
  initialize_veg		Dag Lohmann	 January 1996

  This routine initailizes the vegetation variable array.

  Modifications:
  07-13-98 modified to initialize vegetation structure for all 
           defined elevation bands                                 KAC

**********************************************************************/
{

  int i, j;

  for ( i = 0 ; i < vegetat_type_num ; i++) {
    for ( j = 0 ; j < snowband; j++ ) {
      veg_var[i][j].Wdew = 0.0;
      veg_var[i][j].throughfall = 0.0;
    }
  }
}
