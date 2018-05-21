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
 

energy_bal_struct **make_energy_bal(int nveg, int snowband)
/**********************************************************************
	make_energy_bal	Keith Cherkauer		May 26, 1996

  This routine makes an array of frozen soil data structures, one 
  for each vegetation type and bare soil.

**********************************************************************/
{

  int i, j;
  energy_bal_struct **temp;

  temp = (energy_bal_struct**) lis_calloc(nveg, 
                                          sizeof(energy_bal_struct*),
                                          "make_energy_bal");

  /** Initialize all records to unfrozen conditions */
  for(i = 0; i < nveg; i++) {
    temp[i] = (energy_bal_struct*) lis_calloc(snowband, 
                                              sizeof(energy_bal_struct),
                                              "make_energy_bal");
    for(j = 0; j < snowband; j++) {
      temp[i][j].frozen = FALSE;
      //      if(options.QUICK_FLUX) {
      //        if(options.FULL_ENERGY) {
      //          *Nnodes        = 3;
      //        }
      //        else {
      //          *Nnodes        = 1;
      //        }
      //}
    }
  }

  return temp;
}
