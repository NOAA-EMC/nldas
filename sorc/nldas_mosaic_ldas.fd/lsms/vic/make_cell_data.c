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
 

cell_data_struct **make_cell_data(int veg_type_num, int snowband)
/**********************************************************************
	make_cell_data	Keith Cherkauer		July 9, 1997

  This subroutine makes an array of type cell, which contains soil
  column variables for a single grid cell.

**********************************************************************/
{
  int i;
  cell_data_struct **temp;

  temp = (cell_data_struct**) lis_calloc(veg_type_num,
                                         sizeof(cell_data_struct*),
                                         "make_cell_data");
  for(i=0;i<veg_type_num;i++) {
    temp[i] = (cell_data_struct*) lis_calloc(snowband,
                                             sizeof(cell_data_struct),
                                             "make_cell_data");
/*
  for(j=0;j<snowband;j++)
  { 
     temp[i][j].layer  = (layer_data_struct*) lis_calloc(3, 
                                               sizeof(layer_data_struct),
                                               "make_cell_data"); 
      
  }
*/ 
  }
  return temp;
}
