#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include "vicNl_ldas.h"
#include "preproc.h"
/*
  Author: Justin Sheffield <justin@princeton.edu> July 2001

  Generate a dswrf image from Edas and Goes dswrf images by using 
  Edas if Goes is unavailable. 
  
  Modifications:

*/

void process_dswrf(ldas_index_struct ldas_index, float *grib_data_dswrf, float *grib_data_dswrf1, int var, Stats *stats)
{

  int i;
  int index;

  /* loop through computational cells */
  for(i=0;i<ldas_index.ncell_comp;i++) {

    /* get the array index of the cell */
    index = ldas_index.cell_num[i]-1;
#if DEBUG
    fprintf(stderr, "cell[%d]=%d, ", i, index);
    fprintf(stderr, "data: %f %f, ", grib_data_dswrf[index], grib_data_dswrf1[index]);
#endif

    /* use EDAS if GOES is unavailable. compare using <= in case missing data are -999 or -9999 or -99999 */
    if(grib_data_dswrf[index] <= MISSING) {
      grib_data_dswrf[index] = grib_data_dswrf1[index];
      stats->subs[var]++;
      if(grib_data_dswrf[index] <= MISSING) /* both are zero!!! */
        grib_data_dswrf[index] = 0.0;
    }
#if DEBUG
    fprintf(stderr, "blended: %f\n", grib_data_dswrf[index]);
#endif

  }

}


