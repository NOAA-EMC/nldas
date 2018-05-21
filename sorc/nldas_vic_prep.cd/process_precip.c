#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include "vicNl_ldas.h"
#include "preproc.h"
/*
  Author: Justin Sheffield <justin@princeton.edu> June 2001

  Generate a precip image from Edas and Combo precip images. Blend Edas
  and Combo according to a weighting mask. Use Edas if Combo is unavailable. 
  
  Modifications:

*/

void process_precip(ldas_index_struct ldas_index, float *grib_data_apcp, float *grib_data_apcp1, int var, Stats *stats)
{

  int i;
  int index;
  float weight;

  /* loop through computational cells */
  for(i=0;i<ldas_index.ncell_comp;i++) {

    /* get the array index of the cell */
    index = ldas_index.cell_num[i]-1;
#if DEBUG
    fprintf(stderr, "cell[%d]=%d, ", i, index);
    fprintf(stderr, "data: %f %f, ", grib_data_apcp[index], grib_data_apcp1[index]);
#endif

    /* use EDAS if Combo is unavailable. compare using <= in case missing data are -999 or -9999 or -99999 */
    if(grib_data_apcp[index] <= MISSING) {
      grib_data_apcp[index] = grib_data_apcp1[index];
      stats->subs[var]++;
      if(grib_data_apcp1[index] <= MISSING) /* edas data are missing as well!!! */
        grib_data_apcp[index] = 0.0;
    }
    /* else blend EDAS and Combo */
    else {
      /* get the weighting from the precip weight mask */
      weight = (float)ldas_index.cell_pweight[i];
      grib_data_apcp[index] = (weight/8.0)*grib_data_apcp[index] + ((8.0-weight)/8.0)*grib_data_apcp1[index];
      stats->blended[var]++;
    }

#if DEBUG
    fprintf(stderr, "blended: %f\n", grib_data_apcp[index]);
#endif

  }

}


