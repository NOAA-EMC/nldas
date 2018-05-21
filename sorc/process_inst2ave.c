#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include "vicNl_ldas.h"
#include "preproc.h"

/*
  Author: Justin Sheffield <justin@princeton.edu> Oct 2002

  Convert instantaneous data to time step average
  
  Modifications:

*/

void process_inst2ave(ldas_index_struct ldas_index, float *grib_data, float *grib_data_old, int var, Stats *stats)
{

  int j;
  int index;

  /* loop through computational cells */
  for(j=0;j<ldas_index.ncell_comp;j++) {

    /* get the array index of the cell */
    index = ldas_index.cell_num[j]-1;
#if DEBUG
    fprintf(stderr, "cell[%d]=%d, ", j, index);
    fprintf(stderr, "data: %f %f : ", grib_data[index], grib_data_old[index]);
#endif

    /* linearly interpolate between the 2 values */
    if(grib_data[index] > MISSING && grib_data_old[index] > MISSING) {
      grib_data[index] = (grib_data[index] + grib_data_old[index])/2.0;
      stats->ins2ave[var]++;
    }

#if DEBUG
    fprintf(stderr, "%f\n", grib_data[index]);
#endif

  }
  stats->nins2ave++;

}


