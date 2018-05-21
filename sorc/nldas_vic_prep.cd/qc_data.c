#include <stdio.h>
#include <stdlib.h>

#include "preproc.h"
#include "vicNl_def.h"

/* Justin Sheffield Feb, 2001 
     
   Do quality control on grib data and returns the number of QC failures.
   NOTE: see met_names[][] for FILE pointer name association.

*/

int qc_data(FILE *qcp, float *grib_data, int *comp, int ncomp_cells, float min_bound, float max_bound)
{

  int i, index, nfail;

  nfail = 0;

  /* no QC checks required */
  if(min_bound == MISSING && max_bound == MISSING) {
  
    return (-1);
  
  }
  
  /* Qc checks */
  else {
  
    for(i=0; i<ncomp_cells; i++) {

      index = comp[i]-1;

      if(grib_data[index] < min_bound ||
         grib_data[index] > max_bound) {
       
        nfail++;

#if VERBOSE_QC  
        fprintf(qcp, "%6d %6d %10f %10f %10f\n", i, comp[i], min_bound, grib_data[index], max_bound);
#endif

        /* if it's missing then it's an error, as it should have already been infilled/blended/substituted */ 
        if(grib_data[index] <= MISSING) {
          fprintf(stderr, "%6d %6d %10f %10f %10f\n", i, comp[i], min_bound, grib_data[index], max_bound);
          nrerror("qc_data: missing data found");
        }
        /* else shift any rogue data to within our range */
        else {
          if(grib_data[index] > max_bound) grib_data[index] = max_bound;
          if(grib_data[index] < min_bound) grib_data[index] = min_bound;
        }

      }
    
    }
  }  
  
  return (nfail);
  
}
