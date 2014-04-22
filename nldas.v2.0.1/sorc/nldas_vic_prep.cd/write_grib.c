#include <stdio.h>
#include <stdlib.h>

#include "preproc.h"

/*  
     
   Write the grib data to file.
   NOTE: see met_names[][] for FILE pointer name association.

   Quality control checks inserted.
					JS 
										

*/

void write_grib( float *grib_data, int *cell_num, int *cell_index, int ncell_comp, 
		 int format, FILE *fp, int var, int rec, FILE *qcp, float min, float max, char *fname, Stats *stats )
{

  int icell;
  int nfails;
  int nwritten;

/*
  printf("write_grib\n");
*/

  /* check the data */
  nfails = qc_data(qcp, grib_data, cell_num, ncell_comp, min, max);
  if(nfails == -1) {
    fprintf(stderr, "QC check: rec %d %s no qc done\n", rec, fname);
    fprintf(qcp, "QC check: rec %d %s no qc done\n", rec, fname);
  }
  else if(nfails == 0) {
/*    fprintf(stderr, "QC check: rec %d %s pass\n", rec, fname);*/
    fprintf(qcp, "QC check: rec %d %s pass\n", rec, fname);
  }
  else {
    fprintf(stderr, "QC check: rec %d %s %d failures\n", rec, fname, nfails);
    fprintf(qcp, "QC check: rec %d %s %d failures\n", rec, fname, nfails);
    stats->nqcfails += nfails;
    stats->qcfails[var] += nfails;
  }

  /* write out data in specified format */
  nwritten = 0;
  if(format==1)
    for(icell=0; icell<ncell_comp; icell++) {
      nwritten += fwrite(&grib_data[cell_num[icell]-1], sizeof(float), 1, fp);
#if DEBUG
      fprintf(stderr, "%d %f\n", cell_num[icell], grib_data[cell_num[icell]-1]);
#endif
    }
  else
    for(icell=0; icell<ncell_comp; icell++) {
      fprintf(fp, "%d %f\n", cell_num[icell], grib_data[cell_num[icell]-1]);
      nwritten++;
    }

  fflush(fp);

  if(nwritten != ncell_comp) {
    fprintf(stderr, "# values written != # of cells in mask\n");
    fprintf(stderr, "nwritten = %d\n", nwritten);
    fprintf(stderr, "ncell = %d\n", ncell_comp);
    exit(0);
  }

}
