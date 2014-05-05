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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ftn.h"
#include "vicNl.h"
#define MAXSTR      ( 1024 )
#define MAX_SOILS   (   12 )
#define MAX_FIELDS  (   22 )

enum soil_fields {
  SOIL_ID=0,
  SOIL_DESC,
  SOIL_SAND,
  SOIL_SILT,
  SOIL_SOIL_DENSITY,
  SOIL_BULK_DENSITY,
  SOIL_BUBBLE,
  SOIL_QUARTZ,
  SOIL_POROSITY,
  SOIL_RESID_MOIST,
  SOIL_WCR,
  SOIL_WPWP,
  SOIL_KSAT,
  SOIL_B,
  SOIL_EXPT,
  SOIL_DS,
  SOIL_DSMAX,
  SOIL_WS,
  SOIL_BINFILT,
  SOIL_C,
  SOIL_PHIS,
  SOIL_ROUGH
};

/* local functions */

/****
 * set_soilparam
 * reads in soil type parameters from a look-up-table and 
 * assigns values to cells based on the soil type at that cell
****/
void  FTN(set_soilparam)(char *filename, int *len, int *tex, 
			 int *nch, int *nlayer) {
  
  extern FILE  *open_file(char string[], char type[]);
  extern soil_con_struct *soil_con; 
  char *realfile;
  FILE *soilparam;
  int i, j; 
  char buf[MAXSTR];
  float data[MAX_SOILS][MAX_FIELDS];

  realfile = (char *) lis_malloc(*len+1,"set_soilparam");
  i = 0;
  while(i<*len){
    realfile[i] = filename[i]; 
    i++;
  }
  realfile[*len] = '\0';
    
  printf("FILENAME>>> %s \n", realfile);
  soilparam = open_file(realfile,"r");
  printf("Successfully opened the file..\n");
  free(realfile);

  /* skip the title line */
  fgets(buf, MAXSTR, soilparam); 
  //puts(buf);

  /* read data into temporary array */
  for(i=0; i<MAX_SOILS; i++) {
    fscanf(soilparam, "%*d %*s");
    for(j=2; j<MAX_FIELDS; j++) {
      fscanf(soilparam, "%f", &data[i][j]);
    }
  }

  /* test get_soil_type function */

  /* assign values to each cell */
  for(i=0; i<*nch; i++) {

    soil_con[i].b_infilt = data[tex[i]][SOIL_BINFILT];
    soil_con[i].Ds = data[tex[i]][SOIL_DS];
    soil_con[i].Dsmax = data[tex[i]][SOIL_DSMAX];
    soil_con[i].Ws = data[tex[i]][SOIL_WS];
    soil_con[i].c = data[tex[i]][SOIL_C];
    soil_con[i].rough = data[tex[i]][SOIL_ROUGH]; 

    for(j=0; j<*nlayer; j++) {

      soil_con[i].expt[j] = data[tex[i]][SOIL_EXPT];
      soil_con[i].Ksat[j] = data[tex[i]][SOIL_KSAT];
      //soil_con[i].phi_s[j] = data[tex[i]][SOIL_PHIS];
      soil_con[i].bubble[j] = data[tex[i]][SOIL_BUBBLE];
      soil_con[i].quartz[j] = data[tex[i]][SOIL_QUARTZ];
      soil_con[i].bulk_density[j] = data[tex[i]][SOIL_BULK_DENSITY];
      soil_con[i].soil_density[j] = data[tex[i]][SOIL_SOIL_DENSITY]; 
      soil_con[i].porosity[j] = 1.0 - soil_con[i].bulk_density[j]
	/ soil_con[i].soil_density[j]; 
      soil_con[i].resid_moist[j] = data[tex[i]][SOIL_RESID_MOIST];
      soil_con[i].Wcr[j] = data[tex[i]][SOIL_WCR];
      soil_con[i].Wpwp[j] = data[tex[i]][SOIL_WPWP]; 
    }

    /* some default values: these may be overwritten by parameter maps */
//<kluge orig -- justin's verification code>
    soil_con[i].depth[0] = 0.100; 
    soil_con[i].depth[1] = 1.500; 
    soil_con[i].depth[2] = 0.300; 
//    soil_con[i].depth[0] = 0.100; 
//    soil_con[i].depth[1] = 0.300; 
//    soil_con[i].depth[2] = 0.600; 
//</kluge>
    /* initial moisture content set at critical level =
     * 1000 (m to mm) x critical level fraction x porosity x layer depth
     */
    for(j=0; j<*nlayer; j++) {
      soil_con[i].init_moist[j] = 1000.0 * soil_con[i].Wcr[j] * (1.0 - soil_con[i].bulk_density[j]/soil_con[i].soil_density[j]) * soil_con[i].depth[j];
      soil_con[i].max_moist[j] = soil_con[i].depth[j] * soil_con[i].porosity[j] * 1000.;
      if(soil_con[i].init_moist[j] > soil_con[i].max_moist[j]) { 
	soil_con[i].init_moist[j] = soil_con[i].max_moist[j];
      }
    }
//<kluge orig -- justin's verification code>
    soil_con[i].snow_rough = 0.01; 
    //soil_con[i].annual_prec = 834.8; 
    soil_con[i].elevation = 201.00; 
    soil_con[i].avg_temp = -1.3; 
    soil_con[i].dp = 4.00; 
//    soil_con[i].snow_rough = 0.0025; 
//    //soil_con[i].annual_prec = 550.0; 
//    soil_con[i].elevation = 500.00; 
//    soil_con[i].avg_temp = 3.0; 
//    soil_con[i].dp = 3.0; 
//</kluge>

  }

  fclose(soilparam);

    printf("DBG: set_soilparam -- Done. \n");
}

