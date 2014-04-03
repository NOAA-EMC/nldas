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
#include "vicNl.h"


void read_snowband(int *nch, int *snowband)
/**********************************************************************
  read_snowband		Keith Cherkauer		July 9, 1998

  This routine reads snow elevation band median elevaton, and 
  precipitation fraction for use with the snow model.

**********************************************************************/
{
  /* extern soil_con_struct *soil_con;  */
  int i; 
  i=0; 
}
#if 0   
  for(i=0; i<*nch; i++){
    
    soil_con[i].Tfactor   = (float *)lis_calloc(*snowband,sizeof(float),"read_snowband");
    soil_con[i].Pfactor   = (float *)lis_calloc(*snowband,sizeof(float),"read_snowband");
    soil_con[i].AreaFract = (float *)lis_calloc(*snowband,sizeof(float),"read_snowband");
    
    //    if (*Tfactor == NULL || *Pfactor == NULL || *AreaFract == NULL) 
    //nrerror("Memory allocation failure in read_snowband");
    
    if(*snowband==1) {
      /** If no snow bands, set factors to use unmodified forcing data **/
      soil_con[i].AreaFract[0] = 1.;
      soil_con[i].Tfactor[0]   = 0.;
      soil_con[i].Pfactor[0]   = 1.;
    }
    
    //    else {
    //      sprintf(ErrStr,"Number of snow elevation bands must be > 0 (%i)",Nbands);
    //      nrerror(ErrStr);
    //    }
  }

  if(Nbands>1) {

    /** Find Current Grid Cell in SnowBand File **/
#if !NO_REWIND
    rewind(snowband);
#endif

    fscanf(snowband, "%i", &cell);
    while(cell != gridcell && !feof(snowband)) {
      fgets(ErrStr,MAXSTRING,snowband);
      fscanf(snowband, "%i", &cell);
    }
    if(feof(snowband)) {
      sprintf(ErrStr,"Cannot find current gridcell (%i) in snow band file\n",
	      gridcell);
      nrerror(ErrStr);
    }

    /** Read Area Fraction **/
    total = 0.;
    for(band = 0; band < Nbands; band++) {
      fscanf(snowband, "%lf", &area_fract);
      if(area_fract<0) {
	sprintf(ErrStr,"Negative snow band area fraction (%f) read from file", 
		area_fract);
	nrerror(ErrStr);
      }
      (*AreaFract)[band]  = area_fract;
      total              += area_fract;
    }
    if(total!=1.) {
      fprintf(stderr,"WARNING: Sum of the snow band area fractions does not equal 1 (%f), dividing each fraction by the sum\n",
	      total);
      for(band = 0; band < options.SNOW_BAND; band++) 
	(*AreaFract)[band] /= total;
    }

    /** Read Band Elevation **/
    for(band = 0; band < Nbands; band++) {
      fscanf(snowband, "%lf", &band_elev);
      if(band_elev<0) {
	fprintf(stderr,"Negative snow band elevation (%f) read from file\n", 
		band_elev);
      }
      (*Tfactor)[band] = (elev - band_elev) / 1000. * T_lapse;
    }
    total = 0.;

    /** Read Precipitation Fraction **/
    for(band = 0; band < options.SNOW_BAND; band++) {
      fscanf(snowband, "%lf", &prec_frac);
      if(prec_frac<0) {
	sprintf(ErrStr,"Snow band precipitation fraction (%f) must be between 0 and 1\n", 
		prec_frac);
	nrerror(ErrStr);
      }
      if(prec_frac>0 && (*AreaFract)[band]==0) {
	sprintf(ErrStr,"Snow band precipitation fraction (%f) should be 0 when the area fraction is 0. (band = %i)\n", 
		prec_frac, band);
	nrerror(ErrStr);
      }
      (*Pfactor)[band] = prec_frac;
      total += prec_frac;
    }
    if(total!=1.) {
      fprintf(stderr,"WARNING: Sum of the snow band precipitation fractions does not equal %i (%f), dividing each fraction by the sum\n",
	      1, total);
      for(band = 0; band < options.SNOW_BAND; band++) 
	(*Pfactor)[band] /= total;
    }
    for (band = 0; band < options.SNOW_BAND; band++) {
      if ((*AreaFract)[band] > 0)
	(*Pfactor)[band] /= (*AreaFract)[band];
      else 
	(*Pfactor)[band]  = 0.;
    }
  }
#endif
