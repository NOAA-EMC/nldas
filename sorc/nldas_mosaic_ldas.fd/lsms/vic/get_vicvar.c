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
#include <math.h>
#include "ftn.h"
#include "vicNl.h"
//BOP
//
// !ROUTINE: get_vicvar
//
// !DESCRIPTION: 
//   Extracts specified vic output variable from the output data structure 
//   and prepares it for writing
//
// !INTERFACE:
void FTN(get_vicvar)(int *index,
                     int *nch,
                     float *var){
//EOP
  extern out_data_struct *outdata1; 
  extern atmos_data_struct *atmos; 
  extern dmy_struct dmy; 

  int i; 
  if(*index == 1) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].swnet/outdata1[i].count; 
    }
  } 
  else if(*index == 2) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].lwnet/outdata1[i].count; 
    }
  }
  else if(*index == 3) {
    for(i=0; i<*nch; i++){
      var[i] = (-1.) * outdata1[i].qle/outdata1[i].count; 
    }
  }
  else if(*index == 4) {
    for(i=0; i<*nch; i++){
      var[i] = (-1) * outdata1[i].qh/outdata1[i].count; 
    }
  }
  else if(*index == 5) {
    for(i=0; i<*nch; i++){
      var[i] = (-1) * outdata1[i].qg/outdata1[i].count; 
    }
  }
  else if(*index == 6) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].rainf/outdata1[i].count;
    }
  }
  else if(*index == 7) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].snowf/outdata1[i].count;
    }
  }
  else if(*index == 8) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].evap/outdata1[i].count;
    }
  }
  else if(*index == 9) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].qs/outdata1[i].count;
    }
  }
  else if(*index == 10) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].qsb/outdata1[i].count;
    }
  }
  else if(*index == 11) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].qfz/outdata1[i].count;
    }
  }
  else if(*index == 12) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].snowt+KELVIN;
    }
  }
  else if(*index == 13) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].avgsurft;
    }
  }
  else if(*index == 14) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].radt;
    }
  }
  else if(*index == 15) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].albedo;
    }
  }
  else if(*index == 16) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].soilt[0]+KELVIN;
    }
  }
  else if(*index == 17) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].soilt[1]+KELVIN;
    }
  }
  else if(*index == 18) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].soilt[2]+KELVIN;
    }
  }
  else if(*index == 19) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].moist[0];
    }
  }
  else if(*index == 20) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].moist[1];
    }
  }
  else if(*index == 21) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].moist[2];
    }
  }
  else if(*index == 22) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].tveg/outdata1[i].count;
    }
  }
  else if(*index == 23) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].esoil/outdata1[i].count;
    }
  }
  else if(*index == 24) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].soilwet;
    }
  }
  else if(*index == 25) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].rootmoist;
    }
  }
  else if(*index == 26) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].swe*1000;
    }
  }
  else if(*index == 27) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].qsm/outdata1[i].count;
    }
  }
  else if(*index == 28) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].moist[0]+outdata1[i].moist[1]+
	outdata1[i].moist[2]-outdata1[i].moist_prev; 
    }
  }
  else if(*index == 29) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].swe-outdata1[i].swe_prev; 
    } 
  }
  else if(*index == 30) {
    for(i=0; i<*nch; i++){
      var[i] = outdata1[i].acond; 
    }
  }
  else if(*index == 31) {
    for(i=0; i<*nch; i++){
      var[i] = atmos[i].wind; 
    }
  }
  else if(*index == 32) {
    for(i=0; i<*nch; i++){
      if ( atmos[i].air_temp+KELVIN > 273.15 ) 
      {
         var[i] = atmos[i].prec/(dmy.dt); /* rainfall rate */
      }
      else
      {
         var[i] = 0.0;
      }
    }
  }
  else if(*index == 33) {
    for(i=0; i<*nch; i++){
      if ( atmos[i].air_temp+KELVIN <= 273.15 )
      {
         var[i] = atmos[i].prec/(dmy.dt); /* snowfall rate */
      }
      else
      {
         var[i] = 0.0;
      }
    }
  }
  else if(*index == 34) {
    for(i=0; i<*nch; i++){
      var[i] = atmos[i].air_temp+KELVIN; 
    }
  }
  else if(*index == 35) {
    for(i=0; i<*nch; i++){
      var[i] = .622*atmos[i].vp / (atmos[i].pressure - .378*atmos[i].vp); /* qair */
    }
  }
  else if(*index == 36) {
    for(i=0; i<*nch; i++){
      var[i] = atmos[i].pressure*1000; 
    }
  }
  else if(*index == 37) {
    for(i=0; i<*nch; i++){
      var[i] = atmos[i].shortwave;
    }
  }
  else if(*index == 38) {
    for(i=0; i<*nch; i++){
      var[i] = atmos[i].longwave;
    }
  
  }

}
