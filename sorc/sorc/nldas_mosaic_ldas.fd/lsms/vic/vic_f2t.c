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
#include <math.h>
#include <vicNl.h>
#include "ftn.h"
//BOP
//
// !ROUTINE: vic_f2t
//
// !DESCRIPTION: 
//  Transfers forcing to VIC model tiles.
//
// !REVISION HISTORY: 
// 14 Oct 2003; James Geiger : Initial Specification
//
// !INTERFACE:
void FTN(vic_f2t)(int *t, 
		  float *forcing,
		  int *ts)
//EOP
{
  float uwind, vwind; 
  extern atmos_data_struct *atmos;
  extern dmy_struct dmy; 
  //BOC
#if 0
  dmy.dt = *ts/3600; /* in hours */
#else
  dmy.dt = (float)*ts; /* in secs */
#endif
  atmos[*t-1].air_temp = (float)forcing[0]-KELVIN;
#if 0
  atmos[*t-1].prec = (float)forcing[7]*dmy.dt*3600;
#else
  atmos[*t-1].prec = (float)forcing[7]*dmy.dt;
#endif
  atmos[*t-1].shortwave = (float)forcing[2];
  atmos[*t-1].vp = (float)((forcing[1]*forcing[6])/
			   (0.622+0.378*forcing[1]))/1000;
  atmos[*t-1].vpd = svp(atmos[*t-1].air_temp)-atmos[*t-1].vp;
  atmos[*t-1].longwave = (float)(forcing[3]);
  uwind  = (float)(forcing[4])*(forcing[4]);
  vwind  = (float)(forcing[5])*(forcing[5]);
  atmos[*t-1].wind = sqrt(uwind+vwind);
  atmos[*t-1].density = ((float)forcing[6]-
        0.378*atmos[*t-1].vp*1000)/((float)forcing[0]*287.0423113);
  atmos[*t-1].pressure = (float)forcing[6]/1000;
  //EOC
}
