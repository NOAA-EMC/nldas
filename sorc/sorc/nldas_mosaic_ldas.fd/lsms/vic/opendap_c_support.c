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
#include "ftn.h"
#include "opendap_c_support.h"
OPENDAP_C_STRUCT opendap_vars;

//BOP
//
// !ROUTINE: setup_c_struct
//
// !DESCRIPTION: 
//  Initializes the OPENDAP structure for GDS-based VIC runs. 
//
// !REVISION HISTORY: 
// 14 Oct 2003; James Geiger :Initial Specification
//
// !INTERFACE:
void FTN(setup_c_struct)(int * iam, int * parm_nc, int * parm_nr, 
			 int * parm_slat, int * parm_nlat, 
			 int * parm_wlon, int * parm_elon, 
			 int * tnroffset)
//EOP
{
  //BOC
   opendap_vars.iam       = *iam;
   opendap_vars.parm_nc   = *parm_nc;
   opendap_vars.parm_nr   = *parm_nr;
   opendap_vars.parm_slat = *parm_slat;
   opendap_vars.parm_nlat = *parm_nlat;
   opendap_vars.parm_wlon = *parm_wlon;
   opendap_vars.parm_elon = *parm_elon;
   opendap_vars.tnroffset = *tnroffset;

   return;
   //EOC
}
