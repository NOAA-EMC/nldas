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

/* #include "vicNl_def.h" */
/** vic_module.h **/
soil_con_struct *soil_con; 
veg_con_struct **veg_con; 
dist_prcp_struct *prcp; 
atmos_data_struct *atmos;
dmy_struct dmy; 
veg_lib_struct *veg_lib; 
out_data_struct *outdata1; 
model_state_struct *model_state; 
par_struct par; 
ctile_spmd tspmd; 
/* float *outdata_var;  */
