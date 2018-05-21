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
#include "vicNl.h"
#include "ftn.h"
void make_dist_prcp(int *snowband)
{
  extern ctile_spmd tspmd; 
  extern par_struct par; 
  /* extern veg_con_struct **veg_con;  */
  extern dist_prcp_struct *prcp; 
  
  int i,nveg, j; 
  
  for ( i = 0; i < tspmd.cdi_array[par.rank]; i++ )
  {
    //nveg = veg_con[0][0].vegetat_type_num; 
    nveg = 1; 
    //prcp[i].mu = (float *)lis_calloc(nveg+1, sizeof(float),"make_dist_prcp"); 
  
    //   for(j=0; j<nveg+1; j++){
    //  prcp[i].mu[j]=1; 
    //} 
    prcp[i].snow = make_snow_data(nveg+1, *snowband); 
    prcp[i].energy = make_energy_bal(nveg+1,*snowband); 
    
    for ( j = 0; j < 2; j++ )
    {
      prcp[i].veg_var[j] = make_veg_var(nveg, *snowband); 
      prcp[i].cell[j] = make_cell_data(nveg+1, *snowband); 
    }
  }
}
