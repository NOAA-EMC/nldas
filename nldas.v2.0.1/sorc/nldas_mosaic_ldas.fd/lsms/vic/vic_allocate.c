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
#include "misc.h"
#if (defined SPMD)
#include "mpi.h"
#endif
#include "ftn.h"
#include "vicNl.h"
#include "vic_module.h"
//BOP
//
// !ROUTINE: vic_allocate.c
//
// !DESCRIPTION: 
//  Allocates memory for VIC variables
//
// !INTERFACE:
void FTN(vic_allocate)(int *nch, int *cdi_array, int *cdispls, int *snowbands)
  //EOP
{
  int i; 

  veg_con_struct vc_temp; 
  soil_con_struct sc_temp; 
  veg_lib_struct vl_temp; 
  out_data_struct out_temp; 
  model_state_struct ms_temp; 
  /* atmos_data_struct atmos_temp;  */
#if ( defined SPMD )
  MPI_Aint vc_offsets[2];
  int vc_blkcnts[2]           = {7, 1};
  MPI_Datatype vc_oldtypes[2] = {MPI_FLOAT, MPI_INT};

  MPI_Aint sc_offsets[2];
  int sc_blkcnts[2]           = {110, 1};
  MPI_Datatype sc_oldtypes[2] = {MPI_FLOAT, MPI_INT};

  MPI_Aint vl_offsets[3];
  int vl_blkcnts[3]           = {1, 83, 1};
  MPI_Datatype vl_oldtypes[3] = {MPI_CHAR, MPI_FLOAT, MPI_INT};

  MPI_Aint out_offsets[2];
  int out_blkcnts[2]           = {30, 1};
  MPI_Datatype out_oldtypes[2] = {MPI_FLOAT, MPI_INT};

  MPI_Aint atmos_offsets[1];
  int atmos_blkcnts[1]           = {10};
  MPI_Datatype atmos_oldtypes[1] = {MPI_FLOAT};

  MPI_Aint ms_offsets[3];
  int ms_blkcnts[3]           = {4,4,84};
  MPI_Datatype ms_oldtypes[3] = {MPI_CHAR,MPI_INT,MPI_FLOAT};

  MPI_Comm_rank(MPI_COMM_WORLD,&par.rank);
  MPI_Comm_size(MPI_COMM_WORLD,&par.npes);
#else
  par.rank = 0;
  par.npes = 1;
#endif

  tspmd.cdi_array = (int *) lis_calloc(par.npes, sizeof(int),"vic_allocate"); 
  tspmd.cdispls   = (int *) lis_calloc(par.npes, sizeof(int),"vic_allocate");

  for ( i = 0; i < par.npes; i++ )
  {
    tspmd.cdi_array[i] = cdi_array[i];
    tspmd.cdispls[i]   = cdispls[i];
  }
#if (defined SPMD)  
  MPI_Address(&(vc_temp.zone_depth[0]),&(vc_offsets[0]));
  MPI_Address(&(vc_temp.veg_class),&(vc_offsets[1]));
  for(i=1;i>=0; i--)
    vc_offsets[i] = vc_offsets[i] - vc_offsets[0];
  
  MPI_Type_struct(2,vc_blkcnts,vc_offsets,vc_oldtypes,&(par.MPI_VEG_CON)); 
  MPI_Type_commit(&(par.MPI_VEG_CON)); 

  MPI_Address(&(sc_temp.Ds),&(sc_offsets[0]));
  MPI_Address(&(sc_temp.FS_ACTIVE),&(sc_offsets[1]));
  for(i=1;i>=0; i--)
    sc_offsets[i] = sc_offsets[i] - sc_offsets[0];  

  MPI_Type_struct(2,sc_blkcnts,sc_offsets,sc_oldtypes,&(par.MPI_SOIL_CON)); 
  MPI_Type_commit(&(par.MPI_SOIL_CON)); 

  MPI_Address(&(vl_temp.overstory),&(vl_offsets[0]));
  MPI_Address(&(vl_temp.rarc),&(vl_offsets[1]));
  MPI_Address(&(vl_temp.veg_class),&(vl_offsets[2]));
  for(i=2;i>=0; i--)
    vl_offsets[i] = vl_offsets[i] - vl_offsets[0];  

  MPI_Type_struct(3,vl_blkcnts,vl_offsets,vl_oldtypes,&(par.MPI_VEG_LIB));
  MPI_Type_commit(&(par.MPI_VEG_LIB)); 

  MPI_Address(&(out_temp.swnet),&(out_offsets[0]));
  MPI_Address(&(out_temp.count),&(out_offsets[1]));
  for(i=1;i>=0; i--)
    out_offsets[i] = out_offsets[i] - out_offsets[0];  
  
  MPI_Type_struct(2,out_blkcnts,out_offsets,out_oldtypes,&(par.MPI_OUT_STRUCT));
  MPI_Type_commit(&(par.MPI_OUT_STRUCT));

  //MPI_Address(&(atmos_temp.snowflag),&(atmos_offsets[0]));
  atmos_offsets[0] = 0; 

  MPI_Type_struct(1,atmos_blkcnts,atmos_offsets,atmos_oldtypes,
                  &(par.MPI_ATMOS_STRUCT));
  MPI_Type_commit(&(par.MPI_ATMOS_STRUCT));

  MPI_Address(&(ms_temp.MELTING),&(ms_offsets[0]));
  MPI_Address(&(ms_temp.last_snow),&(ms_offsets[1]));
  MPI_Address(&(ms_temp.moist),&(ms_offsets[2]));
  for(i=2;i>=0; i--)
    ms_offsets[i] = ms_offsets[i] - ms_offsets[0];  

  MPI_Type_struct(3,ms_blkcnts,ms_offsets,ms_oldtypes,
                  &(par.MPI_MS_STRUCT));
  MPI_Type_commit(&(par.MPI_MS_STRUCT));
#endif 
  if ( par.rank == 0 )
  {
    soil_con = (soil_con_struct *) lis_calloc(*nch,
                                              sizeof(soil_con_struct),
                                              "vic_allocate");  
    outdata1 = (out_data_struct *) lis_calloc(*nch,
                                              sizeof(out_data_struct),
                                              "vic_allocate"); 
    atmos = (atmos_data_struct *) lis_calloc(*nch,
                                             sizeof(atmos_data_struct),
                                             "vic_allocate");
    model_state = (model_state_struct *) lis_calloc(*nch, 
                                                    sizeof(model_state_struct),
                                                    "vic_allocate"); 
  /* outdata_var = (float *) lis_calloc(*nch, sizeof(float),"vic_allocate"); */
  }
  else
  {
    atmos = (atmos_data_struct *) lis_calloc(tspmd.cdi_array[par.rank],
                                             sizeof(atmos_data_struct),
                                             "vic_allocate");
    soil_con = (soil_con_struct *) lis_calloc(tspmd.cdi_array[par.rank],
                                              sizeof(soil_con_struct),
                                              "vic_allocate");  
    outdata1 = (out_data_struct *) lis_calloc(tspmd.cdi_array[par.rank],
                                              sizeof(out_data_struct),
                                              "vic_allocate"); 
    model_state = (model_state_struct *) lis_calloc(tspmd.cdi_array[par.rank],
                                                    sizeof(model_state_struct),
                                                    "vic_allocate"); 
  }

  veg_con = (veg_con_struct **) lis_calloc(tspmd.cdi_array[par.rank],
                                           sizeof(veg_con_struct *),
                                           "vic_allocate"); 
  prcp = (dist_prcp_struct *) lis_calloc(tspmd.cdi_array[par.rank],
                                         sizeof(dist_prcp_struct),
                                         "vic_allocate"); 
  make_dist_prcp(snowbands);
  
}
