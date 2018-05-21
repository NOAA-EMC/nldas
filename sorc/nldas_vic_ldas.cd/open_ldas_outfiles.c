#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"
#include "vicNl_ldas.h"

/**** Open the LDAS output data files ****/
void open_ldas_outfiles(ldas_outfiles_struct *ldas_outfiles, 
			filenames_struct     *fname,
			dmy_struct           *dmy,
			char                 *mode )
{

  extern option_struct options;

  char   prefix[BUFSIZ+1];
  char   cindex[2];

  int    i;

  /* 1) Energy Balance Components */
  ldas_outfiles->SWnet     = open_file(make_fname(fname->result_dir,"SWnet_",     "",dmy),mode);
  ldas_outfiles->LWnet     = open_file(make_fname(fname->result_dir,"LWnet_",     "",dmy),mode);
  ldas_outfiles->Qle       = open_file(make_fname(fname->result_dir,"Qle_",       "",dmy),mode);
  ldas_outfiles->Qh        = open_file(make_fname(fname->result_dir,"Qh_",        "",dmy),mode);
  ldas_outfiles->Qg        = open_file(make_fname(fname->result_dir,"Qg_",        "",dmy),mode);
  ldas_outfiles->Qf        = open_file(make_fname(fname->result_dir,"Qf_",        "",dmy),mode);
  ldas_outfiles->SWdown    = open_file(make_fname(fname->result_dir,"SWdown_",    "",dmy),mode);
  ldas_outfiles->LWdown    = open_file(make_fname(fname->result_dir,"LWdown_",    "",dmy),mode);
  
  /* 2) Water Balance Components */
  ldas_outfiles->Snowf     = open_file(make_fname(fname->result_dir,"Snowf_",     "",dmy),mode);
  ldas_outfiles->Rainf     = open_file(make_fname(fname->result_dir,"Rainf_",     "",dmy),mode);
  ldas_outfiles->Evap      = open_file(make_fname(fname->result_dir,"Evap_",      "",dmy),mode);
  ldas_outfiles->Qs        = open_file(make_fname(fname->result_dir,"Qs_",        "",dmy),mode);
  ldas_outfiles->Qsb       = open_file(make_fname(fname->result_dir,"Qsb_",       "",dmy),mode);
  ldas_outfiles->Qsm       = open_file(make_fname(fname->result_dir,"Qsm_",       "",dmy),mode);
  ldas_outfiles->Qgwrec    = open_file(make_fname(fname->result_dir,"Qgwrec_",    "",dmy),mode);
  ldas_outfiles->Qrec      = open_file(make_fname(fname->result_dir,"Qrec_",      "",dmy),mode);
  
  /* 3) Surface State Variables */
  ldas_outfiles->SnowT     = open_file(make_fname(fname->result_dir,"SnowT_",     "",dmy),mode);
  ldas_outfiles->VegT      = open_file(make_fname(fname->result_dir,"VegT_",      "",dmy),mode);
  ldas_outfiles->BareT     = open_file(make_fname(fname->result_dir,"BareT_",     "",dmy),mode);
  ldas_outfiles->AvgSurfT  = open_file(make_fname(fname->result_dir,"AvgSurfT_",  "",dmy),mode);
  ldas_outfiles->RadT      = open_file(make_fname(fname->result_dir,"RadT_",      "",dmy),mode);
  ldas_outfiles->Albedo    = open_file(make_fname(fname->result_dir,"Albedo_",    "",dmy),mode);
  ldas_outfiles->SWE       = open_file(make_fname(fname->result_dir,"SWE_",       "",dmy),mode);
  ldas_outfiles->SurfStor  = open_file(make_fname(fname->result_dir,"SurfStor_",  "",dmy),mode);
  ldas_outfiles->CanopInt  = open_file(make_fname(fname->result_dir,"CanopInt_",  "",dmy),mode);
  
  /* 4) Subsurface State Variables */
  ldas_outfiles->SoilMoistTotal = open_file(make_fname(fname->result_dir,"SoilMoistTotal_", "",dmy),mode);
  ldas_outfiles->SoilMoistRoot  = open_file(make_fname(fname->result_dir,"SoilMoistRoot_",  "",dmy),mode);
  ldas_outfiles->SoilMoist1m    = open_file(make_fname(fname->result_dir,"SoilMoist1m_",    "",dmy),mode);
  ldas_outfiles->SoilWetTotal   = open_file(make_fname(fname->result_dir,"SoilWetTotal_",   "",dmy),mode);
  ldas_outfiles->SoilWetRoot    = open_file(make_fname(fname->result_dir,"SoilWetRoot_",    "",dmy),mode);

  /* 5) Evaporation Components */
  ldas_outfiles->ECanop    = open_file(make_fname(fname->result_dir,"ECanop_",    "",dmy),mode);      
  ldas_outfiles->TVeg      = open_file(make_fname(fname->result_dir,"TVeg_",      "",dmy),mode);
  ldas_outfiles->ESoil     = open_file(make_fname(fname->result_dir,"ESoil_",     "",dmy),mode);
  ldas_outfiles->EWater    = open_file(make_fname(fname->result_dir,"EWater_",    "",dmy),mode);
  ldas_outfiles->SubSnow   = open_file(make_fname(fname->result_dir,"SubSnow_",   "",dmy),mode);
  ldas_outfiles->EPot      = open_file(make_fname(fname->result_dir,"EPot_",      "",dmy),mode);
  ldas_outfiles->ACond     = open_file(make_fname(fname->result_dir,"ACond_",     "",dmy),mode);
  ldas_outfiles->CCond     = open_file(make_fname(fname->result_dir,"CCond_",     "",dmy),mode);
  ldas_outfiles->VGreen    = open_file(make_fname(fname->result_dir,"VGreen_",    "",dmy),mode);
  ldas_outfiles->LAI       = open_file(make_fname(fname->result_dir,"LAI_",       "",dmy),mode);
  ldas_outfiles->MRough    = open_file(make_fname(fname->result_dir,"MRough_",    "",dmy),mode);
  ldas_outfiles->HRough    = open_file(make_fname(fname->result_dir,"HRough_",    "",dmy),mode);

  /* 6) Streamflow */
  /* none */
  
  /* 7) Cold Season Processes */
  ldas_outfiles->SnowDepth = open_file(make_fname(fname->result_dir,"SnowDepth_", "",dmy),mode);
  ldas_outfiles->SnowFrac  = open_file(make_fname(fname->result_dir,"SnowFrac_",  "",dmy),mode);
  ldas_outfiles->SAlbedo   = open_file(make_fname(fname->result_dir,"SAlbedo_",   "",dmy),mode);

  /* Non-LDAS outputs, but are useful anyway */
  ldas_outfiles->Precip   = open_file(make_fname(fname->result_dir,"Precip_",   "",dmy),mode);
  ldas_outfiles->RNet     = open_file(make_fname(fname->result_dir,"RNet_",     "",dmy),mode);
  ldas_outfiles->Air_Temp = open_file(make_fname(fname->result_dir,"AirTemp_",  "",dmy),mode);
  ldas_outfiles->Wind     = open_file(make_fname(fname->result_dir,"Wind_",     "",dmy),mode);
  ldas_outfiles->RHumid   = open_file(make_fname(fname->result_dir,"RHumid_",   "",dmy),mode);
      
  /* open filenames which are required for each layer */
  for(i=0;i<options.Nlayer;i++){
    
    strcpy(prefix,"SoilTemp_");
    sprintf(cindex, "%1.1d", i+1);
    strcat(prefix,cindex);
    strcat(prefix,"_");
    ldas_outfiles->SoilTemp[i] = open_file(make_fname(fname->result_dir,prefix, "", dmy),mode);
      
    strcpy(prefix,"SoilMoist_");
    sprintf(cindex, "%1.1d", i+1);
    strcat(prefix,cindex);
    strcat(prefix,"_");
    ldas_outfiles->SoilMoist[i] = open_file(make_fname(fname->result_dir,prefix, "", dmy),mode);
      
    strcpy(prefix,"LSoilMoist_");
    sprintf(cindex, "%1.1d", i+1);
    strcat(prefix,cindex);
    strcat(prefix,"_");
    ldas_outfiles->LSoilMoist[i] = open_file(make_fname(fname->result_dir,prefix, "", dmy),mode);
  }

}
