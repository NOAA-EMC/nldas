#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include "vicNl_ldas.h"

/**** Open the LDAS met data files 

  Modifications
	092000	Open EDAS precip met file as well	
								JS

****/
void open_ldas_metfiles(ldas_metfiles_struct *ldas_metfiles, 
			filenames_struct     *fname,
			dmy_struct           *dmy,
			char                 *mode)
{
  ldas_metfiles->apcp  = open_file(make_fname(fname->forcing[0],"APCP_", "",dmy),mode);
  ldas_metfiles->dlwrf = open_file(make_fname(fname->forcing[0],"DLWRF_","",dmy),mode);
  ldas_metfiles->dswrf = open_file(make_fname(fname->forcing[0],"DSWRF_","",dmy),mode);
  ldas_metfiles->pres  = open_file(make_fname(fname->forcing[0],"PRES_", "",dmy),mode);
  ldas_metfiles->spfh  = open_file(make_fname(fname->forcing[0],"SPFH_", "",dmy),mode);
  ldas_metfiles->tmp   = open_file(make_fname(fname->forcing[0],"TMP_",  "",dmy),mode);
  ldas_metfiles->ugrd  = open_file(make_fname(fname->forcing[0],"UGRD_", "",dmy),mode);
  ldas_metfiles->vgrd  = open_file(make_fname(fname->forcing[0],"VGRD_", "",dmy),mode);
/** Start of changes -- JS **/
  /* open EDAS precip file */
  //ldas_metfiles->apcp1  = open_file(make_fname(fname->forcing[0],"APCP1_", "",dmy),mode);
/** End of changes -- JS **/
}
