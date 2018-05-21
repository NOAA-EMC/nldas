set -x
date

EXECDIR=../exec

cp -p cpc_precip_convert.fd/cpc_precip_convert  $EXECDIR

cp -p nldas_prep.fd/nldas_prep $EXECDIR

cp -p nldas_precip.fd/precforce  $EXECDIR

cp -p nldas_namrr.fd/mergeforce $EXECDIR

cp -p nldas_nam.fd/namforecast $EXECDIR

cp -p nldas_sac_ldas.fd/nldas_sac_ldas  $EXECDIR

cp -p nldas_rout.fd/nldas_rout  $EXECDIR

cp -p nldas_mosaic_ldas.fd/make/nldas_mosaic_ldas  $EXECDIR

cp -p nldas_noah_ldas.fd/nldas_noah_ldas $EXECDIR

cp -p nldas_vic_prep.cd/nldas_vic_prep  $EXECDIR

cp -p nldas_vic_ldas.cd/nldas_vic_ldas  $EXECDIR

cp -p nldas_vic_post.fd/nldas_vic_post  $EXECDIR

date

