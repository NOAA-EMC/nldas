set -x
date

cd cpc_precip_convert.fd
make clean
make

cd ../nldas_prep.fd
#sh comp.csh
sh compile.sh

cd ../nldas_sac_ldas.fd
make clean
make

cd ../nldas_rout.fd
make clean
make

cd ../nldas_mosaic_ldas.fd
cd make
gmake clean
cd MAKDEP
gmake clean
cd ../..
sh comp.csh


cd ../nldas_noah_ldas.fd
make clean
make

cd ../vic_prep.cd
sh compile.sh


cd ../vic_core.cd
make clean
make

cd ../nldas_vic_post.fd
sh compile.sh

cd  ../nldas_vic_prep.cd
sh compile.sh


cd ..
date

