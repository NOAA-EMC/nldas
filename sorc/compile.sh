set -x
date

module purge
module use .
module load Module_nldas_v2.5.0

cd cpc_precip_convert.fd
make clean
make

cd ../nldas_precip.fd
make clean
make

cd ../nldas_namrr.fd
cd MAKDEP
gmake clean
gmake
cd ..
make clean
make

cd ../nldas_nam.fd
cd MAKDEP
gmake clean
gmake
cd ..
gmake clean
gmake

cd ../nldas_prep.fd
sh compile.sh

cd ../nldas_sac_ldas.fd
make clean
make

cd ../nldas_rout.fd
make clean
make

cd ../nldas_mosaic_ldas.fd
sh compile.sh


cd ../nldas_noah_ldas.fd
make clean
make

cd ../nldas_vic_prep.cd
sh compile.sh


cd ../nldas_vic_ldas.cd
make clean
make

cd ../nldas_vic_post.fd
sh compile.sh

cd ..
date

