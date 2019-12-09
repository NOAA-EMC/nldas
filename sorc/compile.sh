#!/bin/sh
set -x
date

module purge
module use .
module load EnvVars/1.0.2
module load ips/18.0.1.163
module load lsf/10.1
module load impi/18.0.1
module load libpng/1.2.59
module load w3emc/2.3.0
module load w3nco/2.0.6
module load bacio/2.0.2
module load g2/3.1.0
module load zlib/1.2.11
module load jasper/1.900.1
module load sigio/2.1.0
module load sp/2.0.2
module load ip/3.0.1
module load NetCDF/3.6.3
module use -a /gpfs/dell1/nco/ops/nwpara/modulefiles/compiler_prod/ips/18.0.1
module list

echo " nldas.v2.0.6 compiling starts "

cd ./cpc_precip_convert.fd
make clean
make

cd ../nldas_prep.fd
sh compile.sh

cd ../nldas_sac_ldas.fd
mkdir -p -m 775 obj
make clean
make

cd ../nldas_rout.fd
make clean
make

cd ../nldas_mosaic_ldas.fd
sh compile.sh


cd ../nldas_noah_ldas.fd
mkdir -p -m 775 OBJ
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

echo " nldas.v2.0.6 comipling is completed "
date

