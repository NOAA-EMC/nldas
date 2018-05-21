cd ./Libraries/grib_v3.0.0
export GRIB_ENV=`pwd`
./Clean
./Install linux_gcc
cd ../..
make clean
make OS=linux_ia64 LOCATION=PU FORCING=NARR_RETRO
