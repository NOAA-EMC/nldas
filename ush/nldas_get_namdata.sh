#!/bin/ksh
############################################################################
# This script gets nam forcing data (precip, 10-m wind speed, 2-m air
# temperature, 2-m specific humidity, surface pressure, dowrad shoert- and 
# long-wave rdaition) to work directory.
# These data will be directly used for NLDAS forcing to close time lag as they
# are hourly and NLDAS grid and doamin.
# History: original script developed by Yuqiu Zhu for Noah model run, and it
#          was modified by Youlong Xia for nldas data preparation.
# 09 September 2016 ---- Xia      
############################################################################
set -x

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage:  start-date end-date"
  $DATA/err_exit 99
fi

sdate=$1
edate=$2

# set NAM COMING directory
# COMIN_NAM was setup in configuration file

while [ $sdate -le $edate ]; do
   year=`echo $sdate |cut -c1-4`
   today=$sdate
   yearmon=`echo $sdate |cut -c1-6` 

   # set up work directory for nam 
   export fnam=$DATA/ndas

   if [ -s $fnam ]; then
   echo $fnam exists
   else
   mkdir -p $fnam
   fi
# set date for tomorrow 
   tomorrow=`finddate.sh $today d+1`

# set date for the day after tomorrow
   day2fcst=`finddate.sh $today d+2`

#  Get nam data and copy them into correct time
#  from current 00z nam cycle to one forecast (13-23 hour) forecast 
#  to avoid 0-12 hour forecast jump from 12Z nam cycle - to connect NDAS data 
   for hh in 13 14 15 16 17 18 19 20 21 22 23
   do
   
   case $hh in
   13) hh1=13;ymd=$today;;
   14) hh1=14;ymd=$today;;
   15) hh1=15;ymd=$today;;
   16) hh1=16;ymd=$today;;
   17) hh1=17;ymd=$today;;
   18) hh1=18;ymd=$today;;
   19) hh1=19;ymd=$today;;
   20) hh1=20;ymd=$today;;
   21) hh1=21;ymd=$today;;
   22) hh1=22;ymd=$today;;
   23) hh1=23;ymd=$today;;
   esac
   export val_ndas=$COMIN_NAM/nam.${today}/nam.t00z.awldas${hh}.tm00.grib2
#  convert grib2 to grib1
    cp $val_ndas $DATA/ndas/${ymd}${hh1}.ndas
   $CNVGRIB -g21 $DATA/ndas/${ymd}${hh1}.ndas $DATA/ndas/${ymd}${hh1}.NDAS
     rm -rf $DATA/ndas/${ymd}${hh1}.ndas
   done

#  from current 00z nam cycle to one forecast (24-36 hour) forecast

   for hh in 24 25 26 27 28 29 30 31 32 33 34 35 36
   do

   case $hh in
   
   24) hh1=00;ymd=$tomorrow;;
   25) hh1=01;ymd=$tomorrow;;
   26) hh1=02;ymd=$tomorrow;;
   27) hh1=03;ymd=$tomorrow;;
   28) hh1=04;ymd=$tomorrow;;
   29) hh1=05;ymd=$tomorrow;;
   30) hh1=06;ymd=$tomorrow;;
   31) hh1=07;ymd=$tomorrow;;
   32) hh1=08;ymd=$tomorrow;;
   33) hh1=09;ymd=$tomorrow;;
   34) hh1=10;ymd=$tomorrow;;
   35) hh1=11;ymd=$tomorrow;;
   36) hh1=12;ymd=$tomorrow;;
   esac
   export val_ndas=$COMIN_NAM/nam.${today}/nam.t00z.awldas${hh}.tm00.grib2
#  convert grib2 to grib1
    cp $val_ndas $DATA/ndas/${ymd}${hh1}.ndas
   $CNVGRIB -g21 $DATA/ndas/${ymd}${hh1}.ndas $DATA/ndas/${ymd}${hh1}.NDAS
      rm -rf $DATA/ndas/${ymd}${hh1}.ndas
   done

   sdate=`finddate.sh $sdate d+1`
   done

