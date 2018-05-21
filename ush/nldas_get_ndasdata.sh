#!/bin/ksh
############################################################################
# This script gets global cpc gauge precipitation (0.125 deg resolution)
# and ndas & nam forcing data (precip, 10-m wind speed, 2-m air ytemperature,
# 2-m specific humidity, surface pressure, dowrad shoert- and long-wave 
# rdaition), and Stage IV hourly precipitation to work directory.
# These data will be used for generating nldas forcing.
# History: original script developed by Yuqiu Zhu for Noah model run, and it
#          was modified by Youlong Xia for nldas data preparation.      
############################################################################
set -x

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage:  start-date end-date"
  $DATA/err_exit 99
fi

sdate=$1
edate=$2

# set cpc directory in configuration (DCOM_IN,COMINpcpanl)

# set NDAS directory in nldas_para_config of user directory

# use module for WGRIB

while [ $sdate -le $edate ]; do
   year=`echo $sdate |cut -c1-4`
   today=$sdate
   yearmon=`echo $sdate |cut -c1-6` 

   #get cpc globe gauge (0.125 degree) and save into $WORKDIR

   export cpc_globe_in=PRCP_CU_GAUGE_V1.0GLB_0.125deg.lnx.$today.RT
   export val_globe_in=$DCOM_IN/$today/wgrbbul/cpc_rcdas/$cpc_globe_in

#To setup directory for date stored 
   export pglobe=$DATA/CPC_GLOBE

   if [ -s $pglobe ]; then
   echo $pglobe exists
   else
   mkdir -p $pglobe
   fi

   if [ -s $val_globe_in ]; then
   cp  -p $val_globe_in $pglobe/${today}.ll
   else
   echo cpc daily precipitation does not exist
   fi

   #get stage IV precipitation data and put them into $DATA
   export pstage2=$DATA/PRECIP.STAGE2

   if [ -s $pstage2 ]; then
   echo $pstage2 exists
   else
   mkdir -p $pstage2
   fi

#  set up yesterday (12Z yesterday - 12Z today)
   yesterday=`finddate.sh $today d-1`

   for hh in 00 01 02 03 04 05 06 07 08 09 10 11 12
   do
   gunzip -v -c $COMINpcpanl/pcpanl.$today/ST4.${today}${hh}.01h.gz >$pstage2/ST2ml${today}${hh}.Grb
   done

   for hh in 13 14 15 16 17 18 19 20 21 22 23
   do
   gunzip -v -c $COMINpcpanl/pcpanl.$yesterday/ST4.${yesterday}${hh}.01h.gz >$pstage2/ST2ml${yesterday}${hh}.Grb  
   done

   # get ndas forcing for nldas forcing generation
   # set up work directory for ndas 
   export fndas=$DATA/ndas

   if [ -s $fndas ]; then
   echo $fndas exists
   else
   mkdir -p $fndas
   fi

#  set date for yesterday
#  00z run is for 19z,20z,21z,22zand 23z yesterday and 00z today
#  Get 00Z ndas data and copy them into correct time
   for hh in 01 02 03 04 05 06
   do

   case $hh in
   01) hh1=00;ymd=$today;;
   02) hh1=23;ymd=$yesterday;;
   03) hh1=22;ymd=$yesterday;;
   04) hh1=21;ymd=$yesterday;;
   05) hh1=20;ymd=$yesterday;;   
   06) hh1=19;ymd=$yesterday;;
   esac
   export val_ndas=$COMIN_NDAS/nam.${today}/nam.t00z.awldas01.tm${hh}.grib2
#  convert grib2 to grib1
   cp $val_ndas $DATA/ndas/${ymd}${hh1}.ndas   
   $CNVGRIB -g21 $DATA/ndas/${ymd}${hh1}.ndas $DATA/ndas/${ymd}${hh1}.NDAS
   rm -rf $DATA/ndas/${ymd}${hh1}.ndas     
   done

#  06z run is exact today
#  Get 06z ndas data and copy them into correct time
   for hh in 01 02 03 04 05 06
   do

   case $hh in
   01) hh1=06;ymd=$today;;
   02) hh1=05;ymd=$today;;
   03) hh1=04;ymd=$today;;
   04) hh1=03;ymd=$today;;
   05) hh1=02;ymd=$today;;
   06) hh1=01;ymd=$today;;
   esac
   export val_ndas=$COMIN_NDAS/nam.${today}/nam.t06z.awldas01.tm${hh}.grib2
#  convert grib2 to grib1
   cp $val_ndas $DATA/ndas/${ymd}${hh1}.ndas
   $CNVGRIB -g21 $DATA/ndas/${ymd}${hh1}.ndas $DATA/ndas/${ymd}${hh1}.NDAS
   rm -rf $DATA/ndas/${ymd}${hh1}.ndas 
   done

# 12z run is exact today
#  Get 12z ndas data and copy them into correct time
   for hh in 01 02 03 04 05 06
   do

   case $hh in
   01) hh1=12;ymd=$today;;
   02) hh1=11;ymd=$today;;
   03) hh1=10;ymd=$today;;
   04) hh1=09;ymd=$today;;
   05) hh1=08;ymd=$today;;
   06) hh1=07;ymd=$today;;
   esac
   export val_ndas=$COMIN_NDAS/nam.${today}/nam.t12z.awldas01.tm${hh}.grib2
#  convert grib2 to grib1
   cp $val_ndas $DATA/ndas/${ymd}${hh1}.ndas
   $CNVGRIB -g21 $DATA/ndas/${ymd}${hh1}.ndas $DATA/ndas/${ymd}${hh1}.NDAS
   rm -rf $DATA/ndas/${ymd}${hh1}.ndas
   done

# 18z run is exact today
#  Get 18z ndas data and copy them into correct time
   for hh in 01 02 03 04 05 06
   do

   case $hh in
   01) hh1=18;ymd=$today;;
   02) hh1=17;ymd=$today;;
   03) hh1=16;ymd=$today;;
   04) hh1=15;ymd=$today;;
   05) hh1=14;ymd=$today;;
   06) hh1=13;ymd=$today;;
   esac
   export val_ndas=$COMIN_NDAS/nam.${today}/nam.t18z.awldas01.tm${hh}.grib2
#  convert grib2 to grib1
   cp $val_ndas $DATA/ndas/${ymd}${hh1}.ndas
   $CNVGRIB -g21 $DATA/ndas/${ymd}${hh1}.ndas $DATA/ndas/${ymd}${hh1}.NDAS
   rm -rf $DATA/ndas/${ymd}${hh1}.ndas
   done
  
   sdate=`finddate.sh $sdate d+1`
   done

