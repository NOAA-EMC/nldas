#!/bin/ksh
############################################################################
# This script gets cpc gauge precipitation for conus (0.125 deg resolution)
# and globe (0.25 degree resolution), 8 km cpc CMORPH precipitation, and rcads
# forcing data (precip, 10-m wind speed, 2-m air ytemperature, 2-m specific
# humidity, surface pressure, dowrad shoert- and long-wave rdaition) to work
# directory. These data will be used for generating nldas forcing using a method
# developed by Cosgrove et al. (2003).
# History: original script developed by Yuqiu Zhu for Noah model run, and it
#          was modified by Youlong Xia for nldas data preparation.      
############################################################################
set -x

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_noah.sh start-date end-date"
  $DATA/err_exit 99
fi

sdate=$1
edate=$2

# set wgrib command
# NOTE: WGRIB was loaded with the grib_util module

while [ $sdate -le $edate ]; do
   year=`echo $sdate |cut -c1-4`
   today=$sdate
   yearmon=`echo $sdate |cut -c1-6` 

   #get cpc gauge conus and globe and save into $WORKDIR

   export cpc_conus_in=PRCP_CU_GAUGE_V1.0CONUS_0.125deg.lnx.$today.RT
   export cpc_globe_in=PRCP_CU_GAUGE_V1.0GLB_0.50deg.lnx.$today.RT
   export val_conus_in=$DCOM_IN/$today/wgrbbul/cpc_rcdas/$cpc_conus_in
   export val_globe_in=$DCOM_IN/$today/wgrbbul/cpc_rcdas/$cpc_globe_in
   
   export pgm=cpc_precip_convert
   . prep_step

   if [ -s $val_conus_in ]; then
   cp -p $val_conus_in $DATA/cpcops.bin
   fi

   if [ -s $val_globe_in ]; then
   cp -p $val_globe_in $DATA/cpc_global.bin
   fi

   echo 'converting cpc conus and global precip to nldas domain and format'
   startmsg
   $EXECnldas/cpc_precip_convert >> $pgmout 2>errfile
   export err=$?; $DATA/err_chk

   #To setup directory for date stored 
   export pconus=$DATA/PRECIP.PRISM.NEW/$yearmon
   export pglobe=$DATA/PRECIP.HIGGINS.MEX/$yearmon

   if [ -s $pconus ]; then
   echo $pconus exists
   else
   mkdir -p $pconus
   fi

   if [ -s $pglobe ]; then
   echo $pglobe exists
   else
   mkdir -p $pglobe
   fi

   cp  -p cpcldas.bin $pconus/${today}.higgin.prism
   cp  -p cpc_globe.bin $pglobe/${today}.ll

   #get stage II precipitation data and put them into $DATA
   export pstage2=$DATA/PRECIP.STAGE2/$yearmon

   if [ -s $pstage2 ]; then
   echo $pstage2 exists
   else
   mkdir -p $pstage2
   fi

   for hh in 00 01 02 03 04 05 06 07 08 09 10 11 12
   do
   gunzip -v -c $COMINpcpanl/pcpanl.$today/ST2ml${today}${hh}.Grb.gz >$pstage2/ST2ml${today}${hh}.Grb
   done

   for hh in 13 14 15 16 17 18 19 20 21 22 23
   do
   gunzip -v -c $COMINpcpanl/pcpanl.$today/ST2ml${today}${hh}.Grb.gz >$pstage2/ST2ml${today}${hh}.Grb  
   done

   #get CMORPH data
   export dir_cmorph=$DCOM_IN/$today/wgrbbul/cpc_rcdas
   export pcmorph=$DATA/CMORPH/$yearmon

   if [ -s $pcmorph ]; then
   echo $pcmorph exists
   else
   mkdir -p $pcmorph
   fi

   
   hh=00
   while [ $hh -le 23 ]
   do
     # JY - old name 201511: val_cmorph=advt-8km-interp-prim-sat-spat-2lag-2.5+5dovlp8kmIR-${today}${hh}
     val_cmorph=CMORPH_8KM-30MIN_${today}${hh}
     val_cmorph_new=advt-8km-intrp-prim-sat-spat-2lag-2.5+5dovlp8kmIR-${today}${hh}

     if [ -s $dir_cmorph/$val_cmorph ]; then
     cp -p $dir_cmorph/$val_cmorph $pcmorph/$val_cmorph_new
     fi

     hh=`expr $hh + 1`
     if [ $hh -lt 10 ]; then hh=0$hh; fi
   done

   # get rcdas forcing for nldas forcing generation
   # set up work directory for NARR 
   export fnarr=$DATA/NARR/$year/$today

   if [ -s $fnarr ]; then
   echo $fnarr exists
   else
   mkdir -p $fnarr
   fi

   for hh in 00 03 06 09 12 15 18 21
   do
   export val_rcdas=$COMIN_RCDAS/rcdas.${today}/rcdas.t${hh}z.awip32.merged
   export newfile=$DATA/NARR/$year/$today/${today}${hh}.NARR.grb
   cp $val_rcdas $DATA/merged_AWIP32.${today}${hh}
   export bigfile=$DATA/merged_AWIP32.${today}${hh}

   # strip out nldas required records
   $WGRIB $bigfile | grep "kpds5=1:kpds6=109:kpds7=1"  | $WGRIB $bigfile -i -grib -o $newfile
   $WGRIB $bigfile | grep "kpds5=7:kpds6=109:kpds7=1"  | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=11:kpds6=109:kpds7=1" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=51:kpds6=109:kpds7=1" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=33:kpds6=109:kpds7=1" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=34:kpds6=109:kpds7=1" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=221:kpds6=1:kpds7=0"  | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=134:kpds6=1:kpds7=0"  | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=11:kpds6=105:kpds7=2" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=51:kpds6=105:kpds7=2" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=33:kpds6=105:kpds7=10" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=34:kpds6=105:kpds7=10" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=253:kpds6=1:kpds7=0" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=208:kpds6=1:kpds7=0" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=157:kpds6=116:kpds7=46080" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=202:kpds6=1:kpds7=0" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=63:kpds6=1:kpds7=0"  | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=140:kpds6=1:kpds7=0" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=228:kpds6=1:kpds7=0" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=204:kpds6=1:kpds7=0" | $WGRIB $bigfile -i -grib -append -o $newfile
   $WGRIB $bigfile | grep "kpds5=205:kpds6=1:kpds7=0" | $WGRIB $bigfile -i -grib -append -o $newfile

   rm -f $bigfile
   mv -f $newfile $fnarr/

   done

   # Make HPD precipitation directory (used in retrospective rather than realtime) 
   export hpd=$DATA/PRECIP.HPD
   if [ -s $hpd ]; then
   echo $hpd exists
   else
   mkdir -p $hpd
   fi

   sdate=`finddate.sh $sdate d+1`
   done
