#!/bin/ksh
##############################################################
# This script performs the Mosaic land-surface analysis
# Usage:   nldas_mosaic.sh start-date end-date
# History:   2013.05  Charles Alonge and Youlong Xia  Original script 
#                     in C-shell
#            2013.06  Youlong Xia  Converted to Korn Shell
#            2019.03  Youlong Xia  Use modular for mosaic run
##############################################################
set -x

export alert_type=${alert_type:-NLDAS_MOSAIC};

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_mosaic.sh start-date end-date"
  err_exit 99
fi

sdate=$1
edate=$2

while [ $sdate -le $edate ]; do
  year1=`echo $sdate |cut -c1-4`
  year2=`finddate.sh $sdate d+1 |cut -c1-4`

  day1=$sdate
  day2=`finddate.sh $sdate d+1`
  echo 

  # Copy initial conditions for mosaic model
  cp $COM_IN/nldas.$sdate/mosaic.t${cyc}z.MOSrst $DATA/mos.rst
  
  # Copy lis template card
  cp $PARMnldas/lis-template.crd $DATA/lis-template.crd
  
  # BUILD the ldas.crd file for a realtime update
  echo $day1
  export SDAY=$day1
  export EDAY=$day2
  echo $SDAY $EDAY

  export SDA=`echo $SDAY | cut -c 7-8`
  export SMO=`echo $SDAY | cut -c 5-6`
  export SYR=`echo $SDAY | cut -c 1-4`

  export EDA=`echo $EDAY | cut -c 7-8`
  export EMO=`echo $EDAY | cut -c 5-6`
  export EYR=`echo $EDAY | cut -c 1-4`

  #UPDATE TIMES IN TEMPLATE CARD
  sed "s?ZDA?$SDA?" $DATA/lis-template.crd | sed "s?ZMO?$SMO?" | sed "s?ZYR?$SYR?"| sed "s?XDA?$EDA?"| sed "s?XMO?$EMO?"| sed "s?XYR?$EYR?" > $DATA/lis.crd
 
  # link fix fileds to $DATA directory
  ln -s $FIXnldas/mosaic $DATA/mosaic
 
  # copy grib table to $DATA
  cp -p $FIXnldas/mosaic/tables/KPDS_completemos.tbl $DATA/.
 
  # copy nlddas forcing to $DATA directory for $sdate and 
  #rename file name to fit source code

  export asufix=nldasforce-a.grb
  export narrdir=$DATA/NLDASNARR/$year1/$sdate
  mkdir -p $narrdir

  hh=00
  while [ $hh -le 23 ]; do
    export afile1=$COM_IN/nldas.$sdate/nldas.t${cyc}z.force-a.grbf${hh}
    cp $afile1 $narrdir/${sdate}${hh}.${asufix}
    let "hh=hh+1"
    if [ $hh -lt 10 ]; then hh=0$hh; fi
  done
 
  # copy nlddas forcing to $DATA directory for $sdate1 and rename
  # file name to fit source code
  sdate1=`finddate.sh $sdate d+1`
  export narrdir2=$DATA/NLDASNARR/$year2/$sdate1
  mkdir -p $narrdir2

  hh=00
  while [ $hh -le 12 ]; do
    export afile21=$COM_IN/nldas.$sdate1/nldas.t${cyc}z.force-a.grbf${hh}
    cp $afile21 $narrdir2/${sdate1}${hh}.${asufix}
    let "hh=hh+1"
    if [ $hh -lt 10 ]; then hh=0$hh; fi
  done
  
  ## complete nldas foricng setup
  export outforce=$DATA/OUTPUT/EXP888/MOS
  mkdir -p $outforce
   
  # Run mosaic land model with the updated LIS nldas card (ldas.crd)
  export pgm=nldas_mosaic_ldas
  . prep_step

  startmsg
  $EXECnldas/nldas_mosaic_ldas >> $pgmout 2>>errfile
  export err=$?; err_chk
   
  # Archive model initials
  export COMREST=${COMREST:-$COM_OUT}
  export INDIR=$DATA/OUTPUT/EXP888/MOS 
  mv $INDIR/$year2/$sdate1/LIS.E888.${sdate1}12.MOSrst $COMREST/nldas.$sdate1/mosaic.t${cyc}z.MOSrst

   # Copy sdate1 (0-12Z) mosaic model ouput files to /com:

    hh=00
    while [ $hh -le 12 ]; do
      cp -p $INDIR/$year2/$sdate1/${sdate1}${hh}.grb $COM_OUT/nldas.$sdate1/mosaic.t${cyc}z.grbf${hh}
      if [ $SENDDBN = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COM_OUT/nldas.$sdate1/mosaic.t${cyc}z.grbf${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done
  
  #copy sdate nldas frocing file (13z - 23z) to /com

   hh=13
   while [ $hh -le 23 ]; do
     cp  $INDIR/$year1/$sdate/${sdate}${hh}.grb $COM_OUT/nldas.$sdate/mosaic.t${cyc}z.grbf${hh}
     if [ $SENDDBN = YES ]; then
        $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COM_OUT/nldas.$sdate/mosaic.t${cyc}z.grbf${hh}
     fi
     let "hh=hh+1"
     if [ $hh -lt 10 ]; then hh=0$hh; fi
   done
 
  sdate=`finddate.sh $sdate d+1`
done

