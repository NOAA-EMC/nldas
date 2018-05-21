#!/bin/ksh
##############################################################
# This script performs the river routing run for model
# Usage:   nldas_noah_routing.sh start-date end-date
# History:   2013.05  Youlong Xia  Original script in C-shell
#            2013.08  Youlong Xia snd Yuqiu Zhu   Converted to Korn Shell
##############################################################
set -x

export alert_type=${alert_type:-NLDAS_${MODEL}}

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_noah_routing.sh start-date end-date"
  $DATA/err_exit 99
fi

sdate=$1
edate=$2

export model=${model:-noah}
export MODEL=${MODEL:-`echo $model |tr '[a-z]' '[A-Z]'`}

rm -rf  days.txt days0.txt

while [ $sdate -le $edate ]; do
  year=`echo $sdate |cut -c1-4`
  echo $year $sdate > days.txt 

  ydate=`/nwprod/util/ush/finddate.sh $sdate d-1`
  year1=`echo $ydate |cut -c1-4`
  echo $year1 $ydate > days0.txt
  echo

  export pgm=nldas_rout

  . $DATA/prep_step

  # Copy nldas land-sea mask
  echo "copy NLDAS mask to $data directory"
  export land_sea_mask=$FIXnldas/routing/MASKS/UMDunifiedmask19990327.bin
  cp $land_sea_mask .

  #Input files:
  echo "copying river routing mask files"
  export mask_dir=$FIXnldas/routing/MASKS/
  export mask_order=order.bin
  export mask_ingra=internal_hydrographs.bin
  export mask_trgra=transport_hydrographs.bin 

  cp $mask_dir${mask_order} .
  cp $mask_dir${mask_ingra} .
  cp $mask_dir${mask_trgra} .

  echo "linking initial conditions"
  export COM_IN=${COM_IN:-/com/nldas/${envir}}
  export INITDIR=${INITDIR:-$COM_IN}
  export inter_file=${model}.t${cyc}z.internal_runoff.bin
  export trans_file=${model}.t${cyc}z.transport_runoff.bin 
  export inter_in=$INITDIR/nldas.$sdate/${inter_file}
  export trans_in=$INITDIR/nldas.$sdate/${trans_file}

  cp $inter_in .
  cp $trans_in .

  if [ -s directories.txt ]; then
    rm -rf  directories.txt
  fi

  # Build text file for mask and intitials 
   echo $inter_file >>directories.txt
   echo $trans_file >>directories.txt
   echo $mask_order >>directories.txt
   echo $mask_ingra >>directories.txt
   echo $mask_trgra >>directories.txt

  echo "linking input files"
  export INDIR=${INDIR:-$COM_IN}
  export FORT01=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf00
  export FORT02=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf01
  export FORT03=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf02
  export FORT04=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf03
  export FORT05=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf04
  export FORT06=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf05
  export FORT07=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf06
  export FORT08=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf07
  export FORT09=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf08
  export FORT10=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf09
  export FORT11=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf10
  export FORT12=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf11
  export FORT13=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf12
  export FORT14=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf13
  export FORT15=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf14
  export FORT16=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf15
  export FORT17=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf16
  export FORT18=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf17
  export FORT19=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf18
  export FORT20=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf19
  export FORT21=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf20
  export FORT22=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf21
  export FORT23=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf22
  export FORT24=$INDIR/nldas.$sdate/${model}.t${cyc}z.grbf23

  sdate1=`/nwprod/util/ush/finddate.sh $sdate d+1`
  
  # Output files:
  export FORT81=${sdate1}.internal_runoff.bin
  export FORT82=${sdate1}.transport_runoff.bin
  

  # Copy the KPDS table:
  cp $FIXnldas/routing/MASKS/KPDS.tbl .
  
  echo "running river routing model"

  $EXECnldas/nldas_rout
  export err=$?; $DATA/err_chk
  
  chmod 644 *grb
  if [ ${SENDCOM:-YES} = YES ]; then
    # Copy restart files to /com:
    export COMREST=${COMREST:-$COM_OUT}
    mv ${sdate1}.internal_runoff.bin $COMREST/nldas.$sdate1/${model}.t${cyc}z.internal_runoff.bin
    mv ${sdate1}.transport_runoff.bin $COMREST/nldas.$sdate1/${model}.t${cyc}z.transport_runoff.bin
  
    # Copy grib files to /com:
    hh=00
    while [ $hh -le 23 ]; do
      mv ${sdate}${hh}.STRM.grb $COM_OUT/nldas.$sdate/${model}.t${cyc}z.STRM.grbf${hh}
      if [ $SENDDBN = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COM_OUT/nldas.$sdate/${model}.t${cyc}z.STRM.grbf${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done
  fi
 
  echo $sdate >$LOGDIR/nldas_${model}_routing.log
  sdate=`/nwprod/util/ush/finddate.sh $sdate d+1`
done
