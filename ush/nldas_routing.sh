#!/bin/ksh
##############################################################
# This script performs the river routing run for model
# Usage:   nldas_noah_routing.sh start-date end-date
# History:   2013.05  Youlong Xia:  Original script in C-shell
#            2013.08  Youlong Xia and Yuqiu Zhu:   Converted to Korn Shell
#            2018.02  Youlong Xia: Modify the script to ouput all files into
#                     single directory for all 7 days in NLDAS v2.5 (use date
#                     in filename as suggested by NCO staff to identify 
#                     different files)
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

  ydate=`finddate.sh $sdate d-1`
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
  export COM_IN=${COM_IN:-$COMROOT/nldas/${envir}}
  export INITDIR0=${INITDIR0:-$COM_IN}
  export INITDIR=$INITDIR0/nldas.$PDY

  export inter_file=${model}.t${cyc}z.${sdate}.internal_runoff.bin
  export trans_file=${model}.t${cyc}z.${sdate}.transport_runoff.bin 
  export inter_in=$INITDIR/${inter_file}
  export trans_in=$INITDIR/${trans_file}

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
  export INDIR0=${INDIR0:-$COM_IN}
  export INDIR=$INDIR0/nldas.$PDY
  export FORT01=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f00
  export FORT02=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f01
  export FORT03=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f02
  export FORT04=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f03
  export FORT05=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f04
  export FORT06=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f05
  export FORT07=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f06
  export FORT08=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f07
  export FORT09=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f08
  export FORT10=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f09
  export FORT11=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f10
  export FORT12=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f11
  export FORT13=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f12
  export FORT14=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f13
  export FORT15=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f14
  export FORT16=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f15
  export FORT17=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f16
  export FORT18=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f17
  export FORT19=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f18
  export FORT20=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f19
  export FORT21=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f20
  export FORT22=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f21
  export FORT23=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f22
  export FORT24=$INDIR/nldas.t${cyc}z.${sdate}.${model}.grb2f23

  sdate1=`/nwprod/util/ush/finddate.sh $sdate d+1`
  
  # Output files:
  export FORT81=${sdate1}.internal_runoff.bin
  export FORT82=${sdate1}.transport_runoff.bin
  

  # Copy the KPDS table
  cp $FIXnldas/routing/MASKS/KPDS.tbl .
  export FORT30=$DATA/KPDS.tbl 
 
  echo "running river routing model"
  startmsg
  $EXECnldas/nldas_rout >> $pgmout 2>errfile
  export err=$?; $DATA/err_chk
  
  chmod 644 *grb
  if [ ${SENDCOM:-YES} = YES ]; then
    # Copy restart files to /com:
    export COMREST0=${COMREST0:-$COM_OUT}
    export COMT1=$COMREST0/nldas.$PDYp1
# copy 00z routing initials to tomottow diretcory and run
    if [ $sdate1 = $PDYm3 ]; then    
    cp ${sdate1}.internal_runoff.bin $COMT1/${model}.t${cyc}z.${sdate1}.internal_runoff.bin
    cp ${sdate1}.transport_runoff.bin $COMT1/${model}.t${cyc}z.${sdate1}.transport_runoff.bin
    fi

    export COMREST=$COMREST0/nldas.$PDY
    mv ${sdate1}.internal_runoff.bin $COMREST/${model}.t${cyc}z.${sdate1}.internal_runoff.bin
    mv ${sdate1}.transport_runoff.bin $COMREST/${model}.t${cyc}z.${sdate1}.transport_runoff.bin
  
    # Copy grib files to /com:
    hh=00
    while [ $hh -le 23 ]; do
      mv ${sdate}${hh}.STRM.grb  $COMREST/nldas.t${cyc}z.${sdate}.${model}.STRM.grb2f${hh}
      if [ $SENDDBN = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job  $COMREST/nldas.t${cyc}z.${sdate}.${model}.STRM.grb2f${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done
  fi
 
  sdate=`finddate.sh $sdate d+1`
done
