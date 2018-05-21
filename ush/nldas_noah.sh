#!/bin/ksh
##############################################################
# This script performs the Noah land-surface analysis
# Usage:   nldas_noah.sh start-date end-date
# History:   2013.05  Youlong Xia  Original script in C-shell
#            2013.06  Julia Zhu   Converted to Korn Shell
##############################################################
set -x

export alert_type=${alert_type:-NLDAS_NOAH};

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_noah.sh start-date end-date"
  $DATA/err_exit 99
fi

sdate=$1
edate=$2

while [ $sdate -le $edate ]; do
  julday=`$HOMEnldas/util/ush/julian.sh $sdate`
  year=`echo $sdate |cut -c1-4`
  year1=`/nwprod/util/ush/finddate.sh $sdate d-1 |cut -c1-4`

  today=`/nwprod/util/ush/finddate.sh $sdate d+1`
  year2=`echo $today |cut -c1-4`

  echo $julday >julday
  echo $sdate >yesterday
  echo $today >today
  echo 

  export pgm=nldas_noah_ldas
  . $DATA/prep_step

  #Input files:
  echo 'linking land surface characteristics'
  export FORT50=$FIXnldas/noah/UMDunifiedmask19990327.bin
  export FORT51=$FIXnldas/noah/noah_soil_layers19990327.bin
  export FORT53=$FIXnldas/noah/noah_soiltype20040701.bin
  export FORT54=$FIXnldas/noah/UMDpveg20060619.bin
  export FORT55=$FIXnldas/noah/noah_slopetype19990325.bin
  export FORT56=$FIXnldas/noah/noah_soil_depth19990330.bin
  export FORT57=$FIXnldas/noah/noah_soil_lower_T19990325.bin
  export FORT58=$FIXnldas/noah/albedo_ldas.bin
  export FORT59=$FIXnldas/noah/noah_gfrac19990325.bin
  export FORT60=$FIXnldas/noah/max_snow_albedo.bin
  export FORT61=$FIXnldas/noah/precipmask.2deg.bin
  export FORT62=$FIXnldas/noah/noah_soiltype19990808.bin

  echo 'linking initial conditions'
  export RESDIR=${RESDIR:-$COM_IN}
  export FORT69=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_LSTSNW
  export FORT70=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_T1
  export FORT71=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_STC
  export FORT72=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_SMC
  export FORT73=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_SH2O
  export FORT74=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_CMC
  export FORT75=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_SNOWH
  export FORT76=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_SNEQV
  export FORT77=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_CH
  export FORT78=$RESDIR/nldas.$sdate/noah.t${cyc}z.INIT_CM

  #echo 'linking forcing files'
  export COMFORCE=${COMFORCE:-$COM_IN}
  export FORT01=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf00
  export FORT02=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf01
  export FORT03=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf02
  export FORT04=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf03
  export FORT05=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf04
  export FORT06=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf05
  export FORT07=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf06
  export FORT08=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf07
  export FORT09=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf08
  export FORT10=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf09
  export FORT11=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf10
  export FORT12=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf11
  export FORT13=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf12
  export FORT14=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf13
  export FORT15=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf14
  export FORT16=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf15
  export FORT17=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf16
  export FORT18=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf17
  export FORT19=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf18
  export FORT20=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf19
  export FORT21=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf20
  export FORT22=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf21
  export FORT23=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf22
  export FORT24=$COMFORCE/nldas.$sdate/nldas.t${cyc}z.force-a.grbf23
  
  sdate1=`/nwprod/util/ush/finddate.sh $sdate d+1`
  export FORT25=$COMFORCE/nldas.$sdate1/nldas.t${cyc}z.force-a.grbf00
  
  # Output files:
  export FORT80=${sdate1}00.INIT_T1
  export FORT81=${sdate1}00.INIT_STC
  export FORT82=${sdate1}00.INIT_SMC
  export FORT83=${sdate1}00.INIT_SH2O
  export FORT84=${sdate1}00.INIT_CMC
  export FORT85=${sdate1}00.INIT_CH
  export FORT86=${sdate1}00.INIT_CM
  export FORT87=${sdate1}00.INIT_SNOWH
  export FORT88=${sdate1}00.INIT_SNEQV
  export FORT90=${sdate1}00.INIT_LSTSNW
  
  # Build the control file
  cat <<EOF >controlfile
#Time step in seconds
900.0
#Forcing height in meters
6.0
#Ice flag
0
EOF

  # Copy the KPDS table:
  cp $FIXnldas/noah/KPDS.tbl .
  
  echo 'running NOAH model'
  $EXECnldas/nldas_noah_ldas
  export err=$?; $DATA/err_chk
  
  chmod 644 *grb
  if [ ${SENDCOM:-YES} = YES ]; then
    # Copy restart files to /com:
    export COMREST=${COMREST:-$COM_OUT}
    mv ${sdate1}00.INIT_T1 $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_T1
    mv ${sdate1}00.INIT_STC $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_STC
    mv ${sdate1}00.INIT_SMC  $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_SMC
    mv ${sdate1}00.INIT_SH2O  $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_SH2O
    mv ${sdate1}00.INIT_CMC  $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_CMC
    mv ${sdate1}00.INIT_CH  $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_CH
    mv ${sdate1}00.INIT_CM  $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_CM
    mv ${sdate1}00.INIT_SNOWH  $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_SNOWH
    mv ${sdate1}00.INIT_SNEQV  $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_SNEQV
    mv ${sdate1}00.INIT_LSTSNW  $COMREST/nldas.$sdate1/noah.t${cyc}z.INIT_LSTSNW
  
    # Copy grib files to /com:
    mv ${sdate1}00.NOAH.grb  $COM_OUT/nldas.$sdate1/noah.t${cyc}z.grbf00
    if [ $SENDDBN = YES ]; then
       $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COM_OUT/nldas.$sdate1/noah.t${cyc}z.grbf00
    fi
    hh=01
    while [ $hh -le 23 ]; do
      mv ${sdate}${hh}.NOAH.grb $COM_OUT/nldas.$sdate/noah.t${cyc}z.grbf${hh}
      if [ $SENDDBN = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COM_OUT/nldas.$sdate/noah.t${cyc}z.grbf${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done
  fi
 
  echo $sdate >$LOGDIR/nldas_noah_lastprocessed.log
  sdate=`/nwprod/util/ush/finddate.sh $sdate d+1`
done
