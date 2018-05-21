#!/bin/ksh
##############################################################
# This script performs the nldas realtime orcing generation process
# Usage:   nldas_prep_realtime.sh start-date end-date
# History:   2013.05  Youlong Xia  Original script in C-shell
#            2013.06  Youlong Xia  Converted to Korn Shell
#            2016.02  Youlong Xia  Preapred for actual realtime run 
##############################################################
set -x

export alert_type=${alert_type:-NLDAS_PREP_REALTIME}

# using $PDY TO MAKE OUTPUT DIRECTORY

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_prep_ndas.sh start-date end-date"
  $DATA/err_exit 99
fi

sdate=$1
edate=$2

export cyc=${cyc:-12}

while [ $sdate -le $edate ]; do
  year1=`echo $sdate |cut -c1-4`
  year0=`finddate.sh $sdate d-1 |cut -c1-4`

  day1=$sdate
  day0=`finddate.sh $sdate d-1`

  export pgm=nldas_prep_realtime
  . $DATA/prep_step
  
  echo 'Here we start to run nldas actual realtime forcing generation'
  # As forcing generation algorithm needs total two complete days to ensure
  # disaggregated precipitation internal consistency, we need set start date to
  # two days prior $sdate
  
  cp -p $FIXnldas/prep/tables/KPDS.force.ndas_nam.tbl $DATA
  cp -p $FIXnldas/prep/tables/KPDS.force.ndas_nam.grib2.tbl $DATA  
  export FORT60=$DATA/KPDS.force.ndas_nam.tbl
  export FORT70=$DATA/KPDS.force.ndas_nam.grib2.tbl

#  set tempory NAMX directory
   export NAMXD=$DATA/ndas

#  Merge all processed data into a unified forcing data
   export FORT25=$NAMXD/${day0}13.NDAS
   export FORT26=$NAMXD/${day0}14.NDAS
   export FORT27=$NAMXD/${day0}15.NDAS
   export FORT28=$NAMXD/${day0}16.NDAS
   export FORT29=$NAMXD/${day0}17.NDAS
   export FORT30=$NAMXD/${day0}18.NDAS
   export FORT31=$NAMXD/${day0}19.NDAS
   export FORT32=$NAMXD/${day0}20.NDAS
   export FORT33=$NAMXD/${day0}21.NDAS
   export FORT34=$NAMXD/${day0}22.NDAS
   export FORT35=$NAMXD/${day0}23.NDAS
   export FORT36=$NAMXD/${day1}00.NDAS
   export FORT37=$NAMXD/${day1}01.NDAS
   export FORT38=$NAMXD/${day1}02.NDAS
   export FORT39=$NAMXD/${day1}03.NDAS
   export FORT40=$NAMXD/${day1}04.NDAS
   export FORT41=$NAMXD/${day1}05.NDAS
   export FORT42=$NAMXD/${day1}06.NDAS
   export FORT43=$NAMXD/${day1}07.NDAS
   export FORT44=$NAMXD/${day1}08.NDAS
   export FORT45=$NAMXD/${day1}09.NDAS
   export FORT46=$NAMXD/${day1}10.NDAS
   export FORT47=$NAMXD/${day1}11.NDAS
   export FORT48=$NAMXD/${day1}12.NDAS

   export FORT99=$FIXnldas/prep/mask/UMD_Land_Sea_Mask.bin 
#  set up date for forcing output    
   echo $day0 >beforeyesterday
   echo $day1 >yesterday

echo 'merge all processed hourly data into a unfied file'
   export pgm=namforecast
   
   startmsg
   $EXECnldas/namforecast >> $pgmout 2>errfile                #1a
   export err=$?; err_chk
	
  export COMFOR=$COM_OUT/nldas.$PDY
  if [ -s $COMFOR ]; then
   echo $COMFOR exists
   else
   mkdir -p $COMFOR
   fi
    
    export asufix=forcing.grb
    hh=13
    while [ $hh -le 23 ]; do
    cp -rp ${day0}${hh}.${asufix} $COMFOR/nldas.t${cyc}z.${day0}.force-a.grbf${hh}
    cp -rp ${day0}${hh}.${asufix}2 $COMFOR/nldas.t${cyc}z.${day0}.force-a.grb2f${hh}

      if [ $SENDDBN = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR/nldas.t${cyc}z.${day0}.force-a.grbf${hh}
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR/nldas.t${cyc}z.${day0}.force-a.grb2f${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done

  cp -rp ${day1}00.${asufix}  $COMFOR/nldas.t${cyc}z.${day1}.force-a.grbf00
  cp -rp ${day1}00.${asufix}2  $COMFOR/nldas.t${cyc}z.${day1}.force-a.grb2f00

  if [ $SENDDBN = YES ]; then
       $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR/nldas.t${cyc}z.${day1}.force-a.grbf00
       $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR/nldas.t${cyc}z.${day1}.force-a.grb2f00
   fi
    hh=01
    while [ $hh -le 12 ]; do
    cp -rp ${day1}${hh}.${asufix} $COMFOR/nldas.t${cyc}z.${day1}.force-a.grbf${hh}
    cp -rp ${day1}${hh}.${asufix}2 $COMFOR/nldas.t${cyc}z.${day1}.force-a.grb2f${hh}
    
      if [ $SENDDBN = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR/nldas.t${cyc}z.${day1}.force-a.grbf${hh}
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR/nldas.t${cyc}z.${day1}.force-a.grb2f${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done
    
  sdate=`finddate.sh $day1 d+1`
done
