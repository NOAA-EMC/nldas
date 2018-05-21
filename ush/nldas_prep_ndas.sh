#!/bin/ksh
##############################################################
# This script performs the nldas realtime orcing generation process
# Usage:   nldas_prep_realtime.sh start-date end-date
# History:   2013.05  Youlong Xia  Original script in C-shell
#            2013.06  Youlong Xia  Converted to Korn Shell
#            2016.02  Youlong Xia  Preapred for actual realtime run
#            2016.09  Youlong Xia  Prepared for NAMX update 
##############################################################
set -x

export alert_type=${alert_type:-NLDAS_PREP_REALTIME}

# Using $PDY to set up today run directory

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_prep_ndas.sh start-date end-date"
  $DATA/err_exit 99
fi

sdate=$1
edate=$2

export cyc=${cyc:-12}

while [ $sdate -le $edate ]; do
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

#  set up NDAS temporary directory
   export NAMXD=$DATA/ndas

# Temporally disagrregatting CPC daily gauge into hourly data
#  Read NDAS data from temporay directory

   export FORT11=$NAMXD/${day0}13.NDAS
   export FORT12=$NAMXD/${day0}14.NDAS
   export FORT13=$NAMXD/${day0}15.NDAS
   export FORT14=$NAMXD/${day0}16.NDAS
   export FORT15=$NAMXD/${day0}17.NDAS
   export FORT16=$NAMXD/${day0}18.NDAS
   export FORT17=$NAMXD/${day0}19.NDAS
   export FORT18=$NAMXD/${day0}20.NDAS
   export FORT19=$NAMXD/${day0}21.NDAS
   export FORT20=$NAMXD/${day0}22.NDAS
   export FORT21=$NAMXD/${day0}23.NDAS
   export FORT22=$NAMXD/${day1}00.NDAS
   export FORT23=$NAMXD/${day1}01.NDAS
   export FORT24=$NAMXD/${day1}02.NDAS
   export FORT25=$NAMXD/${day1}03.NDAS
   export FORT26=$NAMXD/${day1}04.NDAS
   export FORT27=$NAMXD/${day1}05.NDAS
   export FORT28=$NAMXD/${day1}06.NDAS
   export FORT29=$NAMXD/${day1}07.NDAS
   export FORT30=$NAMXD/${day1}08.NDAS
   export FORT31=$NAMXD/${day1}09.NDAS
   export FORT32=$NAMXD/${day1}10.NDAS
   export FORT33=$NAMXD/${day1}11.NDAS
   export FORT34=$NAMXD/${day1}12.NDAS
#  Read NCEP Stage II hourly radar precipitation data
#  set up stage II input data DIRECTORY
   export temout0=$DATA/PRECIP.STAGE2

   export FORT35=$temout0/ST2ml${day0}13.Grb
   export FORT36=$temout0/ST2ml${day0}14.Grb
   export FORT37=$temout0/ST2ml${day0}15.Grb
   export FORT38=$temout0/ST2ml${day0}16.Grb
   export FORT39=$temout0/ST2ml${day0}17.Grb
   export FORT40=$temout0/ST2ml${day0}18.Grb
   export FORT41=$temout0/ST2ml${day0}19.Grb
   export FORT42=$temout0/ST2ml${day0}20.Grb
   export FORT43=$temout0/ST2ml${day0}21.Grb
   export FORT44=$temout0/ST2ml${day0}22.Grb
   export FORT45=$temout0/ST2ml${day0}23.Grb
   export FORT46=$temout0/ST2ml${day1}00.Grb
   export FORT47=$temout0/ST2ml${day1}01.Grb
   export FORT48=$temout0/ST2ml${day1}02.Grb
   export FORT49=$temout0/ST2ml${day1}03.Grb
   export FORT50=$temout0/ST2ml${day1}04.Grb
   export FORT51=$temout0/ST2ml${day1}05.Grb
   export FORT52=$temout0/ST2ml${day1}06.Grb
   export FORT53=$temout0/ST2ml${day1}07.Grb
   export FORT54=$temout0/ST2ml${day1}08.Grb
   export FORT55=$temout0/ST2ml${day1}09.Grb
   export FORT56=$temout0/ST2ml${day1}10.Grb
   export FORT57=$temout0/ST2ml${day1}11.Grb
   export FORT58=$temout0/ST2ml${day1}12.Grb

#  Read 0.125-degree global CPC gauge daily precipitation 
#  set CPC gauge precipitation directory
   export temout2=$DATA/CPC_GLOBE
#   export FORT59=$temout2/${day1}.ll
# check if cpc daily gauge precipitation exists

   if [ -s $temout2/${day1}.ll ]; then
   cp -p $temout2/${day1}.ll $DATA/cpc_125global.bin
   else
   rm -rf $DATA/cpc_125global.bin
   fi

echo 'Temporally disagregate CPC gauge Precip from daily into houly'
   export pgm=precforce
   startmsg
   $EXECnldas/precforce > $pgmout 2>>errfile                     #1a
   export err=$?; err_chk   
 
   if [ -s ${day1}12.PREC ]; then
   echo daily precipitation is temporally disaggregated!!!
   else
   echo precforce failed!!!!
   echo The program will be terminated!!!!
   export err=9; export err; err_chk
   fi

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

#  diagregated cpc gauge precipitation with Stage II and NAMv4

   export FORT73=${day0}13.PREC
   export FORT74=${day0}14.PREC
   export FORT75=${day0}15.PREC
   export FORT76=${day0}16.PREC
   export FORT77=${day0}17.PREC
   export FORT78=${day0}18.PREC
   export FORT79=${day0}19.PREC
   export FORT80=${day0}20.PREC
   export FORT81=${day0}21.PREC
   export FORT82=${day0}22.PREC
   export FORT83=${day0}23.PREC
   export FORT84=${day1}00.PREC
   export FORT85=${day1}01.PREC
   export FORT86=${day1}02.PREC
   export FORT87=${day1}03.PREC
   export FORT88=${day1}04.PREC
   export FORT89=${day1}05.PREC
   export FORT90=${day1}06.PREC
   export FORT91=${day1}07.PREC
   export FORT92=${day1}08.PREC
   export FORT93=${day1}09.PREC
   export FORT94=${day1}10.PREC
   export FORT95=${day1}11.PREC
   export FORT96=${day1}12.PREC

   export FORT99=$FIXnldas/prep/mask/UMD_Land_Sea_Mask.bin
 
#  set up date for forcing output    
   echo $day0 >beforeyesterday
   echo $day1 >yesterday

echo 'merge all processed hourly data into a unfied file'
   startmsg
   $EXECnldas/mergeforce >> $pgmout 2>errfile                     #1b
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
    cp -rp  ${day0}${hh}.${asufix}2 $COMFOR/nldas.t${cyc}z.${day0}.force-a.grb2f${hh}

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
