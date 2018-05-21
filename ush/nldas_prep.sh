#!/bin/ksh
##############################################################
# This script performs the nldas forcing generation process
# Usage:   nldas_prep.sh start-date end-date
# History:   2013.05  Youlong Xia  Original script in C-shell
#            2013.06  Youlong Xia  Converted to Korn Shell
#            2018.02  Youlong Xia  Modified for NLDAS v2.5 
##############################################################
set -x

export alert_type=${alert_type:-NLDAS_PREP}

#using $PDY to output directory as 5-day run is put in today run

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_prep.sh start-date end-date"
  $DATA/err_exit 99
fi

sdate=$1
edate=$2

export cyc=${cyc:-12}

while [ $sdate -le $edate ]; do
  year1=`echo $sdate |cut -c1-4`
  year2=`finddate.sh $sdate d+1 |cut -c1-4`

  day1=$sdate
  day2=`finddate.sh $sdate d+1`

  export pgm=nldas_prep
  . $DATA/prep_step
  
  echo 'Here we start to run nldas forcing generation'
  # BUILD the ldas.crd file for a realtime update
  # As forcing generation algorithm needs total two complete days to ensure
  # disaggregated precipitation internal consistency, we need set start date to
  # two days prior $sdate
  cp -p $PARMnldas/force-template.crd $DATA
  cp -p $PARMnldas/TEMP $DATA
  cp -p $FIXnldas/prep/tables/KPDS_completemosaic.force.narr.tbl $DATA
  cp -p $FIXnldas/prep/tables/KPDS_completemosaic.force.narr.grib2.tbl $DATA  

  export SDAY=`finddate.sh $sdate d-2`
  export SDA=`echo $SDAY | cut -c 7-8`
  export SMO=`echo $SDAY | cut -c 5-6`
  export SYR=`echo $SDAY | cut -c 1-4`

  export EDAY=$day2
  export EDA=`echo $EDAY | cut -c 7-8`
  export EMO=`echo $EDAY | cut -c 5-6`
  export EYR=`echo $EDAY | cut -c 1-4`

  #UPDATE TIMES IN TEMPLATE CARD
  sed "s?ZDA?$SDA?" $DATA/force-template.crd  | sed "s?ZMO?$SMO?" | sed "s?ZYR?$SYR?"| sed "s?XDA?$EDA?"| sed "s?XMO?$EMO?"| sed "s?XYR?$EYR?" > $DATA/ldas.crd
  
  # COMLETE TO UPDATE TIMES IN TEMPLATE CARD
  export temout=$DATA/OUTPUT/EXP888/FOR
  if [ -s $temout ]; then
  echo $temout exists
  else
  mkdir -p $temout
  fi

  # Make directory MERGEA to store temporary merging data   
   export merdir=$DATA/MERGEA
   if [ -s $merdir ]; then
   echo $merdir exists
   else
   mkdir -p $merdir
   fi

  # Link radiation monthly bias correction ratio (derived from GOES data) 
  ln -s $FIXnldas/prep/RAD_BC $DATA/RAD_BC

  # link the parameters to generate NLDAS-2 forcing 
  ln -s  $FIXnldas/prep/BCS $DATA/BCS

  # Run the code with the updated LIS nldas card (ldas.crd)
  startmsg
  $EXECnldas/nldas_prep >> $pgmout 2>>errfile
  export err=$?; $DATA/err_chk

  sdate1=`finddate.sh $sdate d+1`

  if [ ${SENDCOM:-YES} = YES ]; then

  # Copy sdate1 ouput files to /com:

  export COM_OUT=${COM_OUT:-$COMROOTp1/nldas/${envir}}
  # copy 00Z a and b nldas forcing file
  export outdir1=$DATA/OUTPUT/EXP888/FOR/$year2/$sdate1
  export asufix=nldasforce-a.grb
  export bsufix=nldasforce-b.grb   
 
   # Make directory OUTPUT to store data
   
   if [ -s $outdir1 ]; then
   echo $outdir1 exists
   else
   mkdir -p $outdir1
   fi
	
  export COMFOR1=$COM_OUT/nldas.$PDY
  if [ -s $COMFOR1 ]; then
   echo $COMFOR1 exists
   else
   mkdir -p $COMFOR1
   fi
# copy the forcing to tomorrow run directory
   export COMFOR2=$COM_OUT/nldas.$PDYp1
   if [ -s $COMFOR2 ]; then
   echo $COMFOR2 exists
   else
   mkdir -p $COMFOR2
   fi

  mv $outdir1/${sdate1}00.${asufix}  $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grbf00
  mv $outdir1/${sdate1}00.${asufix}2 $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grb2f00
  mv $outdir1/${sdate1}00.${bsufix}2 $COMFOR1/nldas.t${cyc}z.${sdate1}.force-b.grb2f00

  cp $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grbf00 $COMFOR2/nldas.t${cyc}z.${sdate1}.force-a.grbf00
  cp $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grb2f00 $COMFOR2/nldas.t${cyc}z.${sdate1}.force-a.grb2f00
  cp  $COMFOR1/nldas.t${cyc}z.${sdate1}.force-b.grb2f00 $COMFOR2/nldas.t${cyc}z.${sdate1}.force-b.grb2f00

  if [ $SENDDBN = YES ]; then
   $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grb2f00
   $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR1/nldas.t${cyc}z.${sdate1}.force-b.grb2f00
    fi
    hh=01
    while [ $hh -le 12 ]; do
    mv $outdir1/${sdate1}${hh}.${asufix} $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grbf${hh}
    mv $outdir1/${sdate1}${hh}.${asufix}2 $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grb2f${hh}
    mv $outdir1/${sdate1}${hh}.${bsufix}2 $COMFOR1/nldas.t${cyc}z.${sdate1}.force-b.grb2f${hh}

    cp $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grbf${hh} $COMFOR2/nldas.t${cyc}z.${sdate1}.force-a.grbf${hh}
    cp $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grb2f${hh} $COMFOR2/nldas.t${cyc}z.${sdate1}.force-a.grb2f${hh}
    cp $COMFOR1/nldas.t${cyc}z.${sdate1}.force-b.grb2f${hh} $COMFOR2/nldas.t${cyc}z.${sdate1}.force-b.grb2f${hh}
      
    if [ $SENDDBN = YES ]; then
    $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR1/nldas.t${cyc}z.${sdate1}.force-a.grb2f${hh}
    $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR1/nldas.t${cyc}z.${sdate1}.force-b.grb2f${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done

  #copy sdate nldas frocing file to /com
  export outdir=$DATA/OUTPUT/EXP888/FOR/$year1/$sdate
  # Make directory OUTPUT to store data

   if [ -s $outdir ]; then
   echo $outdir exists
   else
   mkdir -p $outdir
   fi

  mv $outdir/${sdate}00.${asufix}  $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grbf00
  mv $outdir/${sdate}00.${asufix}2  $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grb2f00
  mv $outdir/${sdate}00.${bsufix}2  $COMFOR1/nldas.t${cyc}z.${sdate}.force-b.grb2f00
 
  cp $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grbf00 $COMFOR2/nldas.t${cyc}z.${sdate}.force-a.grbf00
  cp $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grb2f00 $COMFOR2/nldas.t${cyc}z.${sdate}.force-a.grb2f00
  cp $COMFOR1/nldas.t${cyc}z.${sdate}.force-b.grb2f00 $COMFOR2/nldas.t${cyc}z.${sdate}.force-b.grb2f00

  if [ $SENDDBN = YES ]; then
       $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grb2f00
       $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR1/nldas.t${cyc}z.${sdate}.force-b.grb2f00
   fi
    hh=01
    while [ $hh -le 23 ]; do
    mv $outdir/${sdate}${hh}.${asufix} $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grbf${hh}
    mv $outdir/${sdate}${hh}.${asufix}2 $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grb2f${hh}
    mv $outdir/${sdate}${hh}.${bsufix}2 $COMFOR1/nldas.t${cyc}z.${sdate}.force-b.grb2f${hh}
    cp $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grbf${hh}  $COMFOR2/nldas.t${cyc}z.${sdate}.force-a.grbf${hh}
    cp $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grb2f${hh} $COMFOR2/nldas.t${cyc}z.${sdate}.force-a.grb2f${hh} 
    cp $COMFOR1/nldas.t${cyc}z.${sdate}.force-b.grb2f${hh} $COMFOR2/nldas.t${cyc}z.${sdate}.force-b.grb2f${hh}

      if [ $SENDDBN = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR1/nldas.t${cyc}z.${sdate}.force-a.grb2f${hh}
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COMFOR1/nldas.t${cyc}z.${sdate}.force-b.grb2f${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done
    fi
 
  sdate=`finddate.sh $sdate1 d+1`
done
