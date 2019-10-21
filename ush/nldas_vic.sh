#!/bin/ksh
############################################################################
# This script performs the VIC land-surface analysis
# The Variable Infitration Capacity (VIC) was jointly developed by the 
# University of Washington and Princeton University as a hydrological model. 
# Usage:   nldas_vic.sh start-date end-date
# History:   2013.05  Justin Sheffield  Original script 
#            2013.06  Youlong Xia Converted to Korn Shell to fit 
#                     NCEP operational requirments   
#            2014.01  Fixed a typo "ENDEAR" to "ENDYEAR" in the substitution 
#                       for the global_retro_v4.0.3.txt file
#            2019.03  Use modular for VIC run    
############################################################################
#       Control script for the LDAS version of the VIC model.
#       Running daily simulations on Linux X86 machines at NCEP.
#       The script:
#               1. search for and extract the daily forcing data
#               2. update the VIC state initialisation file
#               3. update the VIC global input file with today's info
#               4. run the VIC pre-processor
#               5. run the VIC model
#               6. run the VIC post-processor
#               7. back-up the log, state and output data
############################################################################
set -x

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_vic.sh start-date end-date"
  err_exit 99
fi

sdate=$1
edate=$2

# log file names

export LOG_FILE=$DATA/VIC_RUN/log
export SUMMARY_LOG_FILE=$DATA/VIC_RUN/summary.log
#export VIC_LOG_FILE=$DATA/VIC_RUN/vic.log
export QC_LOG_FILE=$DATA/VIC_RUN/qc.log

# some global vars
TRUE=0
FALSE=1
FLAG_OK=${FALSE}
FORCE_OK=${FALSE}
RESTART_OK=${FALSE}
EXIT_STATUS=0
day=0
export alert_type=${alert_type:-NLDAS_VIC}

# forcing and vic output file extensions
  FORCE_FILE_EXT=FORCING.GRB
  VIC_FILE_EXT=VIC.grb

while [ $sdate -le $edate ]; do
  year1=`echo $sdate |cut -c1-4`
  year0=`finddate.sh $sdate d-1 |cut -c1-4`
  year2=`finddate.sh $sdate d+1 |cut -c1-4`

  day1=$sdate
  day0=`finddate.sh $sdate d-1`
  day2=`finddate.sh $sdate d+1`

  export pgm=nldas_prep
  . prep_step

  # copy initial conditions for vic model
  export RESDIR=${RESDIR:-$COM_IN}
  export viclog=$DATA/VIC_RUN
  if test ! -d ${viclog}; then mkdir -p ${viclog}; fi
  cp $RESDIR/nldas.$sdate/vic.t${cyc}z.VICrst $viclog/state.in
  
  # link fix fileds to $DATA directory
  ln -s $FIXnldas/vic $DATA/VIC_PARAM

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
 
  # copy last hour grib forcing file for pre-processor
  export narrd0=$DATA/NLDASNARR/$year0/$day0
  mkdir -p $narrd0
  
  export afile23=$COM_IN/nldas.$day0/nldas.t${cyc}z.force-a.grbf23
  cp $afile23 $narrd0/${day0}23.${asufix}
 
 # copy the first hour grib forcing file  

  sdate1=`finddate.sh $sdate d+1`
  
  ## comlete nldas foricng setup 
  echo "Here we build global parameter text file"
  # BUILD the vic global parameter file for a realtime update
  cp $PARMnldas/global_retro_v4.0.3.txt $DATA/.
  export GLOBAL_FILE=global_retro_v4.0.3.txt

  # extract the date elements for the start day
  export SDAY=$day1
  export HOUR=00
  export DAY=`echo $SDAY | cut -c 7-8`
  export MONTH=`echo $SDAY | cut -c 5-6`
  export YEAR=`echo $SDAY | cut -c 1-4`

  #UPDATE TIMES and DATES IN GLOBAL FILE FOR START DATE
  echo "start date:  $YEAR/$MONTH/$DAY:$HOUR"

  # update the start dates in the VIC global file
    YEAR="STARTYEAR    ${YEAR}          # year ncep data starts"
    MONTH="STARTMONTH   ${MONTH}        # month ncep data starts"
    DAY="STARTDAY     ${DAY}            # day ncep data starts"
    HOUR="STARTHOUR    ${HOUR}          # hour ncep data starts"
    sed -e '/^STARTYEAR/ s/.*/'"$YEAR"'/' \
        -e '/^STARTMONTH/ s/.*/'"$MONTH"'/' \
        -e '/^STARTDAY/ s/.*/'"$DAY"'/' \
        -e '/^STARTHOUR/ s/.*/'"$HOUR"'/' ${GLOBAL_FILE} > temp1.txt
 
  # extract the date elements for the end day
    export EDAY=$day1
    export HOUR=23
    export DAY=`echo $EDAY | cut -c 7-8`
    export MONTH=`echo $EDAY | cut -c 5-6`
    export YEAR=`echo $EDAY | cut -c 1-4`

  #UPDATE TIMES and DATES IN GLOBAL FILE END DATE
  echo "end date:  $YEAR/$MONTH/$DAY:$HOUR"

  # update the end dates in the VIC global file
    YEAR="ENDYEAR    ${YEAR}          # year ncep data ends"
    MONTH="ENDMONTH   ${MONTH}        # month ncep data ends"
    DAY="ENDDAY     ${DAY}            # day ncep data endss"
    HOUR="ENDHOUR    ${HOUR}          # hour ncep data ends"
    sed -e '/^ENDYEAR/ s/.*/'"$YEAR"'/' \
        -e '/^ENDMONTH/ s/.*/'"$MONTH"'/' \
        -e '/^ENDDAY/ s/.*/'"$DAY"'/' \
        -e '/^ENDHOUR/ s/.*/'"$HOUR"'/' temp1.txt > temp2.txt

  # update the date elements for model output and restart storage
    export SDAY=$day1
    export HOUR=23
    export DAY=`echo $SDAY | cut -c 7-8`
    export MONTH=`echo $SDAY | cut -c 5-6`
    export YEAR=`echo $SDAY | cut -c 1-4`

    echo "restart dump date: $YEAR/$MONTH/$DAY:$HOUR"

    YEAR="STATEYEAR    ${YEAR}          # year to write out state file"
    MONTH="STATEMONTH   ${MONTH}        # month to write out state file"
    DAY="STATEDAY     ${DAY}            # day to write out state file"
    HOUR="STATEHOUR    ${HOUR}          # hour to write out state file"
    sed -e '/^STATEYEAR/ s/.*/'"$YEAR"'/' \
        -e '/^STATEMONTH/ s/.*/'"$MONTH"'/' \
        -e '/^STATEDAY/ s/.*/'"$DAY"'/' \
        -e '/^STATEHOUR/ s/.*/'"$HOUR"'/' temp2.txt > temp3.txt
    cp -p  temp3.txt ${GLOBAL_FILE}
    rm -rf temp*.txt
   

   # make meteorological forcing directory after vic pre-processor
   export metdir=$DATA/VIC_MET
   mkdir -p $metdir

  # make wtiting directory for VIC after vic pre-processor
  export vicout=$DATA/VIC_OUT
  mkdir -p $vicout

  ########
  # VIC pre-processor
  ########
  echo "VIC pre-processor" 
  # run the VIC post-processor
  echo "Pre-processor is running..." 
  export pgm=nldas_vic_prep
  . prep_step
  startmsg
  $EXECnldas/nldas_vic_prep ${GLOBAL_FILE}   
  export err=$?; err_chk

  cp -p ${QC_LOG_FILE} ${viclog}/qc.log.${day1} >> ${LOG_FILE} 2>&1

  ################
  # RUN VIC model
  ###############
  echo "VIC model" 
  echo "------------------" 
  echo "VIC model is running..." 
  export pgm=nldas_vic_nldas
  . prep_step
  startmsg
  $EXECnldas/nldas_vic_ldas  -g ${GLOBAL_FILE}
  export err=$?; err_chk

  # make vic grib output directory
  export OUTPUT_GRIB_DIR=$DATA/VIC_OUT_grib_retro_v4.0.3
  mkdir -p $OUTPUT_GRIB_DIR
  
  #####################
  # VIC post-processor
  #####################

  echo "VIC post-processor" 
  echo "------------------" 

  # run the VIC post-processor
  # remove tabs from the global file for the fortran preprocessor


   sed 's/     / /g' ${GLOBAL_FILE} > ${GLOBAL_FILE}.tmp
  
   export pgm=nldas_vic_post
   . prep_step
   startmsg
   $EXECnldas/nldas_vic_post ${GLOBAL_FILE}.tmp
   export err=$?; err_chk

   # Archive model initials
   export COMREST=${COMREST:-$COM_OUT}
   mv  $viclog/state.out $COMREST/nldas.$sdate1/vic.t${cyc}z.VICrst

   # Copy sdate vic model ouput files to /com:
   # copy VIC output file form 00Z-23z on day1
   hh=00
   while [ $hh -le 23 ]; do
     mv $OUTPUT_GRIB_DIR/${sdate}${hh}.VIC.grb $COM_OUT/nldas.$sdate/vic.t${cyc}z.grbf${hh}
     if [ $SENDDBN = YES ]; then
        $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COM_OUT/nldas.$sdate/vic.t${cyc}z.grbf${hh}
     fi
     let "hh=hh+1"
     if [ $hh -lt 10 ]; then hh=0$hh; fi
  done
   
  sdate=`finddate.sh $sdate d+1`
done
