#!/bin/ksh
##############################################################
# This script performs the SAC land-surface analysis
# Usage:   nldas_sac.sh start-date end-date
# History:   2013.05  Youlong Xia:  Original script in C-shell
#            2013.06  Youlong Xia and Yuqiu Zhu:  Converted to Korn 
#                     Shell script for product implementation
##############################################################
set -x

export alert_type=${alert_type:-NLDAS_SAC}

cd $DATA

if [ $# -lt 2 ]; then
  echo "Usage: nldas_sac.sh start-date end-date"
  err_exit 99
fi

sdate=$1
edate=$2

while [ $sdate -le $edate ]; do
  julday=`$HOMEnldas/util/ush/julian.sh $sdate`
  year=`echo $sdate |cut -c1-4`
  year1=`finddate.sh $sdate d-1 |cut -c1-4`

  today=`finddate.sh $sdate d+1`
  year2=`echo $today |cut -c1-4`

  echo $julday >julday
  echo $sdate >yesterday
  echo $today >today 

  export pgm=nldas_sac
  . prep_step

  #Input files:
  echo 'linking Noah land surface characteristics'
  export FORT65=$FIXnldas/sac/UMDunifiedmask19990327.bin
  export FORT51=$FIXnldas/noah/noah_soil_layers19990327.bin
  export FORT53=$FIXnldas/noah/noah_soiltype19990808.bin
  export FORT54=$FIXnldas/sac/UMDpveg19990811.bin
  export FORT55=$FIXnldas/noah/noah_slopetype19990325.bin
  export FORT57=$FIXnldas/noah/noah_soil_lower_T19990325.bin
  export FORT58=$FIXnldas/noah/noah_albedo19990420.bin
  export FORT59=$FIXnldas/noah/noah_gfrac19990325.bin
  export FORT60=$FIXnldas/sac/precipmask.2deg.bin
  export FORT61=$FIXnldas/sac/PE_LDAS.bin
  export FORT62=$FIXnldas/sac/LDAS_ELEV.bin 
  export FORT63=$FIXnldas/sac/PE_RATIO.bin
  export FORT64=$FIXnldas/sac/PEADJ_LDAS.bin


  echo 'linking initial conditions for SAC model'
  export RESDIR=${RESDIR:-$COM_IN}
  export FORT72=$RESDIR/nldas.${sdate}/sac.t${cyc}z.INIT_SMC
  export FORT74=$RESDIR/nldas.${sdate}/sac.t${cyc}z.INIT_SNOWCO

  echo 'linking sac parameters'
  # snow model parameters
  export FORT71=$FIXnldas/sac/SNOW_PAR

  export FORT81=$FIXnldas/sac/nldas_sac_adimp.bin
  export FORT82=$FIXnldas/sac/nldas_sac_lzfpm.bin
  export FORT83=$FIXnldas/sac/nldas_sac_lzfsm.bin
  export FORT84=$FIXnldas/sac/nldas_sac_lzpk.bin
  export FORT85=$FIXnldas/sac/nldas_sac_lzsk.bin
  export FORT86=$FIXnldas/sac/nldas_sac_lztwm.bin
  export FORT87=$FIXnldas/sac/nldas_sac_pctim.bin
  export FORT88=$FIXnldas/sac/nldas_sac_pfree.bin
  export FORT89=$FIXnldas/sac/nldas_sac_rexp.bin
  export FORT90=$FIXnldas/sac/nldas_sac_riva.bin
  export FORT91=$FIXnldas/sac/nldas_sac_rserv.bin
  export FORT92=$FIXnldas/sac/nldas_sac_side.bin
  export FORT93=$FIXnldas/sac/nldas_sac_uzfwm.bin
  export FORT94=$FIXnldas/sac/nldas_sac_uzk.bin
  export FORT95=$FIXnldas/sac/nldas_sac_uztwm.bin
  export FORT96=$FIXnldas/sac/nldas_sac_zperc.bin

  echo 'linking nldas prep forcing files'
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

  sdate1=`finddate.sh $sdate d+1`
  export FORT25=$COMFORCE/nldas.$sdate1/nldas.t${cyc}z.force-a.grbf00

  
  echo 'linking Noah output (PE) as SAC forcing'

  export COMNOAH=${COMNOAH:-$COM_IN}
  export FORT26=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf00
  export FORT27=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf01
  export FORT28=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf02
  export FORT29=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf03
  export FORT30=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf04
  export FORT31=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf05
  export FORT32=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf06
  export FORT33=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf07
  export FORT34=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf08
  export FORT35=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf09
  export FORT36=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf10
  export FORT37=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf11
  export FORT38=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf12
  export FORT39=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf13
  export FORT40=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf14
  export FORT41=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf15
  export FORT42=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf16
  export FORT43=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf17
  export FORT44=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf18
  export FORT45=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf19
  export FORT46=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf20
  export FORT47=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf21
  export FORT48=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf22
  export FORT49=$COMNOAH/nldas.$sdate/noah.t${cyc}z.grbf23

  export FORT50=$COMNOAH/nldas.$sdate1/noah.t${cyc}z.grbf00

  # Output SAC initials files:
  export FORT77=${sdate1}00.INIT_SMC
  export FORT78=${sdate1}00.INIT_SNOWCO
  
  # Build the SAC control file
  cat <<EOF >cntl_sac
#Time step in seconds
3600.0
#Forcing height in meters
6.0
#Snow albedo
0.75
#Ice flag
0
#NOAH ETP flag: 1-using Noah PE, 2-using PE climatology from OBS 
1
EOF

  # Copy the KPDS table:
  cp $FIXnldas/sac/KPDS.tbl .
  
  echo 'running SAC model'
  startmsg
  $EXECnldas/nldas_sac_ldas >> $pgmout 2>errfile
  export err=$?; err_chk
  
  chmod 644 *grb
  if [ ${SENDCOM:-YES} = YES ]; then
    # Copy restart files to /com:
    export COMREST=${COMREST:-$COM_OUT}
    mv ${sdate1}00.INIT_SMC $COMREST/nldas.$sdate1/sac.t${cyc}z.INIT_SMC
    mv ${sdate1}00.INIT_SNOWCO $COMREST/nldas.$sdate1/sac.t${cyc}z.INIT_SNOWCO
   
    # Copy grib files to /com:
    mv ${sdate1}00.SAC.grb  $COM_OUT/nldas.$sdate1/sac.t${cyc}z.grbf00
    if [ $SENDDBN = YES ]; then
       $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COM_OUT/nldas.$sdate1/sac.t${cyc}z.grbf00
    fi
    hh=01
    while [ $hh -le 23 ]; do
      mv ${sdate}${hh}.SAC.grb $COM_OUT/nldas.$sdate/sac.t${cyc}z.grbf${hh}
      if [ $SENDDBN = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL ${alert_type} $job $COM_OUT/nldas.$sdate/sac.t${cyc}z.grbf${hh}
      fi
      let "hh=hh+1"
      if [ $hh -lt 10 ]; then hh=0$hh; fi
    done
  fi
 
  sdate=`finddate.sh $sdate d+1`
done
