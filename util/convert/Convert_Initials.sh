#!/bin/ksh
set -x

#### prepare today test run for nldas2.5 
pday=`date -u +%Y%m%d`

#### if prepare today run, comment out pday=`........ $pday d+1`

#### preapre tomorrow test run for nldas v2.5 
pday=`/nwprod/util/ush/finddate.sh $pday d+1`

sdate=`/nwprod/util/ush/finddate.sh $pday d-5`
edate=`/nwprod/util/ush/finddate.sh $pday d-4`

export indir=/com/nldas/prod/nldas.${edate}
export outd=/land/noscrub/Youlong.Xia/com/nldas.v2.5.0/dev/nldas.${pday}

mkdir -p $outd

cp -rp $indir/mosaic.t12z.MOSrst $outd/mosaic.t12z.${edate}.MOSrst
cp -rp $indir/mosaic.t12z.internal_runoff.bin $outd/mosaic.t12z.${edate}.internal_runoff.bin
cp -rp $indir/mosaic.t12z.transport_runoff.bin $outd/mosaic.t12z.${edate}.transport_runoff.bin

cp -rp $indir/noah.t12z.INIT_CH $outd/noah.t12z.${edate}.INIT_CH
cp -rp $indir/noah.t12z.INIT_CM $outd/noah.t12z.${edate}.INIT_CM
cp -rp $indir/noah.t12z.INIT_CMC $outd/noah.t12z.${edate}.INIT_CMC
cp -rp $indir/noah.t12z.INIT_LSTSNW $outd/noah.t12z.${edate}.INIT_LSTSNW
cp -rp $indir/noah.t12z.INIT_SH2O $outd/noah.t12z.${edate}.INIT_SH2O
cp -rp $indir/noah.t12z.INIT_SMC $outd/noah.t12z.${edate}.INIT_SMC
cp -rp $indir/noah.t12z.INIT_SNEQV $outd/noah.t12z.${edate}.INIT_SNEQV
cp -rp $indir/noah.t12z.INIT_SNOWH $outd/noah.t12z.${edate}.INIT_SNOWH
cp -rp $indir/noah.t12z.INIT_STC $outd/noah.t12z.${edate}.INIT_STC
cp -rp $indir/noah.t12z.INIT_T1 $outd/noah.t12z.${edate}.INIT_T1

cp -rp $indir/noah.t12z.internal_runoff.bin $outd/noah.t12z.${edate}.internal_runoff.bin
cp -rp $indir/noah.t12z.transport_runoff.bin $outd/noah.t12z.${edate}.transport_runoff.bin

cp -rp $indir/sac.t12z.INIT_SMC $outd/sac.t12z.${edate}.INIT_SMC
cp -rp $indir/sac.t12z.INIT_SNOWCO $outd/sac.t12z.${edate}.INIT_SNOWCO

cp -rp $indir/sac.t12z.internal_runoff.bin $outd/sac.t12z.${edate}.internal_runoff.bin
cp -rp $indir/sac.t12z.transport_runoff.bin $outd/sac.t12z.${edate}.transport_runoff.bin

cp -rp $indir/vic.t12z.VICrst $outd/vic.t12z.${edate}.VICrst
cp -rp $indir/vic.t12z.internal_runoff.bin $outd/vic.t12z.${edate}.internal_runoff.bin
cp -rp $indir/vic.t12z.transport_runoff.bin $outd/vic.t12z.${edate}.transport_runoff.bin
