onvert_forcing_model.sh
set -x

#### prepare today test run for nldas2.5
pday=`date -u +%Y%m%d`

#### if prepare for today run, 
#### comment out pday=`/nwprod/util/ush/finddate.sh $pday d+1`

#### preapre tomorrow test run for nldas v2.5
pday=`/nwprod/util/ush/finddate.sh $pday d+1`

sdate=`/nwprod/util/ush/finddate.sh $pday d-5`
edate=`/nwprod/util/ush/finddate.sh $pday d-4`

export indir=/com/nldas/prod
export outd=/land/noscrub/Youlong.Xia/com/nldas.v2.5.0/dev/nldas.${pday}/

mkdir -p $outd
cd $outd

cp -rp $indir/nldas.${sdate}/nldas.t12z.force-a.grbf* .

for file in `ls nldas.t12z.force-a.grbf*`
do
  hh=`echo $file |awk -F"." '{print $4}'|cut -c5-6`
  type=`echo $file |awk -F"." '{print $3}'|cut -c1-7`
  mv $file nldas.t12z.${sdate}.${type}.grbf${hh}
done

cp -rp $indir/nldas.${edate}/*grb* .

for file in `ls nldas.t12z.force-a.grbf*`
do
  hh=`echo $file |awk -F"." '{print $4}'|cut -c5-6`
  type=`echo $file |awk -F"." '{print $3}'|cut -c1-7`
  mv $file nldas.t12z.${edate}.${type}.grbf${hh}
done

for file in `ls nldas.t12z.force*grb2f*`
do
  hh=`echo $file |awk -F"." '{print $4}'|cut -c6-7`
  type=`echo $file |awk -F"." '{print $3}'|cut -c1-7`
  mv $file nldas.t12z.${edate}.${type}.grb2f${hh}
done

for file in `ls mosaic.t12z.grbf*`
do
  type=`echo $file |awk -F"." '{print $1}'|cut -c1-6`
  hh=`echo $file |awk -F"." '{print $3}'|cut -c5-6`
  mv $file nldas.t12z.${edate}.${type}.grb2f${hh}
done

mv noah.t12z.grbf00 nldas.t12z.${edate}.noah.grb2f00
mv sac.t12z.grbf00 nldas.t12z.${edate}.sac.grb2f00


