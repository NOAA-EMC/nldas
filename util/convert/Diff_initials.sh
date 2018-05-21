#!/bin/ksh
set -x

export day0=20180402
export day1=20180406
export dir1=/land/noscrub/Youlong.Xia/com/nldas.v2.5.2/dev/nldas.${day1}
export dir2=/com/nldas/prod/nldas.${day0}

diff $dir1/mosaic.t12z.$day0.MOSrst $dir2/mosaic.t12z.MOSrst
diff $dir1/mosaic.t12z.$day0.internal_runoff.bin $dir2/mosaic.t12z.internal_runoff.bin
diff $dir1/mosaic.t12z.$day0.transport_runoff.bin $dir2/mosaic.t12z.transport_runoff.bin

diff $dir2/noah.t12z.INIT_CH $dir1/noah.t12z.$day0.INIT_CH
diff $dir2/noah.t12z.INIT_CM $dir1/noah.t12z.$day0.INIT_CM
diff $dir2/noah.t12z.INIT_CMC $dir1/noah.t12z.$day0.INIT_CMC
diff $dir2/noah.t12z.INIT_LSTSNW $dir1/noah.t12z.$day0.INIT_LSTSNW
diff $dir2/noah.t12z.INIT_SH2O $dir1/noah.t12z.$day0.INIT_SH2O
diff $dir2/noah.t12z.INIT_SMC $dir1/noah.t12z.$day0.INIT_SMC
diff $dir2/noah.t12z.INIT_SNEQV $dir1/noah.t12z.$day0.INIT_SNEQV
diff $dir2/noah.t12z.INIT_SNOWH $dir1/noah.t12z.$day0.INIT_SNOWH
diff $dir2/noah.t12z.INIT_STC $dir1/noah.t12z.$day0.INIT_STC
diff $dir2/noah.t12z.INIT_T1 $dir1/noah.t12z.$day0.INIT_T1

diff $dir2/noah.t12z.internal_runoff.bin $dir1/noah.t12z.$day0.internal_runoff.bin
diff $dir2/noah.t12z.transport_runoff.bin $dir1/noah.t12z.$day0.transport_runoff.bin

diff $dir2/sac.t12z.INIT_SMC $dir1/sac.t12z.$day0.INIT_SMC
diff $dir2/sac.t12z.INIT_SNOWCO $dir1/sac.t12z.$day0.INIT_SNOWCO

diff $dir2/sac.t12z.internal_runoff.bin $dir1/sac.t12z.$day0.internal_runoff.bin
diff $dir2/sac.t12z.transport_runoff.bin $dir1/sac.t12z.$day0.transport_runoff.bin

diff $dir2/vic.t12z.VICrst $dir1/vic.t12z.$day0.VICrst
diff $dir2/vic.t12z.internal_runoff.bin $dir1/vic.t12z.$day0.internal_runoff.bin
diff $dir2/vic.t12z.transport_runoff.bin $dir1/vic.t12z.$day0.transport_runoff.bin
