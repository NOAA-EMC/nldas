sorc=convert_cpcops_cpcus_cpcglobe.f
pgrm=cpc_precip_convert

lib=/nwprod/lib/libip_4.a /nwprod/lib/libsp_4.a /nwprod/lib/libbacio_4.a /nwprod/lib/libw3nco_4.a -openmp 
 
FC=ifort
FFLAGS=-extend_source -convert big_endian

${pgrm}:	${sorc}
	${FC} ${FFLAGS} -o ${pgrm} ${sorc} ${lib}

clean:
	rm ${pgrm}
